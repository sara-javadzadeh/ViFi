#!/bin/bash

VIFI_DIR=`pwd`

#Get data repos
if [ ! -d "data_repo" ]; then
    echo "Downloading the data_repo"
    wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/master/gdown.pl
    perl gdown.pl "https://drive.google.com/file/d/1il10KUxJ5Q5JvR5pHJB4GUMBlBPgTjrj/view?usp=sharing" data_repo.tar.gz
    tar -zxvf data_repo.tar.gz
    rm data_repo.tar.gz
fi
if [ ! -d "viral_data" ]; then
    echo "Uncompressing the HMM models"
    tar -xzvf viral_data.tar.gz
fi

#Set up environmental variables
echo "Set environmental variables"
echo export VIFI_DIR=$VIFI_DIR >> ~/.bashrc
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
echo export REFERENCE_REPO=$PWD/viral_data >> ~/.bashrc

VIFI_DIR=$VIFI_DIR
AA_DATA_REPO=$PWD/data_repo
REFERENCE_REPO=$PWD/viral_data

source ~/.bashrc

#Pull the Docker file
echo "Getting the dockerized version of ViFi"
docker pull docker.io/namphuon/vifi

#Set up reference for alignment
HUMAN_REF="GRCh38"
HUMAN_REF_FILE_NAME="hg38full.fa"
#for virus in "hpv_655" "hbv_2012"; do
for virus in "hpv" "hbv" "hcv" "ebv"; do
    HUMAN_VIRAL_REF="grch38_${virus}.fas"
    echo "Building the ${HUMAN_REF}+${virus} reference"
    cat $AA_DATA_REPO//${HUMAN_REF}/${HUMAN_REF_FILE_NAME} $REFERENCE_REPO/${virus}/${virus}.unaligned.fas > $REFERENCE_REPO/${virus}/${HUMAN_VIRAL_REF}
    docker run -v $REFERENCE_REPO/${virus}/:/home/${virus}/ docker.io/namphuon/vifi bwa index /home/${virus}/${HUMAN_VIRAL_REF}

    #Build reduced list of HMMs for testing
    echo "building the list of hmms for testing in $VIFI_DIR"
    ls $VIFI_DIR/viral_data/${virus}/hmms/*.hmmbuild > $VIFI_DIR/viral_data/${virus}/hmms/hmms.txt
    ls $VIFI_DIR/viral_data/${virus}/hmms/*.[0-9].hmmbuild > $VIFI_DIR/viral_data/${virus}/hmms/partial_hmms.txt

done

#Run ViFi under docker mode on HPV test dataset on reduced HMM list set
virus="hpv"
echo "Running test for ViFi"
python $VIFI_DIR/scripts/run_vifi.py --cpus 2 --hmm_list $VIFI_DIR/viral_data/${virus}/hmms/partial_hmms.txt -f $VIFI_DIR/test/data/test_R1.fq.gz -r $VIFI_DIR/test/data/test_R2.fq.gz -o $VIFI_DIR/tmp/docker/ --docker --viral_reference_dir $REFERENCE_REPO --hg_data_dir $AA_DATA_REPO --vifi_dir $VIFI_DIR
