#!/bin/bash
INPUT_FASTA=$1
PREFIX=$2
HMMS=$3

OUTPUT_DIR=$REFERENCE_REPO/$PREFIX
mkdir -p $OUTPUT_DIR
OUTPUT_DIR=`realpath $OUTPUT_DIR`

#Create hg19+virus fasta file
HUMAN_REF="grch38"
HUMAN_REF_FILE_NAME="hg38full.fa"
HUMAN_REF_DIR="GRCh38"
cat $AA_DATA_REPO//${HUMAN_REF_DIR}/${HUMAN_REF_FILE_NAME} $REFERENCE_REPO/${virus}/${virus}.unaligned.fas > $REFERENCE_REPO/${virus}/${HUMAN_REF}_${virus}.fas

#Build index
docker run -v $REFERENCE_REPO/$PREFIX/:/home/$PREFIX/ docker.io/namphuon/vifi bwa index /home/$PREFIX/${HUMAN_REF}_${PREFIX}.fas

