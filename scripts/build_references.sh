#!/bin/bash

if [ "$#" -ne 4]; then
    echo "Missing at least one of the four required arguments: <path to the input fasta file> <output directory> <prefix for generated files> <path to ViFi scripts directory>"
    exit 1
fi

INPUT_FASTA=$1
OUTPUT_DIR=$2
PREFIX=$3
SRC_DIR=$4

mkdir -p $OUTPUT_DIR
OUTPUT_DIR=`realpath $OUTPUT_DIR`
INPUT_DIR=`dirname $INPUT_FASTA`
if [ "$INPUT_DIR" == "." ]
then
  INPUT_DIR=`pwd`
fi 
INPUT_NAME=`basename $INPUT_FASTA`

#Build alignment/tree, use 4GB for alignment, via Docker
docker run -v $INPUT_DIR/:/data/ -v $OUTPUT_DIR/:/output/ smirarab/pasta run_pasta.py  --max-mem-mb=120000 -d dna -j $PREFIX -i $INPUT_NAME  -o /output/ --merger=muscle

cd $OUTPUT_DIR

PASTA_LOG="$OUTPUT_DIR/${PREFIX}.out.txt"
echo $PASTA_LOG
OUT_ALN=`grep "Writing resulting alignment" $PASTA_LOG | awk '{print $7}'`
echo "OUT_ALN: $OUT_ALN"
OUT_ALN=`basename $OUT_ALN`
OUT_TREE=`grep "Writing resulting tree" $PASTA_LOG | awk '{print $7}'`
echo "OUT_TREE: $OUT_TREE"
OUT_TREE=`basename $OUT_TREE`

#Decompose alignment/tree into HMMs via Docker
docker run -v $OUTPUT_DIR/:/output/ -v $SRC_DIR/build_hmms.py:/home/scripts/build_hmms.py docker.io/namphuon/vifi python "scripts/build_hmms.py" --tree_file /output/$OUT_TREE --alignment_file /output/$OUT_ALN --prefix $PREFIX --output_dir /output/ --keep_alignment

#Build HMM list
ls $OUTPUT_DIR/*.hmmbuild > $OUTPUT_DIR/hmm_list.txt
