#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo -e "Incorrect number of parameters! Usage:\n    demux_fastq.sh <fastq> <out_dir> <prefix>" >&2
    exit 1
fi

fastq="$1"
outDir="$2"
prefix="$3"

mkdir -p $outDir

echo "Input file path:" $fastq
echo "Output directory:" $outDir
echo "Prefix:" $prefix

if [[ $fastq =~ \.gz$ ]]; then
    zcat $fastq | awk -v prefix="$prefix" -v outdir="$outDir" -f scripts/demux_mca.awk
else
    cat $fastq | awk -v prefix="$prefix" -v outdir="$outDir" -f scripts/demux_mca.awk
fi

echo "Demultiplexing complete."
