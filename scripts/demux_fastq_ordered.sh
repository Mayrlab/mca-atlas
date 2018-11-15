#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo -e "Incorrect number of parameters! Usage:\n    demux_fastq.sh <fastq> <out_dir> <prefix> <threads> <buffer_size>" >&2
    exit 1
fi

fastq="$1"
outDir="$2"
prefix="$3"
threads="$4"
bufferSize="$5"

mkdir -p $outDir || exit 1

echo "Input file path:" $fastq
echo "Output directory:" $outDir
echo "Prefix:" $prefix
echo "Threads:" $threads
echo "Buffer Size:" $bufferSize

if [[ $fastq =~ \.gz$ ]]; then
    zcat $fastq
else
    cat $fastq
fi | awk -f scripts/collapse_fastq.awk | sort -k1,1 --parallel $threads -S $bufferSize |
    awk -v prefix="$prefix" -v outdir="$outDir" -f scripts/demux_mca_collapsed.awk

## Generate batch file
ls $outDir/*.fastq | awk -v prefix="$prefix" -v outdir="$outDir" '
{
   match($0, /(.*\/)([^\/]+).fastq/, cell); 
   print cell[2]"\t"cell[1]cell[2]".umi\t"$0 > outdir"/"prefix".dat"
}'

echo "Demultiplexing complete."
