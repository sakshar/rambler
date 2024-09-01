#!/bin/bash

# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# minimap2, seqkit (channel: bioconda), biopython (channel: conda-forge)
# commands to install them:
# conda install -c bioconda minimap2
# conda install -c bioconda seqkit
# conda install -c conda-forge biopython


# Parse command-line options using getopts
while getopts "g:n:b:e:o:s:" opt; do
  case $opt in
    g)
      REF_FILE="$OPTARG"
      ;;
    n)
      CONTIG_ID="$OPTARG"
      ;;
    b)
	  start="$OPTARG"
	  ;;
	e)
	  end="$OPTARG"
	  ;;
	o)
	  OUT_DIR="$OPTARG"
	  ;;
	s)
	  REF_SIZE="$OPTARG"
	  ;;
    \?)
      echo "Usage: $0 -g draft_genome_file -n contig_id -b start_position -e end_position -o output_directory -s contig_length"
      exit 1
      ;;
  esac
done

# Shift the positional parameters to exclude the parsed options
shift "$((OPTIND - 1))"

# Check if the required options were provided
if [ -z "$REF_FILE" ]; then
  echo "Error: -g option is required."
  exit 1
fi
if [ -z "$CONTIG_ID" ]; then
  echo "Error: -n option is required."
  exit 1
fi
if [ -z "$start" ]; then
  echo "Error: -b option is required."
  exit 1
fi
if [ -z "$end" ]; then
  echo "Error: -e option is required."
  exit 1
fi
if [ -z "$OUT_DIR" ]; then
  echo "Error: -o option is required."
  exit 1
fi
if [ -z "$REF_SIZE" ]; then
  echo "Error: -s option is required."
  exit 1
fi

ASSEMBLY=$OUT_DIR/assembly
FINAL=$OUT_DIR/final

mkdir -p $FINAL

if [ "$start" -eq "-1" ]; then
    seqkit subseq -r $end:-1 --chr $CONTIG_ID $REF_FILE > $FINAL/end.fasta
    cat $ASSEMBLY/rambler.fasta $FINAL/end.fasta > $FINAL/combined.fasta
elif [ "$end" -eq "-1" ]; then
    seqkit subseq -r 1:$start --chr $CONTIG_ID $REF_FILE > $FINAL/start.fasta
    cat $FINAL/start.fasta $ASSEMBLY/rambler.fasta > $FINAL/combined.fasta
else
    seqkit subseq -r 1:$start --chr $CONTIG_ID $REF_FILE > $FINAL/start.fasta
    seqkit subseq -r $end:-1 --chr $CONTIG_ID $REF_FILE > $FINAL/end.fasta
    cat $FINAL/start.fasta $ASSEMBLY/rambler.fasta $FINAL/end.fasta > $FINAL/combined.fasta
fi

minimap2 -x ava-pb $FINAL/combined.fasta $FINAL/combined.fasta > $FINAL/merge_overlaps.paf

python3 ./from_rambler_to_merged.py $OUT_DIR $REF_SIZE
