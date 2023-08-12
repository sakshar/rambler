#!/bin/bash

# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# hifiasm, minimap2 (channel: bioconda), biopython (channel: conda-forge)
# commands to install them:
# conda install -c bioconda hifiasm
# conda install -c bioconda minimap2
# conda install -c conda-forge biopython


# Default values
k_len=21
tolerance=15
threshold=15

# Parse command-line options using getopts
while getopts "r:u:o:s:k:l:h:" opt; do
  case $opt in
    r)
      READ_FILE="$OPTARG"
      ;;
    u)
      KMER_FILE="$OPTARG"
      ;;
	o)
	  OUT_DIR="$OPTARG"
	  ;;
	s)
	  REF_SIZE="$OPTARG"
	  ;;
	k)
	  k_len="$OPTARG"
	  ;;
	l)
	  tolerance="$OPTARG"
	  ;;
	h)
	  threshold="$OPTARG"
	  ;;
    \?)
      echo "Usage: $0 -r read_file -u unikmer_file -o output_directory -s assembly_length [-k kmer_length] [-l tolerance] [-h threshold]"
      exit 1
      ;;
  esac
done

# Shift the positional parameters to exclude the parsed options
shift "$((OPTIND - 1))"

# Check if the required options were provided
if [ -z "$READ_FILE" ]; then
  echo "Error: -r option is required."
  exit 1
fi
if [ -z "$KMER_FILE" ]; then
  echo "Error: -u option is required."
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


INTERMEDIATES=$OUT_DIR/intermediates
CLUSTERS=$OUT_DIR/clusters
ASSEMBLY=$OUT_DIR/assembly

mkdir -p $OUT_DIR
mkdir -p $INTERMEDIATES
mkdir -p $CLUSTERS
mkdir -p $ASSEMBLY

# command for converting the read ids to the expected format
awk '/^>/{print ">S1_" ++i; next}{print}' < $READ_FILE > $INTERMEDIATES/reads.fasta

python3 ./from_read_unikmer_to_clusters.py $READ_FILE $KMER_FILE $OUT_DIR $REF_SIZE $k_len $tolerance $threshold

rm $INTERMEDIATES/reads.fasta

j=0
for a in $CLUSTERS/*.fasta;
do
	hifiasm -o $ASSEMBLY/hifiasm/"$j" -t 20 $a
	awk '/^S/{print ">"$2;print $3}' $ASSEMBLY/hifiasm/"$j".bp.p_ctg.gfa > $ASSEMBLY/hifiasm/"$j".asm.fasta
	((j=j+1))
done	
cat $ASSEMBLY/hifiasm/*.asm.fasta >> $ASSEMBLY/raw.asm.fasta
awk '/^>/{print ">" ++i; next}{print}' < $ASSEMBLY/raw.asm.fasta > $ASSEMBLY/asm.fasta

minimap2 -x ava-pb $ASSEMBLY/asm.fasta $ASSEMBLY/asm.fasta > $ASSEMBLY/overlaps.paf

python3 ./from_clusters_to_assembly.py $OUT_DIR $REF_SIZE