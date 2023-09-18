#/bin/csh
# rm *.amb *.ann *.bwt *.pac *.sa *.fai
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# bwa, samtools, minimap2, pysam (channel: bioconda), numpy, pandas, matplotlib, seaborn (channel: conda-forge)
# commands to install them:
# conda install -c bioconda bwa
# conda install -c bioconda samtools
# conda install -c bioconda minimap2
# conda install -c bioconda pysam
# conda install -c conda-forge numpy
# conda install -c conda-forge pandas
# conda install -c conda-forge matplotlib
# conda install -c conda-forge seaborn

# It requires the following GitHub repository available for visualizing the coverage plot
# https://github.com/mrvollger/NucFreq

# ------- COMMAND
NUCPLOT="path_to_NucFreq_GitHub_repository"

# ------- DATA
READS="path_to_the_read_file"
REF="path_to_the_reference_genome"

# ------- CONFIGURATION
THREADS=20

# ------- OUTPUT
OUT_DIR="path_to_output_directory"
Q=$OUT_DIR/WGS/hifi #HIFI

mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/WGS
mkdir -p $OUT_DIR/contigs

# command for converting the contig ids to the expected format (contig$i)
awk '/^>/{print ">contig" ++i; next}{print}' < $REF > $OUT_DIR/WGS/genome_renamed.fasta

# ------- INDEX THE REFERENCE
bwa index $REF
samtools faidx $REF

minimap2 -ax map-pb -t $THREADS $OUT_DIR/WGS/genome_renamed.fasta $READS | samtools view -b -o $Q.bam -
samtools sort -m10G -@ 20 -T /tmp/aln.sorted -o $Q.sorted.bam $Q.bam
samtools view -F 2308 -b $Q.sorted.bam -o $Q.sorted.nosecondary.bam # REMOVES SECONDARY ALIGNMENTS
samtools index $Q.sorted.bam
samtools index $Q.sorted.nosecondary.bam

WGS_BAM=$OUT_DIR/WGS/hifi.sorted.bam
CONTIG_COUNT="Put the number of contigs/chromosomes present in the genome"

for ((i=1; i<=CONTIG_COUNT; i++));
do
    H=$OUT_DIR/contigs/$i
    mkdir -p $H
    $SAM view -b $WGS_BAM "contig"$i > $H"/contig"$i.bam
    $SAM view -F 2308 -b $H"/contig"$i.bam -o $H/nosecondary."contig"$i.bam
    $SAM index $H"/contig"$i.bam
    $SAM index $H/nosecondary."contig"$i.bam
    python3 $NUCPLOT/NucPlot.py -w 32 $H/nosecondary."contig"$i.bam $H"/contig"$i.png
done

# run the following commands only after executing all the commands above
# from the NucFreq plot determine the location of unresolved repeats in the reference genome
# indicated by the peaks, extract the reads mapping to the peaks and their neighborhoods within
# 50 Kb upstream and downstream whenever possible (unless you reach an end)


# # ------- SPECIFICATION
# CONTIG="id_of_the_contig_containing_unresolved_repeats"
# START="start_index_in_the_reference_sequence_for_read_extraction"
# END="end_index_in_the_reference_sequence_for_read_extraction"

# # ------- OUTPUT
# EXTRACTED_BAM="path_to_bam_file_for_selected_region_with_repeats"
# EXTRACTED_FASTQ="path_to_fastq_file_with_reads_from_selected_region_with_repeats"
# EXTRACTED_FASTA="path_to_fasta_file_with_reads_from_selected_region_with_repeats"

# # extracting and then converting bam to fastq and finally to fasta
# samtools index $H.sorted.bam
# samtools view -b $H.sorted.bam "$CONTIG:$START-$END" > $EXTRACTED_BAM.bam
# samtools fastq $EXTRACTED_BAM.bam > $EXTRACTED_FASTQ.fq
# cat $EXTRACTED_FASTQ.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $EXTRACTED_FASTA.fasta