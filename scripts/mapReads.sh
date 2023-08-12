#/bin/csh
# rm *.amb *.ann *.bwt *.pac *.sa *.fai
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# bwa, samtools, minimap2 (channel: bioconda)
# commands to install them:
# conda install -c bioconda bwa
# conda install -c bioconda samtools
# conda install -c bioconda minimap2

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
H=$OUT_DIR/hifi #HIFI

# ------- INDEX THE REFERENCE
bwa index $REF
samtools faidx $REF

minimap2 -ax map-pb -t $THREADS $REF $READS | samtools view -b -o $H.bam -
samtools sort -m10G -@ 20 -T /tmp/aln.sorted -o $H.sorted.bam $H.bam
samtools view -F 2308 -b $H.sorted.bam -o $H.sorted.nosecondary.bam # REMOVES SECONDARY ALIGNMENTS
samtools index $H.sorted.bam
samtools index $H.sorted.nosecondary.bam
python3 $NUCPLOT/NucPlot.py -w 32 $H.sorted.nosecondary.bam $H.png


# run the following commands only after executing all the commands above
# from the NucFreq plot determine the location of unresolved repeats in the reference genome
# indicated by the peaks, extract the reads mapping to the peaks and their neighborhoods within
# 50 Kb upstream and downstream whenever possible (unless you reach an end)


# # ------- SPECIFICATION
# CONTIG="id_of_the_contig_containing_unresolved_repeats"
# START="start_index_in_the_sequence_for_read_extraction"
# END="end_index_in_the_sequence_for_read_extraction"

# # ------- OUTPUT
# EXTRACTED_BAM="path_to_bam_file_for_selected_region_with_repeats"
# EXTRACTED_FASTQ="path_to_fastq_file_with_reads_from_selected_region_with_repeats"
# EXTRACTED_FASTA="path_to_fasta_file_with_reads_from_selected_region_with_repeats"

# # extracting and then converting bam to fastq and finally to fasta
# samtools index $H.sorted.bam
# samtools view -b $H.sorted.bam "$CONTIG:$START-$END" > $EXTRACTED_BAM.bam
# samtools fastq $EXTRACTED_BAM.bam > $EXTRACTED_FASTQ.fq
# cat $EXTRACTED_FASTQ.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $EXTRACTED_FASTA.fasta