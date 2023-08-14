#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# jellyfish (channel: bioconda), xlsxwriter (channel: conda-forge)
# command to install jellyfish:
# conda install -c bioconda jellyfish
# conda install -c conda-forge xlsxwriter

# ------- SPECIFICATION
# READ_FILE="path_to_HiFi_reads_fasta_file"
# k="kmer_length"
# OUT_COUNT="path_to_output_from_count_command"
# OUT_HISTO="path_to_output_from_histo_command"
# OUT_XCEL="path_to_output_from_freq_distribution_kmers.py"
# UNIKMERS="path_to_output_from_dump_command"
READ_FILE="../data/10000_5_100.fasta"
k=21
OUT_COUNT="../data/test.jellyfish"
OUT_HISTO="../data/test.histo"
OUT_XCEL="../data/test.xlsx"
UNIKMERS="../data/test.kmers"

# the count and histo commands together generate the frequency distribution of kmers
# in the $OUT_HISTO file 
jellyfish count -m $k -C -o $OUT_COUNT -c 3 -s 10000000 --disk -t 16 $READ_FILE
jellyfish histo -o $OUT_HISTO -v $OUT_COUNT

# the python script generates an excel (.xlsx) file to visualize the distribution
# using this excel file, you can find out the [mean - 3*sd, mean + 3*sd] window
# containing the expected unikmers
python3 ./freq_distribution_kmers.py $OUT_HISTO $OUT_XCEL


# only after calculating the mean and sd from the excel file, run the following command
# to get the final list of expected unikmers

# lower="put round(mean - 3*sd) here as expects an integer value"
# upper="put round(mean + 3*sd) here as expects an integer value"
# jellyfish dump -c -t -L $lower -U $upper $OUT_COUNT | awk '{print $1}' > $UNIKMERS
lower=1
upper=1
jellyfish dump -c -t -L $lower -U $upper $OUT_COUNT | awk '{print $1}' > $UNIKMERS