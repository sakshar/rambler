#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# jellyfish (channel: bioconda), xlsxwriter (channel: conda-forge)
# command to install them:
# conda install -c bioconda jellyfish
# conda install -c conda-forge xlsxwriter

# ------- SPECIFICATION
READ_FILE="path_to_HiFi_reads_fasta_file"
k="kmer_length"
OUT_COUNT="path_to_output_from_count_command"
OUT_HISTO="path_to_output_from_histo_command"
OUT_XCEL="path_to_output_from_freq_distribution_kmers.py"
UNIKMERS="path_to_output_from_dump_command"

# the count and histo commands together generate the frequency distribution of kmers
# in the $OUT_HISTO file 
jellyfish count -m $k -C -o $OUT_COUNT -c 3 -s 10000000 --disk -t 16 $READ_FILE
jellyfish histo -o $OUT_HISTO -v $OUT_COUNT

# the python script generates an excel (.xlsx) file to visualize the distribution
# using this excel file, you can find out the [mean - 3*sd, mean + 3*sd] window
# containing the expected unikmers
# the python script also outputs the calculated mean, sd, lower bound and upper bound
# values from the input distribution file
OUTPUT=$(python3 ./freq_distribution_kmers.py $OUT_HISTO $OUT_XCEL)

echo $OUTPUT

# extract the lower and upper bound values from the output
lower=$(echo "$OUTPUT" | sed -n '2p' | cut -d':' -f2 | xargs)
upper=$(echo "$OUTPUT" | sed -n '3p' | cut -d':' -f2 | xargs)

# extract the unikmers using the lower and upper bound values
jellyfish dump -c -t -L $lower -U $upper $OUT_COUNT | awk '{print $1}' > $UNIKMERS
