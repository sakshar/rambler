# RAmbler

This repository is a stand-alone copy of RAmbler. Anyone can run RAmbler locally in their machines just by cloning this.

# Pre-requisites
You will need conda installed in your machine to proceed. If you don't have conda installed, please follow [conda installation](https://conda.io/docs/user-guide/install/). If you already had it installed, please make sure it is updated.

Create a conda environment first:

`conda create -n rambler python=3.8 numpy pandas matplotlib seaborn xlsxwriter biopython`

This will create a conda environment named rambler with the afforementioned packages. Then activate the environment with:

`conda activate rambler`

Finally, install the following packages inside the newly created environment using the commands below:

<pre>
conda install -c bioconda jellyfish
conda install -c bioconda minimap2
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda pysam
conda install -c bioconda hifiasm
conda install -c bioconda seqkit
</pre>

Currently, RAmbler only takes reads in `.fasta` format. If you have the reads in a different file format, please convert them to `.fasta` before executing the assembly pipeline.

# Scripts

`mapReads.sh` contains the required commands to perform step A (Determine the reads corresponding to repetitive regions)  
`extractUnikmers.sh` contains the required commands to perform step B (Determine unikmers)  
`run.sh` contains the required commands to execute RAmbler's steps C, D, E and F (the assembly pipeline)  
`stitch.sh` contains the required commands to stitch the repeat resolved assembly back to the original draft/reference genome

If you already have the HiFi reads corresponding to the repetitive regions and the unikmers, directly execute `run.sh` to obtain the repeat resolved assembly.  
If you don't have them, then first execute `mapReads.sh` by providing the current draft genome and the complete set of HiFi reads to it. Follow the detailed steps mentioned in `mapReads.sh` to extract reads mapping to the repetitive regions. For executing `mapReads.sh`, you will also need a local copy of [NucFreq](https://github.com/mrvollger/NucFreq). Next, for extracting the unikmers, execute `extractUnikmers.sh`. Follow the detailed instructions mentioned in `extractUnikmers.sh`.  
If you have already executed `run.sh` and obtained the repeat resolved assembly, execute `stitch.sh` to merge the newly reconstructed repeats with the draft genome.

# Execution

Open the terminal and follow the steps below to run RAmbler from scratch:
<pre>
  cd ~
  git clone https://github.com/sakshar/rambler.git
  git clone https://github.com/mrvollger/NucFreq.git
</pre>

We have provided a test dataset inside `data` directory:  
`data/10.fasta` contains the HiFi reads generated with a 10x coverage depth  
`data/10.kmers` contains the unikmers  
`data/10000_5_100.fasta` contains the true reference genome for the above HiFi reads and unikmers. It contains 5 copies of a 10 Kb repeat unit with a mutation rate of 1 per 100 bases along with 50 Kb upstream and downstream flanking regions  

Always run RAmbler from inside the `rambler/scripts` directory. So, to run RAmbler directly when you have the required HiFi reads and unikmers, execute the following command:
<pre>
  cd ~
  cd rambler/scripts
  bash run.sh -r read_file -u unikmer_file -o output_directory -s assembly_length [-k kmer_length] [-l tolerance] [-h threshold]
  required:
  -r path to the fasta file containing the HiFi reads extracted in step A
  -u path to the file containing the list of unikmers extracted in step B
  -o path to the directory where RAmbler will put all the generated outputs
  -s estimated length of the repeat resolved part of the genome in bp
  optional:
  -k length of kmers used (default: 21)
  -l tolerance parameter (default: 15)
  -h threshold parameter (default: 15)
</pre>

The repeat resolved assembly will be generated in a file named `rambler.fasta` inside `output_directory/assembly` directory.

Next, execute the following command to merge `rambler.fasta` with the draft reference genome for obtaining the final assembly:
<pre>
  bash stitch.sh -g draft_genome_file -n contig_id -b start_position -e end_position -o output_directory -s contig_length
  required:
  -g path to the fasta file containing the draft/reference genome
  -n id of the contig that is being re-assembled
  -b approximate start position of the repeat region in the contig sequence from the draft genome (if stitching at the end, put -1)
  -e approximate end position of the repeat region in the contig sequence from the draft genome (if stitching at the beginning, put -1)
  -o path to the directory where RAmbler will put all the generated outputs (always put the same one as run.sh)
  -s estimated length of the repeat resolved contig in bp
</pre>

# Example run

`bash run.sh -r ../data/10.fasta -u ../data/10.kmers -o ../output -s 150000`  

This will create an output directory `rambler/output` and the reconstruced repeats will be in `rambler/output/assembly/rambler.fasta`.

`bash stitch.sh -g ../data/10000_5_100.fasta -n draft -b 50000 -e 100000 -o ../output -s 150000`

This will create a sub-directory `final` inside the previously created output directory `rambler/output` and the repeat reconstructed contig will be in `rambler/output/final/rambler_merged.fasta`.

# Citation
- Sakshar Chakravarty, Glennis Logsdon, Stefano Lonardi. RAmbler: de novo genome assembly of complex repetitive regions. bioRxiv 2023.05.26.542525; doi: https://doi.org/10.1101/2023.05.26.542525  
or
- Sakshar Chakravarty, Glennis Logsdon, and Stefano Lonardi. 2023. RAmbler: de novo genome assembly of complex repetitive regions. In 14th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics (BCB ’23), September 3–6, 2023, Houston, TX, USA. ACM, New York, NY, USA, 1 page. https://doi.org/10.1145/3584371.3612971 

# FAQs
1. How to calculate the lower and upper bounds of the window for extracting unikmers from the generated excel file?

While calculating for the mean and the standard deviation from the excel file, please ignore the first 5 rows (which are mostly k-mers due to sequencing errors) and also ignore the k-mers that appear way higher than the expected coverage depth (suppose, if the expected coverage depth is 100x, ignore rows after row 250). For example, with a read set having a coverage depth of 100x, try to calculate the mean and the standard deviation with values from rows 6 to 250. This is a trick required to discard k-mers coming from sequencing errors and the k-mers which appear way higher than expected that can throw the values of mean and standard deviation way off.

Note: the scripts, `freq_distribution_kmers.py` and `extractUnikmers.sh`, have been updated to fully automate step B for the extraction of unikmers

# TODO
A stand-alone conda package for RAmbler with a single command execution feature
