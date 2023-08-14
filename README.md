# RAmbler

This repository is a stand-alone copy of RAmbler. Anyone can run RAmbler locally in their machines just by cloning this.

# Pre-requisites

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
</pre>

Currently, RAmbler only takes reads in `.fasta` format. If you have the reads in a different file format, please convert them to `.fasta` before executing the assembly pipeline.

# Scripts

`mapReads.sh` contains the required commands to perform step A (Determine the reads corresponding to repetitive regions)  
`extractUnikmers.sh` contains the required commands to perform step B (Determine unikmers)  
`run.sh` contains the required commands to execute RAmbler's steps C, D, E and F (the assembly pipeline)  

If you already have the HiFi reads corresponding to the repetitive regions and the unikmers, directly execute `run.sh` to obtain the repeat resolved assembly.  
If you don't have them, then first execute `mapReads.sh` by providing the current draft genome and the complete set of HiFi reads to it. Follow the detailed steps mentioned in `mapReads.sh` to extract reads mapping to the repetitive regions. For executing `mapReads.sh`, you will also need a local copy of [NucFreq](https://github.com/mrvollger/NucFreq). Next, for extracting the unikmers, execute `extractUnikmers.sh`. Follow the detailed instructions mentioned in `extractUnikmers.sh`.

# Execution

