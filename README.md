# Nanopore Workshop

In this workshop, you will learn how to use Snakemake, a workflow management tool that enhances reprodicibility and scalability, to create an assembly of XX genome using the Nanopore Minion platform.

This workshop assumes that you have:
* A Mac/OS system
* Basic unix command knowledge using Terminal on a Mac
* Set of folders containing matching /fast5 and /fastq files

## Installing Dependencies

Open terminal on your Mac and create a snakemake directory.
```
mkdir snakemake
cd snakemake
```

Install Mambaforge.
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```

Install Snakemake.
```
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda deactivate
```

Test to see if Snakemake was installed.
```
conda activate snakemake
snakemake --help
conda deactivate
```

Install Atom using the link below.

https://atom.io/

Atom is a text editor that will allow you to edit the Snakefile as well as the env file to run Snakemake while having a Terminal window open at the same time.

## Running a Snakemake Workflow


