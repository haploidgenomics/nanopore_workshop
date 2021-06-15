# Nanopore Workshop

In this workshop, you will learn how to use Snakemake, a workflow management tool that enhances reprodicibility and scalability, to create an assembly of XX genome using the Nanopore Minion platform.

This workshop assumes that you have:
* A Mac/OS system
* Basic unix command knowledge using Terminal on a Mac
* A folder containing /fastq files basecalled from the MinION
* For this workshop, fastq files were generated from fast5 files base-called with Guppy 4.0 (there might be newer versions of Guppy available, which is updated by ONT on a regular basis, so check to use the latest version for best performance/accuracy)

## Installing Dependencies

Open terminal on your Mac and create a snakemake directory.
```
mkdir snakemake
cd snakemake
```

**Install Mambaforge.**
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```
Answer "yes" to first question about license.  
Answer "no" to second question about activation (This prevents anaconda from activating without you knowing). If you accidentally answered "yes" you can undo it by typing the command below.
```
conda config --set auto_activate_base false
```
Next, close the terminal window.  
Then, reopen a new terminal window to activate conda and install Snakemake.

**Install Snakemake.**
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

**Install Atom using the link below.**

https://atom.io/

Atom is a text editor that will allow you to edit the Snakefile as well as the env file to run Snakemake while having a Terminal window open at the same time. Open the Atom editor after installation.  

You will also have to install a **Terminal** package in Atom. To do this, click on "Install a Package" and search for "Terminal", and install the package called "terminal-tab 0.6.0". Now you are ready to run a Snakemake Workflow!

## Running a Snakemake Workflow

First, you'll have to download the folder containing the the fastq files from the strain that you've been assigned to assemble. Put this folder under the /snakemake directory. In this example we will be working with the /s42 strain. There are three general steps:
* Creating a filtered set of fastq files for assembly (using filtlong)
* Creating a draft assembly (using flye)
* Polishing the draft assembly (using racon and medaka)
* Evaluating the assembly (using quast)

You will follow along to see how this Snakemake Workflow is implemented so that you can use it for other applications as well.

**Create Project in Atom**
Open a new window on Atom, and click "Add Project Folder", then choose the /s42 folder (or the folder you're assigned to). After this, you should also launch a Terminal window in Atom by clicking "Packages" and then "Terminal". Terminal will appear at the bottom of the screen.

Next, you will need to create two files to run snakemake. Using Atom, create a new file called "Snakefile". Then create a new folder called "env" and a file called "nanopore.yaml" within the env folder. Click on the Snakefile and nanopore.yaml file so you can edit them in Atom.

Add these dependencies in the nanopore.yaml file, and hit save.
```
channels:
  - conda-forge
  - bioconda
dependencies:
 - bwa =0.7.17
 - samtools =1.10
 - filtlong =0.2.0
 - minimap2 =2.17
 - bwa =0.7.17
 - biopython =1.77
 - python =3.6
 - flye =2.8.3
 - racon =1.4.20
 - medaka =1.0.3
```

**Concatenate and Filter Fastq File**
Add these commands to the Snakefile, and hit save.
```
FASTQ, = glob_wildcards("data/fastq/{sample}.fastq")

rule fastq_merge:
    input:
        expand("data/fastq/{sample}.fastq", sample=FASTQ)
    output:
        "output/merged.fastq"
    shell:
        "cat {input} > {output}"

rule filtlong_reads:
    input:
        "output/merged.fastq"
    output:
        "output/merged_filtlong.fastq"
    conda:
        "env/nanopore.yaml"
    shell:
        "filtlong --min_length 3000 --min_mean_q 30 --length_weight 10 -p 90 "
        "--target_bases 500000000 {input} > {output}"
```
You will then perform a dry run, followed by regular run using the snakemake command.
```
conda activate snakemake
snakemake --use-conda --cores 2 output/merged_filtlong.fastq -n
```
Run the command without the "-n" command next.

**Create a Draft Assembly using flye**
Add these commands to the Snakefile, and hit save.
```
rule flye:
    input:
        "output/merged_filtlong.fastq"
    output:
        dir=directory("output/flye"),
        assembly="output/flye/assembly.fasta"
    conda:
        "env/nanopore.yaml"
    shell:
        "flye --nano-raw {input} --out-dir {output.dir} --genome-size 3m --threads 8 --iterations 5"
```

You will then perform a dry run, followed by regular run using the snakemake command.
```
conda activate snakemake
snakemake --use-conda --cores 2 output/fyle -n
```
Run the command without the "-n" command next.



