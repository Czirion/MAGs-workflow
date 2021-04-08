# Obtaining MAGs (Metagenome-Assembled-Genomes) from shotgun metagenomic sequencing data
This repository contains a pipeline to obtain Metagenome-Assembled-Genomes(MAGs) from shotgun metagenomic sequencing data.
You will need:
  1. Access to a high-performance computing cluster.
  2. A computer with administrator permission where you can access a UNIX terminal.
  3. Shotgun metagenomic sequencing reads whose quality you already approve.

## Instructions for the use of this pipeline
The *mags_pipeline.sh* file contains a collection of scripts to be run in the cluster and additional commands to be run in the cluster and in the local computer, it is not a script that runs automatically. The scripts are meant for you to copy them in a text editor in the cluster you are using but in this pipeline I show them to you mixed with more commands that you will need to run both in your computer and in the cluster in order. If you use the same cluster I used (Mazorka Langebio Cinvestav Mexico) the required software is already installed (latest versions on March 2021) and the scripts are ready to be pasted in a text editor in the cluster and run them; but if you are using a different cluster make sure to adapt the scripts to the way your cluster runs them. 

## The output you will have at the end of the pipeline

  


### Software that must be installed in the cluster
*The numbers show the version I used, if you use a different version make sure it is compatible with the rest of the pipeline* 

[SPAdes]() 3.10.1

[kraken]() 2.0.7

[Braken]() 2.0

[bowtie2]() 2.3.5.1

[samtools]() 1.9

[bwa]() 0.7.15

[bamtools]() 2.4.1

[minimap2]() 2.12

[CheckM]() 1.1.3

[hmmer]() 3.1b2

[Prodigal]() 2.6.2


### Software that must be installed in the local computer 
[VAMB](https://github.com/RasmussenLab/vamb)





	
