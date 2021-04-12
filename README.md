# Obtaining MAGs (Metagenome-Assembled Genomes) from shotgun metagenomic sequencing data
On this [link](https://czirion.github.io/MAGs-workflow/) you will find a pipeline to obtain Metagenome-Assembled-Genomes(MAGs) from shotgun metagenomic sequencing data and a quality evaluation.
You will need:
  1. Access to a high-performance computing cluster.
  2. A computer with administrator permission where you can access a UNIX terminal.
  3. Shotgun metagenomic sequencing reads whose quality you already approve.

## Instructions for the use of this pipeline
The pipeline contains a script that should be run in the cluster and additional commands to be run in the cluster and in the local computer. When you run the script it creates a colection of scripts adapted to the information of your sample. If you use the same cluster I used (Mazorka Langebio Cinvestav Mexico) the required software is already installed (latest versions on March 2021), but if you are using a different machine make sure it has the required software and adapt the script to the way it runs jobs. 

### Software that must be installed in the cluster
*The numbers show the version I used, if you use a different version make sure it is compatible with the rest of the pipeline.* 

[SPAdes](https://github.com/ablab/spades) 3.10.1

[kraken](https://github.com/DerrickWood/kraken2) 2.0.7

[Braken](https://github.com/jenniferlu717/Bracken) 2.0

[bowtie2](https://github.com/BenLangmead/bowtie2) 2.3.5.1

[samtools](https://github.com/samtools/samtools) 1.9

[bwa](https://github.com/lh3/bwa) 0.7.15

[bamtools](https://github.com/pezmaster31/bamtools) 2.4.1

[minimap2](https://github.com/lh3/minimap2) 2.12

[CheckM](https://github.com/Ecogenomics/CheckM) 1.1.3

[hmmer](https://github.com/EddyRivasLab/hmmer) 3.1b2

[Prodigal](https://github.com/hyattpd/Prodigal) 2.6.2


### Software that must be installed in the local computer 
[VAMB](https://github.com/RasmussenLab/vamb)





	
