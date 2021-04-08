---
title: "Creating MAGs"
objectives:
- "Assemble a metagenome with metaSPAdes"
- "Bin the contigs of your metagenome"
- "Reassemble the reads from each bin to create the MAGs"
- "Check the completeness and contamination of your MAGs"
---

> ## Important notes
> The scripts are a template to create a script addapted to your sample. You have to run this scripts puting certain arguments in the command line in order (each scripts tells you what arguments). The output will be a script with the information of your sample ready to be run by the cluster (`qsub runScript_mysample.sh`). 


# Script to run metaSPAdes (in the cluster)
This script will run the assembly of the metagenome and copy the resulting scaffolds to the working directory. The output will be in the METASPADES folder.
```bash
nano metaspades.sh 
```
Copy the next code and paste it inside `metaspades.sh`
```bash
#!/bin/bash
FILE1=$1 # Forward reads, a fastq.gz file
FILE2=$2 # Reverse reads, a fastq.gz file
prefix=$3 # Sample ID, it will be appended to the beginning of the files
root=$(pwd) # Your working directory must be the one where you have the required files and you want the output folders to be created  

cat > runMetaspades_${prefix}.sh <<EOF
#PBS -N metaspades_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=40g,vmem=40g,walltime=100:00:00
#PBS -e ${root}/LOGS/metaspades_${prefix}.error
#PBS -o ${root}/LOGS/metaspades_${prefix}.output
#PBS -V

module load SPAdes/3.10.1

cd $root

mkdir LOGS/

metaspades.py --pe1-1 $FILE1 --pe1-2 $FILE2 -o METASPADES

scp METASPADES/scaffolds.fasta ./${prefix}_scaffolds.fasta 2>>/dev/null 

EOF
```
Now run the script with the required arguments in order and submit the resulting script
```bash
sh metaspades.sh READS_R1.fastq.gz READS_R2.fastq.gz sample01 #The output is the runMetaspades_sample01.sh
qsub runMetaspades_sample01.sh
```
Once it has finished running and you have your assembly, run the next script.

# Script to run minimap2 (in the cluster)
This script will perform the alignment of the reads to the assembled metagenome. We need this step to do the binning. After it does the alignment it converts the SAM file to a BAM file (te one we need) and erases the SAM file, because it is very heavy. 
```bash
nano minimap.sh
```
```bash
#!/bin/bash
FILE1=$1 # The metagenome scaffolds fasta file
FILE2=$2 # Forward reads, a fastq.gz file
FILE3=$3 # Reverse reads, a fastq.gz file
prefix=$4 # Sample ID, it will be appended to the beginning of the files
root=$(pwd) 

cat > runMinimap$_{prefix}.sh <<E0F

#PBS -N minimap_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e ${root}/LOGS/minimap_${prefix}.error
#PBS -o ${root}/LOGS/minimap_${prefix}.output
#PBS -V
 
module load samtools/1.9
module load minimap2/2.12 

cd $root
minimap2 -ax sr $FILE1 $FILE2 $FILE3 > $prefix.sam 
samtools view -S -b $prefix.sam > $prefix.bam
rm $prefix.sam

E0F
```
Run it

```bash
sh minimap.sh sample01_scaffolds.fasta READS_R1.fastq.gz READS_R2.fastq.gz sample01
```
```bash
qsub runMinimap_sample01.sh
```
Once it is finished you have to move the BAM file and the scaffolds file to your local computar to run VAMB (unless your cluser can run VAMB, Mazorka can´t do it). It will take a while. 

# Running VAMB in the local computer
Once you have the BAM and the scaffolds in your local computer run:

```bash
vamb --outdir sample01/ --fasta sample01_scaffolds.fasta --bamfiles sample01.bam --minfasta 200000
```
It will take a while.
Your bins will be in `sample01/bins/` and will have names like this: `4657.fna`, `15777.fna`
Now you can erase the BAM file if you want because it is somewhat heavy.
Now create a folder called `VAMB/` in your working directory in the cluster and copy your bins there (the `fna` files, not the `bins/` folder)

# Script for extracting the reads for each bin, reassembling them and perform a taxonomic assignation
This script will use bowtie2, samtools and bamtools to map the reads (corrected by metaSPAdes) to each bin and deliver a FASTQ file with the reads. This reads are going to be used by SPAdes to assemble the MAGs. It will also give you a taxonomic assignation for each bin. 

```bash
nano reassembly.sh
```
```bash
#!/bin/bash
FILE1=$1 # Forward reads, a fastq.gz file
FILE2=$2 # Reverse reads, a fastq.gz file
prefix=$3 # Sample ID, it will be appended to the beginning of the files
root=$(pwd) # Your working directory must be the one where you have the required files and you want the output folders to be created 
sign='$'    

cat > runReassembly_${prefix}.sh << E0F 

#PBS -N reassembly_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e $root/LOGS/reassembly_${prefix}.error
#PBS -o $root/LOGS/reassembly_${prefix}.output
#PBS -V

module load bowtie2/2.3.5.1
module load samtools/1.9
module load bwa/0.7.15
module load bamtools/2.4.1
module load SPAdes/3.10.1
module load kraken/2.0.7
module load Braken/2.0

cd $root

mkdir MAP_REASSEMBLY
mkdir MAP_REASSEMBLY/INDEX
mkdir MAP_REASSEMBLY/FILE
mkdir TAXONOMY_MAGS


ls VAMB/*.fna | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f2);
forw=${sign}(echo $FILE1| cut -d'.' -f1); rev=${sign}(echo $FILE2| cut -d'.' -f1);

kraken2 --db kraken-db --threads 12 -input ${sign}line --output TAXONOMY_MAGS/${sign}file-kraken.kraken --report TAXONOMY_MAGS/${sign}file-kraken.report;
bracken -d kraken-db -i TAXONOMY_MAGS/${sign}file-kraken.report -o TAXONOMY_MAGS/${sign}file.bracken; 

bowtie2-build ${sign}line MAP_REASSEMBLY/INDEX/${sign}file; 
bowtie2 --threads 12 --sensitive-local -x  MAP_REASSEMBLY/INDEX/${sign}file -1 METASPADES/corrected/${sign}forw* -2 METASPADES/corrected/${sign}rev* -S  MAP_REASSEMBLY/FILE/${sign}file.sam;
samtools view -F 4 -bS  MAP_REASSEMBLY/FILE/${sign}file.sam >  MAP_REASSEMBLY/FILE/${sign}file.bam; 
bamtools convert -in MAP_REASSEMBLY/FILE/${sign}file.bam -format fastq >  MAP_REASSEMBLY/FILE/${sign}file.fastq;

mkdir SPADES_MAGS/${sign}file;
mkdir MAGS/;
spades.py -s MAP_REASSEMBLY/FILE/${sign}file.fastq -o SPADES_MAGS/${sign}file;
scp SPADES_MAGS/${sign}file/scaffolds.fasta MAGS/${sign}file.scaffolds.fasta 2>>/dev/null; done

E0F
```
```bash
sh reassembly.sh READS_R1.fastq.gz READS_R2.fastq.gz sample01
```
```bash
qsub runReassembly_sample01.sh
```
Once it is finished you can erase the SAM file:
```bash
rm MAP_REASSEMBLY/FILE/*.sam 
```
In the SPADES_MAGS/ folder you will have the assembled MAGs scaffolds and a folder for each MAG with all of the SPAdes outputs.
In the TAXONOMY_MAGS/ folder you will have the Kraken and Bracken reports with the taxonomy

# Script for running CheckM on the MAGs 
This script will run the lineage workflow of CheckM, that will give you an estimate for the completeness, the contamination and the heterogeneity of your MAGs.
The result you want is a table that is diplayed in the output log `LOGS/checkm_sample01.output`.
```bash
nano checkm.sh 
```
```bash
#!/bin/bash
prefix=$3 # Sample ID, it will be appended to the beginning of the files
root=$(pwd) # # Your working directory must be the one where you have the required files and you want the output folders to be created 

cat > runCheckm_${prefix}.sh << E0F
#PBS -N runCheckm_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=48g,vmem=48g,walltime=100:00:00
#PBS -e $root/LOGS_MAGS/checkm_${prefix}.error
#PBS -o $root/CHECKM/checkm_${prefix}.output
#PBS -V

module load CheckM/1.1.3
module load hmmer/3.1b2
module load Prodigal/2.6.2

cd $root/
mkdir CHECKM
checkm lineage_wf -x fasta -r MAGS/ CHECKM/

E0F
```
```bash
sh checkm.sh sample01
```
```bash
qsub runCheckm_sample01.sh
```
# You finished!

