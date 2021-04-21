
## Instructions to use the pipeline
What the script you will find here does is to create several scripts adapted to the information of your sample (you should run it several times if you have more than one sample). You will also find the instructions of the order in which you should run the created scripts and additional steps. 
Unfortunately VAMB is not able to run in the cluster that I used, so the pipeline is not as automatic as it could be, and you will need to run VAMB in your local computer. This involves to wait for some heavy files to copy from the cluster to your computer. 

# Script
You should have a folder for each of your samples with the files of your forward and reverse reads in a gzip-compressed format.
The folder of the sample you are working with will be your working directory, in there you will create the script: 
```bash
nano mags_pipeline.sh 
```
And paste inside it the following code:
```bash
#!/bin/bash
FILE1=$1 # Forward reads, a fastq.gz file
FILE2=$2 # Reverse reads, a fastq.gz file
prefix=$3 # Sample ID, it will be appended to the beginning of the files
root=$(pwd) # Your working directory must be the one where you have the required files and where the outputs will be created  
sign='$'

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

cat > runMinimap_${prefix}.sh <<E0F1

#PBS -N minimap_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e ${root}/LOGS/minimap_${prefix}.error
#PBS -o ${root}/LOGS/minimap_${prefix}.output
#PBS -V
 
module load samtools/1.9
module load minimap2/2.12 

cd $root
minimap2 -ax sr ${prefix}_scaffolds.fasta $FILE1 $FILE2 > $prefix.sam 
samtools view -S -b $prefix.sam > $prefix.bam
rm $prefix.sam

E0F1


cat > runTaxonomy_${prefix}.sh <<E0F2 

#PBS -N taxonomy_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e $root/LOGS/taxonomy_${prefix}.error
#PBS -o $root/LOGS/taxonomy_${prefix}.output
#PBS -V

module load kraken/2.0.7
module load Braken/2.0

cd $root

mkdir TAXONOMY_MAGS

ls VAMB/*.fna | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f2);
kraken2 --db kraken-db --threads 12 -input ${sign}line --output TAXONOMY_MAGS/${sign}file_kraken.kraken --report TAXONOMY_MAGS/${sign}file_kraken.report;
bracken -d kraken-db -i TAXONOMY_MAGS/${sign}file_kraken.report -o TAXONOMY_MAGS/${sign}file.bracken; done

E0F2

cat > runCheckm_${prefix}.sh <<E0F3
#PBS -N runCheckm_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=48g,vmem=48g,walltime=100:00:00
#PBS -e $root/LOGS/checkm_${prefix}.error
#PBS -o $root/LOGS/checkm_${prefix}.output
#PBS -V

module load CheckM/1.1.3
module load hmmer/3.1b2
module load Prodigal/2.6.2

cd $root/
mkdir CHECKM
checkm lineage_wf -x fasta -r MAGS/ CHECKM/
checkm qa CHECKM/lineage.ms CHECKM/ --file CHECKM/quality_${prefix}.tsv --tab_table -o 2

E0F3
```
As you can see in lines 2, 3 and 4, you need to run the script giving 3 arguments in order: the forward reads, the reverse reads and the sample ID:
```bash
sh mags_pipeline.sh READS_R1.fastq.gz READS_R2.fastq.gz sample01 
```
This created our 4 scripts:
```bash
runCheckm_sample01.sh
runMetaspades_sample01.sh
runMinimap_sample01.sh
runTaxonomy_sample01.sh
```
And now we will run them in order with a lot of waiting in the middle. 

# Running metaSPAdes
This script will run the assembly of the metagenome (the output will be in the `METASPADES/` folder) and copy the resulting scaffolds to the working directory.
```bash
qsub runMetaspades_sample01.sh
```
Now you should have all of this in your sample folder:
```bash
LOGS/
    metaspades_sample01.output
    metaspades_sample01.error
mags_pipeline.sh
METASPADES/
    assembly_graph.fastg
    assembly_graph.gfa
    before_rr.fasta
    contigs.fasta
    contigs.paths
    corrected/
    dataset.info
    first_pe_contigs.fasta
    input_dataset.yaml
    K21/
    K33/
    K55/
    misc/
    params.txt
    scaffolds.fasta
    scaffolds.paths
    spades.log
    tmp/
runCheckm_sample01.sh
runMetaspades_sample01.sh
runMinimap_sample01.sh
runReassembly_sample01.sh
sample01_R1.fastq.gz
sample01_R2.fastq.gz
sample01_scaffolds.fasta

```
# Running minimap2
This script will perform the alignment of the reads to the assembled metagenome. We need this step to do the binning. After it does the alignment it converts the SAM file to a BAM file (the one we need) and erases the SAM file because it is very heavy. 
```bash
qsub runMinimap_sample01.sh
```
Once it is finished you have to move the BAM file and the scaffolds file from the cluster to a folder in your local computer to run VAMB (unless your cluser can run VAMB). It will take a while. 

# Running VAMB in the local computer
I installed VAMB using Conda and the version I used is 3.0.2.
Once you have the BAM and the scaffolds in your local computer run:

```bash
vamb --outdir sample01/ --fasta sample01_scaffolds.fasta --bamfiles sample01.bam --minfasta 200000
```
With the `--outdir` flag you specify the name of the folder where it is going to put the output. With the `--fasta` and `--bamfiles` flags you give the program the assembly and the mapping files. The `--minfasta` option is to tell the program to give you FASTA files of the bins that are at least 200,000 bases long, if you don't use this options VAMB will give you information about the bins it identified but not the bins themselves. 
It will take a while.

You will have an output like this: 
```bash
sample01/
    bins/
        145.fna
        62132.fna
    clusters.tsv
    latent.npz
    lengths.npz
    log.txt
    mask.npz
    model.pt
    rpkm.npz
    tnf.npz
```
Your bins are the `145.fna`, `62132.fna` files (`.fna` is an extension equivalent to `.fasta`).
Now create a folder called `VAMB/` in your working directory in the cluster and copy your bins there (the `fna` files, not the `bins/` folder)

# Taxonomic assignation 
This script will give you a taxonomic assignation for the contigs of each bin. 

```bash
qsub runTaxonomy_sample01.sh
```
In the `TAXONOMY_MAGS/` folder you will have the Kraken and Bracken reports with the taxonomy.
```bash
145.fna.bracken
145.fna_kraken_bracken.report
145.fna_kraken.kraken
145.fna-kraken.report
62132.fna.bracken
62132.fna_kraken_bracken.report
62132.fna_kraken.kraken
62132.fna_kraken.report
```

# Running CheckM on the MAGs 
This script will run the lineage workflow of CheckM, that will give you an estimate for the completeness, the contamination and the heterogeneity of your MAGs, and other statistics like genome size and GC %.

```bash
qsub runCheckm_sample01.sh
```
The results are in the `CHECKM/quality_sample01.tsv`file. 


# You finished!
