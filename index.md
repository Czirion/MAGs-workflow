
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
root=$(pwd) # Your working directory must be the one where you have the required files and wherer the outputs will be created  
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


cat > runReassembly_${prefix}.sh << E0F2 

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


ls VAMB/*.fna | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f
2);
forw=${sign}(echo $FILE1| cut -d'.' -f1); rev=${sign}(echo $FILE2| cut -d'.' -f1
);

kraken2 --db kraken-db --threads 12 -input ${sign}line --output TAXONOMY_MAGS/${
sign}file-kraken.kraken --report TAXONOMY_MAGS/${sign}file-kraken.report;
bracken -d kraken-db -i TAXONOMY_MAGS/${sign}file-kraken.report -o TAXONOMY_MAGS
/${sign}file.bracken; 

bowtie2-build ${sign}line MAP_REASSEMBLY/INDEX/${sign}file; 
bowtie2 --threads 12 --sensitive-local -x  MAP_REASSEMBLY/INDEX/${sign}file -1 M
ETASPADES/corrected/${sign}forw* -2 METASPADES/corrected/${sign}rev* -S  MAP_REA
SSEMBLY/FILE/${sign}file.sam;
samtools view -F 4 -bS  MAP_REASSEMBLY/FILE/${sign}file.sam >  MAP_REASSEMBLY/FI
LE/${sign}file.bam; 
bamtools convert -in MAP_REASSEMBLY/FILE/${sign}file.bam -format fastq >  MAP_RE
ASSEMBLY/FILE/${sign}file.fastq;

mkdir SPADES_MAGS/${sign}file;
mkdir MAGS/;
spades.py -s MAP_REASSEMBLY/FILE/${sign}file.fastq -o SPADES_MAGS/${sign}file;
scp SPADES_MAGS/${sign}file/scaffolds.fasta MAGS/${sign}file.scaffolds.fasta 2>>
/dev/null; done

E0F2

cat > runCheckm_${prefix}.sh << E0F3
#PBS -N runCheckm_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=48g,vmem=48g,walltime=100:00:00
#PBS -e $root/LOGS/checkm_${prefix}.error
#PBS -o $root/CHECKM/checkm_${prefix}.output
#PBS -V

module load CheckM/1.1.3
module load hmmer/3.1b2
module load Prodigal/2.6.2

cd $root/
mkdir CHECKM
checkm lineage_wf -x fasta -r MAGS/ CHECKM/

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
runReassembly_sample01.sh
```
And now we will run them in order with a lot of waiting in the middle. 

# Running metaSPAdes
This script will run the assembly of the metagenome (the output will be in the `METASPADES/` folder) and copy the resulting scaffolds to the working directory.
```bash
qsub runMetaspades_sample01.sh
```
Now ypu should have all of this in your sample folder:
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
Once it is finished you have to move the BAM file and the scaffolds file from the cluster to the same folder in your local computar to run VAMB (unless your cluser can run VAMB). It will take a while. 

# Running VAMB in the local computer
Once you have the BAM and the scaffolds in your local computer run:

```bash
vamb --outdir sample01/ --fasta sample01_scaffolds.fasta --bamfiles sample01.bam --minfasta 200000
```
It will take a while.
Your bins will be in `sample01/bins/` and will have names like this: `171479.fna`, `171530.fna`
Now create a folder called `VAMB/` in your working directory in the cluster and copy your bins there (the `fna` files, not the `bins/` folder)

# Extracting the reads for each bin, reassembling them and performing a taxonomic assignation
This script will use bowtie2, samtools and bamtools to map the reads (corrected by metaSPAdes) to each bin and deliver a FASTQ file with the reads. This reads are going to be used by SPAdes to assemble the MAGs. It will also give you a taxonomic assignation for the contigs of each bin. 

```bash
qsub runReassembly_sample01.sh
```
Once it is finished you can erase the SAM file:
```bash
rm MAP_REASSEMBLY/FILE/*.sam 
```
In the `SPADES_MAGS/` folder you will have a folder for each MAG with all of the SPAdes outputs:
```bash
171799.fna/
    assembly_graph.fastg
    assembly_graph.gfa
    before_rr.fasta
    contigs.fasta
    contigs.paths
    corrected/
    dataset.info
    input_dataset.yaml
    K21/
    K33/
    K55/
    K77/
    misc/
    params.txt
    scaffolds.fasta
    scaffolds.paths
    spades.log
    tmp/
    warnings.log

```
The scaffolds FASTA files will also be in the folder called `MAGS/`.

In the `TAXONOMY_MAGS/` folder you will have the Kraken and Bracken reports with the taxonomy.
```bash
171479.fna.bracken
171479.fna-kraken_bracken.report
171479.fna-kraken.kraken
171479.fna-kraken.report
171530.fna.bracken
171530.fna-kraken_bracken.report
171530.fna-kraken.kraken
171530.fna-kraken.report
171732.fna.bracken
171732.fna-kraken_bracken.report
171732.fna-kraken.kraken
171732.fna-kraken.report
171799.fna.bracken
171799.fna-kraken_bracken.report
171799.fna-kraken.kraken
171799.fna-kraken.report
```

# Running CheckM on the MAGs 
This script will run the lineage workflow of CheckM, that will give you an estimate for the completeness, the contamination and the heterogeneity of your MAGs.

```bash
qsub runCheckm_sample01.sh
```
The result you want is a table that is diplayed in the output log `LOGS/checkm_sample01.output`. It should look something like this: 
```bash
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                    Marker lineage      # genomes   # markers   # marker sets   0     1    2   3   4   5+   Completeness   Contamination   Strain heterogeneity  
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  171732.fna.scaffolds    k__Archaea (UID2)        207         145           103        36   102   7   0   0   0       75.42            3.62              14.29          
  171799.fna.scaffolds   k__Bacteria (UID203)      5449        104            58        96    6    2   0   0   0       12.93            2.59              100.00         
  171530.fna.scaffolds   k__Bacteria (UID203)      5449        104            58        84    3    6   3   2   6       10.33           20.88              21.35          
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


```

# You finished!
