---
title: "Creating MAGs"
objectives:
- ""
- ""
- ""
keypoints:
- ""
- ""
---

 <a href="{{ page.root }}/fig/bioinformatic workflow.png">
  <img src="{{ page.root }}/fig/bioinformatic workflow.png" alt="Cog Metagenome" />
</a>

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

~~~

~~~
{: .bash}

> ## Important notes
>
> 
{: .callout}

<span style="color:blue">some *blue* text</span>.



## "Pipelines for binning a single sample with VAMB and performig CheckM on its resulting bins"
#Copy the code for each of the 3 scpripts in different files and run them where it corresponds (minimap and checkm in mazorka, vamb in local computer)

~~~
#!/bin/bash
#Script to be run in mazorka, it is to make the script to run minimap on one sample

<span style="color:blue">some *blue* text</span>.

FILE1=$1 #This one will be the metagenome scaffolds fasta file
FILE2=$2 #This one will be the raw forward reads, it can be compressed.gz
FILE3=$3 #This one will be the raw reverse reads, it can be compressed.gz
prefix=$4 #Sample name
root=$(pwd) # Sample directory with the 3 required files (for me is: /LUSTRE/usuario/czirion/zamia-dic2020/Zf_##)

cat > runMinimap${prefix}.sh <<E0F

#PBS -N minimap_${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e ${root}/LOGS_MAGS/minimap_${prefix}.error
#PBS -o ${root}/LOGS_MAGS/minimap_${prefix}.output
#PBS -V
 
module load samtools/1.9
module load minimap2/2.12 

cd $root
minimap2 -ax sr $FILE1 $FILE2 $FILE3 > $prefix.sam 
samtools view -S -b $prefix.sam > $prefix.bam
rm $prefix.sam

E0F
~~~
{: .code}

#once you run it run:
#qsub runMinimap${prefix}.sh 
~~~
mv runMinimap${prefix}.sh LOGS_MAGS

<span style="color:blue">some *blue* text</span>.
~~~
{: .code}
## The next script is to be run in the local computer once you have the bam file in mazorka ready, it moves the bam files to our local computer and runs vamb 

#!/bin/bash
prefix=$1 #sample name
root=$(pwd) # the folder were you want the direcory with vamb results, you must have the metagenome_scaffolds.fasta here (/home/claudia/Documetos/genomas/zamia-dic2020/vamb_ensam_mags/)

cat > runVAMB${prefix}.sh <<E0F

cd $root
scp czirion@mazorka.langebio.cinvestav.mx:/LUSTRE/usuario/czirion/zamia-dic2020/$prefix/$prefix.bam .
vamb --outdir #prefix/ --fasta ${prefix}_scaffolds.fasta --bamfiles $prefix.bam --minfasta 200000

E0F

##When vamb is finished the bam files must be compressed or erased because the use a lot of space:
#gzip $prefix.bam
##to move the resulting bins from the local computer to mazorka:
#cd $prefix/bins/
#scp *.fna czirion@mazorka.langebio.cinvestav.mx:/LUSTRE/usuario/czirion/zamia-dic2020/$prefix/vamb/



## Next script is for extracting the reads and reassembling the MAGs made with VAMB

#!/bin/bash

FILE1=$1 #forward reads
FILE2=$2 #reverse reads
prefix=$3 #sample name
root=$(pwd) #Sample folder with raw reads, METASPADES folder, and VAMB folder 
sign='$'

cat > runReassembly${prefix}.sh << E0F 

#PBS -N reassembly${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=32g,vmem=32g,walltime=100:00:00
#PBS -e $root/LOGS_MAGS/reassembly${prefix}.error
#PBS -o $root/LOGS_MAGS/reassembly${prefix}.output
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


# Enlista todos los bins.fna y le asigna a la variable file el nombre del archivo .fna
ls VAMB/*.fna | while read line; do file=${sign}(echo ${sign}line | cut -d'/' -f2);

# Le asigna a la variable forw el nombre del archivo de reads R1 pero sin la extension .fastq.gz y lo mismo para rev y lor reads R2
forw=${sign}(echo $FILE1| cut -d'.' -f1); rev=${sign}(echo $FILE2| cut -d'.' -f1);

# Corre craken para los contigs del bin
kraken2 --db kraken-db --threads 12 -input ${sign}line --output TAXONOMY_MAGS/${sign}file-kraken.kraken --report TAXONOMY_MAGS/${sign}file-kraken.report;

# Corre bracken
bracken -d kraken-db -i TAXONOMY_MAGS/${sign}file-kraken.report -o TAXONOMY_MAGS/${sign}file.bracken; 

# Hace los indices para le mapeo de los reads a los contigs del bin
bowtie2-build ${sign}line MAP_REASSEMBLY/INDEX/${sign}file; 

# Hace le mapeo
bowtie2 --threads 12 --sensitive-local -x  MAP_REASSEMBLY/INDEX/${sign}file -1 METASPADES/corrected/${sign}forw* -2 METASPADES/corrected/${sign}rev* -S  MAP_REASSEMBLY/FILE/${sign}file.sam;

# Convierte el mapeo a bam
samtools view -F 4 -bS  MAP_REASSEMBLY/FILE/${sign}file.sam >  MAP_REASSEMBLY/FILE/${sign}file.bam; 

# Convierte el mapeo a reads fastq
bamtools convert -in MAP_REASSEMBLY/FILE/${sign}file.bam -format fastq >  MAP_REASSEMBLY/FILE/${sign}file.fastq;

# Hace el directorio para los resultados del reensamble
mkdir SPADES_MAGS/${sign}file;

# Corre spades con los reads de un bin
spades.py -s MAP_REASSEMBLY/FILE/${sign}file.fastq -o SPADES_MAGS/${sign}file;

# Copia los scaffolds ensamblados, de su carpeta interna a la carpeta SPADES_MAGS, y le pone al principio el nombre de la muestra
scp SPADES_MAGS/${sign}file/scaffolds.fasta SPADES_MAGS/${sign}file.scaffolds.fasta 2>>/dev/null; done


E0F

#borrar los archivos .sam de MAP_REASSEMBLY/FILE
rm MAP_REASSEMBLY/FILE/*.sam 


## The next script is for running in mazorka CheckM on the vamb bins of one sample

#!/bin/bash

prefix=$1 #sample name
root=$(pwd) # directory above VAMB/ (where you put the bins for this sample in mazorka) (/LUSTRE/usuario/czirion/zamia-dic2020/$prefix/)

cat > runCheckm${prefix}.sh << E0F
#PBS -N runCheckm${prefix}
#PBS -q default
#PBS -l nodes=1:ppn=8,mem=48g,vmem=48g,walltime=100:00:00
#PBS -e $root/LOGS_MAGS/checkm${prefix}.error
#PBS -o $root/CHECKM/checkm${prefix}.output
#PBS -V

module load CheckM/1.1.3
module load hmmer/3.1b2
module load Prodigal/2.6.2

cd $root/
mkdir CHECKM
checkm lineage_wf -r VAMB/ CHECKM/

E0F

#once all your .fna files are in the vamb directory in mazorka you can run this script:
#qsub runCheckM${prefix}.sh 
#mv runCheckm${prefix}.sh LOGS_MAGS
