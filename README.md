# Applied Bioinformatics in Agriculture and Medicine

Website: http://www.napaamericas.org/bioinformatics.php

Organized by [Association of Nepalese Agricultural Professionals of Americas (NAPA)](http://www.napaamericas.org) and [Agriculture and Forestry University](http://afu.edu.np/).

July 19-23, 2020
____
Instructors: 
- **Dr. Ananta Acharya** (ananta.acharya at gmail dot com)
- **Sishir Subedi** (subedisishir at gmail dot com)
- **Dr. Dev Paudel** (merodev at gmail dot com)
- **Dr. Saroj Parajuli** (sarose97 at gmail dot com)
____
**DATA AND SCRIPTS FOR TRAINING NEWBIES IN BIOINFORMATICS**

### Software used
For the training session, following software will be used. These are pre-installed on the instance of Google Cloud that we will be using. For your personal use later, you will need to download and install these software into your computer.

- SRA-toolkit: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
- BWA: http://bio-bwa.sourceforge.net/
- Fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Samtools: http://www.htslib.org/download/
- IGV: http://software.broadinstitute.org/software/igv/download
- MUSCLE: https://www.drive5.com/muscle/

### Computing session
For this session, we will be using [Google Cloud Platform](https://cloud.google.com/).

### Practise bash
https://repl.it/languages/bash

https://linuxconfig.org/bash-scripting-tutorial-for-beginners

### Commonly used commands for practice
https://bioinformaticsworkbook.org/Appendix/Unix/UnixCheatSheet.html#gsc.tab=0

```
ls       # list the files
ls -l    # detailed list
ls -ltr  # list in reverse chronological order
mkdir corona # make directory
cd corona # change directory
pwd # present working directory
cd ..
ls
cd
less filename.txt # view file
head filename.txt # head of file
head -n 5 covid.fasta # top 5 lines
tail covid.fasta # bottom of file
```

```
# Copy source destination
cp /home/guest/Downloads/covid_samples/covid_samples.fasta covid.fasta
```

Redirect
```
grep ">" covid_samples.fasta > covid_sample_names.txt
```
Append it

```
grep ">" covid_samples.fasta >> covid_sample_names.txt
```
Count number of lines

```
wc covid_sample_names.txt # line, word, character
wc -l covid_sample_names.txt # lines
```

### Download data from NCBI

```
wget https://sra-pub-sars-cov2.s3.amazonaws.com/sra-src/SRR11177792/WHV-Nepal-61-TW_1.fastq.gz
fastq-dump --outdir fastq --gzip -F --split-files SRR11177792
```
### Variant / SNP calling pipeline
![SNP calling pipeline](https://github.com/dpaudel/NepalBioinformatics/blob/master/snp_calling_pipeline.PNG?raw=true)

# Seminar Day 4: Commands used

## Login directions for using MobaXterm
https://1drv.ms/b/s!AnZvwUbCCvF77XErrqD1HaawPfd8?e=mPKRAU

## Login directions for Linux

Install putty

```
brew install putty
```
Convert .ppk file to openssh

```
puttygen key.ppk -O private-openssh -o key.pem
```
Use ssh to connect:

```
ssh -i key.pem afu@104.154.53.116
```
If permission error is seen then use the following code to fix prior to ssh:

```
chmod 400 key.pem
```

Use given passphrase to login to the google cloud.

## Start the analysis

### Go to your directory

```
cd name12345
```

### Download raw data

```
# Do not run this as it takes a lot of storage and time. Filter only 100,000 reads as shown in the next command.
# fastq-dump --outdir fastq --gzip -F --split-files SRR11177792 
```
### Download only 100,000 reads

```
fastq-dump SRR11177792 --split-files -N 10000 -X 110000
```
### Copy reference genome to your directory from parent directory

```
cp /data/rawdata/NC_045512.2.fasta .
```
### Index the reference (genome) sequence

```
bwa index NC_045512.2.fasta
```
### Check quality of raw data

```
fastqc SRR11177792_1.fastq
fastqc SRR11177792_2.fastq 
```
View sample of fastqc report: https://drive.google.com/file/d/15MVcmcjlLPH82NZDpYSQdnLuWIqy5rrN/view?usp=sharing

### Align raw reads with bwa-mem to the reference genome

```
bwa mem /data/rawdata/NC_045512.2.fasta SRR11177792_1.fastq SRR11177792_2.fastq > SRR11177792.sam
```

### Check the alignment (sam) file

```
head SRR11177792.sam
```
### Convert sam to binary bam file to save storage space

```
samtools view -S -b SRR11177792.sam > SRR11177792.bam
```

### Sort bam file
```
samtools sort SRR11177792.bam -o SRR11177792.sorted.bam
```
### Index bam file to view in IGV viewer

```
samtools index SRR11177792.sorted.bam
```
### View statistics of bam file

```
samtools flagstat SRR11177792.sorted.bam
```

### Download IGV viewer

http://software.broadinstitute.org/software/igv/download

### Call variants

```
bcftools mpileup -f /data/rawdata/NC_045512.2.fasta SRR11177792.sorted.bam | bcftools call -mv -Ov -o SRR11177792.variants.txt
```

```
samtools view /data/rawdata/covid_samples.sorted.bam NC_045512.2:21563-25384 > Sgen.sorted.bam
```

# Multiple alignment

```
muscle -in /data/rawdata/Sgen.fa -out Sgen.2.mfa
```
### View as html

```
 muscle -in /data/rawdata/Sgen.fa -out Sgen.2.mfa -html
 ```
 ### Make cluster using MUSCLE
 
 ```
 muscle -maketree -in /data/rawdata/Sgen.mfa -out covid.phy -cluster neighborjoining
 ```
 
 View phylogeny: http://etetoolkit.org/treeview/
 
