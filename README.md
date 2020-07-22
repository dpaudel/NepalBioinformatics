# Applied Bioinformatics in Agriculture and Medicine
Organized by [Association of Nepalese Agricultural Professionals of Americas (NAPA)](http://www.napaamericas.org) and [Agriculture and Forestry University](http://afu.edu.np/).

July 19-23, 2020
Instructors: Dr. Ananta Acharya, Sishir Subedi, Dr. Dev Paudel, Dr. Saroj Parajuli

**Data and scripts for training newbies in bioinformatics**

### Computing session
For this session, we will be using [Google Cloud Platform](https://cloud.google.com/).

### Practise bash
https://repl.it/languages/bash

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
fastq-dump -I --split-files SRR11177792
```

# Seminar Day 4: Commands used

## Login directions for using MobaXterm
https://1drv.ms/b/s!AnZvwUbCCvF77XErrqD1HaawPfd8?e=mPKRAU

### Go to your directory

```
cd name12345
```

### Download raw data

```
fastq-dump SRR11177792 --split-files
```

### Copy reference genome to your directory from parent directory

```
cp /data/rawdata/NC_045512.2.fasta .
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

### Download IGV viewer

http://software.broadinstitute.org/software/igv/download
