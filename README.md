# Applied Bioinformatics in Agriculture and Medicine
Organized by [Association of Nepalese Agricultural Professionals of Americas (NAPA)](http://www.napaamericas.org) and [Agriculture and Forestry University](http://afu.edu.np/).

July 19-23, 2020

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
# aCtivate conda

```
 conda activate bioinformatics
 ```
