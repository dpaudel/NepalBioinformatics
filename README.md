# Applied Bioinformatics in Agriculture and Medicine
Organized by [Association of Nepalese Agricultural Professionals of Americas (NAPA)](http://www.napaamericas.org) and [Agriculture and Forestry University](http://afu.edu.np/).

July 19-23, 2020

**Data and scripts for training newbies in bioinformatics**

### Computing session
For this session, we will be using [Google Cloud Platform](https://cloud.google.com/).

### Practise bash
https://repl.it/languages/bash

```
ls
ls -l
ls -ltr
mkdir corona
cd corona
pwd
cd ..
ls
cd
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
