# HTSmetadata
Introspect a High Throughput Sequenceing (HTS) set of files and generate tsv and json datasets 
# Usage 
## metadataFromHTSfileOfFiles.py [-h] [-i [I]] [-l [L]] [-t [T]] [-j [J]]

### Collects a dataset of metadata directed by a file that contains paths of hts(bam/cram) files one path perline.
## Example:
python metadataFromHTSfileOfFiles.py -i ~/tmp/paths.txt -t ~/tmp/test.tsv -j ~/tmp/jsonTest.out -l 'testCmd'
## optional arguments:
##  -h, --help  
show this help message and exit.
##  -i [I]      
input a file of HTS file paths, one full bam path, line by line in a text file. Must be named with the convention ID.*.bam (e.g. <123>.merged.bam). Must have an index in the same directory.
##  -l [L]      
group label defaults to <test>.
##  -t [T]      
  output a tsv metadata File: defaults to stdout.
##  -j [J]      
  output a json metadata File: defaults to nothing; error if -o is also specified for stdout or equal to -j.
# Requirements:
## python 3+
numpy,scipy,pysam,pandas,os,re,math,datetime,sys,argparse,types
## htslib v1.9+, samtools v1.9+
