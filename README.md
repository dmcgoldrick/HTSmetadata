# HTSmetadata
Introspect a High Throughput Sequenceing (HTS) set of files and generate tsv and json datasets 
# Usage 
## metadataFromHTSfileOfFiles.py [-h] [-i [I]] [-l [L]] [-t [T]] [-j [J]]

### Collects a dataset of metadata directed by a file that contains paths of hts(bam) files one path per line.
## Example:
python metadataFromHTSfileOfFiles.py -i ~/tmp/paths.txt -t ~/tmp/test.tsv -j ~/tmp/jsonTest.out -l 'testCmd'
## optional arguments:
##  -h, --help  
show this help message and exit.
##  -i [I]      
input a file of HTS file paths, one full bam path, line by line in a text file. Must be named with the convention ID.*.bam (e.g. <123>.merged.bam). Must have an index in the same directory but do not list the index files in the text.
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
# Example output
### column 1: count	
### column 2: group
### column 3: limsID
### column 4: Date
### column 5: flowcell_suffix
### column 6: target
### column 7: sequencer	
### column 8: bwa_version
### column 9: read_len
### column 10: sample_sex
### column 11: centers
### column 12: libraries
### column 13: samples
### column 14: assembly
  
1	testCmd	102859	2014-08-09T01:20:44-0700	['ACXX', 'ADXX']	nimblegen_solution_V2refseq_2010	['HiSeq_2000', 'HiSeq_2500']	0.7.10-r789	50	F	['University_of_Washington_Genome_Sciences']	['27-27365_B:2']	['102859']	hs37d5

2	testCmd	112610	2014-08-27T08:47:52-0700	['ACXX', 'ADXX']	nimblegen_solution_V2refseq_2010	['HiSeq_2000', 'HiSeq_2500']	0.7.10-r789	50	F	['University_of_Washington_Genome_Sciences']	['27-28750_F:10', '27-28697_F:10']	['112610']	hs37d5

3	testCmd	112625	2014-10-17T03:37:37-0700	['ANXX']	nimblegen_solution_V2refseq_2010	['HiSeq_2500']	0.7.10-r789	50	F	['University_of_Washington_Genome_Sciences']	['27-30015_B:5', '27-30021_B:5']	['112625']	hs37d5
