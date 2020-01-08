# HTSmetadata
Introspect a High Throughput Sequenceing (HTS) set of files and generate tsv and json datasets 
# Usage 
## metadataFromHTSfileOfFiles.py [-h] [-i [I]] [-l [L]] [-t [T]] [-j [J]]

### Collects a dataset of metadata directed by a file that contains paths of hts(bam) files one path per line.
## Example:
python metadataFromHTSfileOfFiles.py -i ~/tmp/paths.txt -t ~/tmp/test.tsv -j ~/tmp/jsonTest.out -l 'testCmd'
## optional arguments:
     -h, --help  
     show this help message and exit.
     -i [I]      
     input a file of HTS file paths, one full bam path, line by line in a text file. Must be named with the convention ID.*.bam (e.g. <123>.merged.bam). Must have an index in the same directory but do not list the index files in the text.
     -l [L]      
     group label defaults to <test>.
     -t [T]      
     output a tsv metadata File: defaults to stdout.
     -j [J]      
     output a json metadata File: defaults to do nothing; -j <file>
# Requirements:
## python 3+
numpy,scipy,pysam,pandas,os,re,math,datetime,sys,argparse,types
## htslib v1.9+, samtools v1.9+
# Output file key:
    column 1: count	
    column 2: group
    column 3: limsID
    column 4: Date
    column 5: flowcell_full
    column 6: flowcell_suffix
    column 7: target
    column 8: sequencer (type not id)	
    column 9: aligner_version
    column 10: read_len
    column 11: sample_sex
    column 12: centers
    column 13: libraries
    column 14: samples
    column 15: assembly (requires AS tag)

# Flowcell Suffix to Machine Dictionary:
## (Machines are identified by mapping the suffix of the flowcell)
    'MiSeq': 'A[AZ]XX$',
    'HiSeq_2000': 'A[BC]XX$',
    'HiSeq_2500': 'A[DN]XX$',
    'HiSeq_4000': 'B[BC]X[XY]$',
    'HiSeq_X': 'ALXX$|CCX[2XY]$',
    'NextSeq': 'AFX[XY]$|BGX[0-9A-Z]$',
    'NovaSeq': 'D[A-Z]XX$'

