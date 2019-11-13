# HTSmetadata
Introspect a High Throughput Sequenceing set of files for samples and generate tsv and json datasets 
# Usage 
## metadataFromHTSfileOfFiles.py [-h] [-i [I]] [-l [L]] [-t [T]] [-j [J]]

### Collects metadata from a file that contains paths of bam files one path perline.

## optional arguments:
##  -h, --help  ### show this help message and exit
##  -i [I]      ###input a file of HTS file paths, one full bam path, line by line in a text file. Must be named with the convention ID.*.bam
###              (e.g. <123>.merged.bam). Must have an index in the same directory
##  -l [L]      ### group label defaults to <test>
##  -t [T]      ### output tsv metadata File: defaults to stdout
##  -j [J]      ### output json metadata File: defaults to nothing; error if -o is also specified for stdout or equal
