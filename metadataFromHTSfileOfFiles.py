#!/usr/bin/python
import numpy as np
import scipy as sp
import pysam
import pandas as pd
import os
import re
import math
from datetime import date
import sys
import argparse
import types

sequencer_re_dict={
    'MiSeq': re.compile('A[AZ]XX$'),
    'HiSeq_2000': re.compile('A[BC]XX$'),
    'HiSeq_2500': re.compile('A[DN]XX$'),
    'HiSeq_4000': re.compile('B[BC]X[XY]$'),
    'HiSeq_X': re.compile('ALXX$|CCX[2XY]$'),
    'NextSeq': re.compile('AFX[XY]$|BGX[0-9A-Z]$'),
    'NovaSeq': re.compile('D[A-Z]XX$')
}

target_re=re.compile('Agilent|agilent|RNA_SEQ|Twist|twist|WHOLE_GENOME|nimblegen|fgen|nextera|rapidCapture|exome3|xgen')

def flowcell2sequencer(suffix_4):
    seqs=[]
    for k in sequencer_re_dict:
        v=sequencer_re_dict[k]
        if re.search(v,suffix_4) is not None:
            seqs.append(k)
    return sorted(list(set(seqs)))

class htsFile(pysam.libcalignmentfile.AlignmentFile):
    def about(self):
        return "class <htsfile> extends <pysam.libcalignmentfile.AlignmentFile>\nAs an extension of pysam this provides additional methods for an htsLib file [BAM/CRAM]\n"
    def contigs(self):
        headerObj=self.header.to_dict()
        contigs=[]
        for rg in headerObj['SQ']:
            try:
                ctg=rg['SN']
                contigs.append(ctg)
                contigs=list(set(contigs()))
            except:
                pass
        return contigs
    def assembly(self):
        headerObj=self.header.to_dict()
        assembly='.'
        try:
            seqs=self.contigs()
            assembly=seqs.pop()
        except:
            pass
        return assembly
    def centers(self):
        headerObj=self.header.to_dict()
        centers=[]
        for rg in headerObj['RG']:
            try:
                ctr=rg['CN']
                centers.append(ctr)
                centers=list(set(centers))
            except:
                pass
        if len(centers)==0:
            centers='.'
        return centers
    def libraries(self):
        headerObj=self.header.to_dict()
        libraries=[]
        for rg in headerObj['RG']:
            try:
                lbr=rg['LB']
                libraries.append(lbr)
                libraries=list(set(libraries))
            except:
                pass
        if len(libraries)==0:
            libraries='.'
        return libraries
    def samples(self):
        headerObj=self.header.to_dict()
        samples=[]
        for rg in headerObj['RG']:
            try:
                smp=rg['SM']
                samples.append(smp)
                samples=list(set(samples))
            except:
                pass
        if len(samples)==0:
            samples='.'
        return samples
    def readLen(self):
        read_len='.'
        for i,r in enumerate(self):
            if i > 10:
                break
            read_len=len(r.seq)
            if read_len != '.':
                break
        return read_len
    def idxStats(self):
        return self.get_index_statistics()
    def date(self):
        headerObj=self.header.to_dict()
        return headerObj['RG'][0]['DT']
    def sex(self,**kwargs):
        logC=kwargs.get('logC',-4.0)
        xi=kwargs.get('Xidx',22)
        yi=kwargs.get('Yidx',23)
        X=self.idxStats()[xi][1]
        Y=self.idxStats()[yi][1]
        r=Y/X
        logR=math.log(r)
        sx='.'
        if logR > logC:
            sx='M'
        elif logR <=logC:
            sx='F'
        return sx
    def sexStats(self,**kwargs):
        logC=kwargs.get('logC',-4.0)
        xi=kwargs.get('Xidx',22)
        yi=kwargs.get('Yidx',23)
        X=self.idxStats()[xi][1]
        Y=self.idxStats()[yi][1]
        r=Y/X
        logR=math.log(r)
        sx='.'
        if logR > logC:
            sx='M'
        elif r <=logC:
            sx='F'
        return (sx,X,Y,logR)
    def idxStats(self):
        return self.get_index_statistics()
    def flowCells(self):       
        headerObj=self.header.to_dict()
        suffix=[]
        for rg in headerObj['RG']:
            suff=rg['ID'].split('.')[0][-4:]
            for k2 in sequencer_re_dict:
                regExp=sequencer_re_dict[k2]
                if re.search(regExp,suff) is not None: 
                    suffix.append(suff[-4:])
        suffix=list(set(suffix))
        return suffix   
    def sequencers(self):
        sequencer=[]
        suffixes=self.flowCells()
        for suff in suffixes:
 #           print(suff)
            seqncr=flowcell2sequencer(suff)
            sequencer.extend(seqncr)
        sequencer=list(set(sequencer))
        return sequencer
    def target(self):
        headerObj=self.header.to_dict()
        target='.'
        for pg in headerObj['PG']:
            if re.search(re.compile('GATK'),pg['ID']) is not None:
                try:
                    for q1 in pg['CL'].split(' '):
                        if re.search('target',q1) is not None:
                            q2=q1.split('=')
                            q3=q2[1].split('/')
                            target='.'
                            for v in q3:
                                if re.search(target_re,v) is not None:
                                    target=v
                except:
                    target='.'
            return target
    def bwa_version(self):
        headerObj=self.header.to_dict()
	version='.'
        for pg in headerObj['PG']:
            if ((re.match(re.compile('bwa'),pg['ID']) is not None) | (re.match(re.compile('STAR'),pg['ID']) is not None:)):
                try:
                    version=pg['VN']
                except:
                    version='.'
        return version

def metadataFromHTSfileOfFiles(fof,**kwargs):
    grp=kwargs.get('Group','group')
    result=[]
    count=0
    json=[]
    paths = open(fof,"r")
    for path in paths:
        count+=1
        limsId=os.path.basename(path).split('.')[0]
#        print(str(count) + "\t" + path + "\t" + limsId)
        path=path.strip()
        f=open(path,'rb')
        htsObj=htsFile(f,'rb')
        attr_dict={}
        attr_dict['count']=count
        attr_dict['group']=grp
        attr_dict['limsID']=limsId
        attr_dict['Date']=htsObj.date()
        attr_dict['flowcell_suffix']=htsObj.flowCells()
        attr_dict['target']=htsObj.target()
        attr_dict['sequencer']=htsObj.sequencers()
        attr_dict['bwa_version']=htsObj.bwa_version()
        attr_dict['read_len']=htsObj.readLen()
        attr_dict['sample_sex']=htsObj.sex()
        attr_dict['centers']=htsObj.centers()
        attr_dict['libraries']=htsObj.libraries()
        attr_dict['samples']=htsObj.samples()
        attr_dict['contigs']=htsObj.contigs()
        attr_dict['assembly']=htsObj.assembly()

        json.append(attr_dict)        
        f.close()
    paths.close()
    df=pd.DataFrame.from_dict(json, orient='columns')
    df=df[['count','group','limsID','Date','flowcell_suffix','target','sequencer','bwa_version','read_len','sample_sex','centers','libraries','samples','assembly']]
    return (json,df)

def main():
    parser = argparse.ArgumentParser(prog='metadataFromHTSfileOfFiles.py',description='Collects metadata from a file that contains paths of bam files one path per line.')
    parser.add_argument('-i', nargs='?',default=sys.stdin,type=str, help='input a file of HTS file paths, one full bam path, line by line in a text file. Must be named with the convention <ID>.*.bam (e.g. <123>.merged.bam). Must have an index in the same directory')
    parser.add_argument('-l', nargs='?',default='test', type=str, help='group label defaults to <test>')
    parser.add_argument('-t',  nargs='?',default=sys.stdout, type=str, help='output tsv metadata File: defaults to stdout')
    parser.add_argument('-j',  nargs='?',default=False, type=str, help='output json metadata File: defaults to nothing; error if -o is also specified for stdout')

    args = parser.parse_args()
    if len(sys.argv) < 1:
        parser.print_help()
        sys.exit(1)
    f=args.i
    l=args.l
    o=args.t
    j=args.j

    if o==j:
        parser.print_help()
        sys.exit(1)
    result=metadataFromHTSfileOfFiles(f,Group=l)
    json=result[0]
    df=result[1]
    if o==sys.stdout:
        df.to_csv(o,sep="\t",header=True,index=False)
    else:
        try:
            df.to_csv(o,sep="\t",header=True,index=False)
        except:
            sys.stderr.write("ERROR: could not write " + o + "\n")	
    if j==sys.stdout:
        jsonStr=df.to_json(orient='index')
        sys.stdout.write(jsonStr)
    elif j!=False:
        try:
            out=open(j,"w")
            jsonStr=df.to_json(orient='index')
            out.write(jsonStr)
            out.close()
        except:
            sys.stderr.write("ERROR: could not write " + j + "\n")	
            out.close()
    else:
        pass

if __name__ == '__main__':
#	main()
	try:
		main()
	except OSError:
		sys.stdout.write("RUN ERROR: Are the arguments correct? Try runnning metadataFromHTSfileOfFiles.py  -h\n")
