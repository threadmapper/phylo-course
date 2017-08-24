import sys
import re
from collections import defaultdict
from itertools import groupby
  

def strMUT(text, dic):
    """ Replaces keys of dic with values of dic in text. 2005-02 by Emanuel Rumpf """
    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
    return re.sub(pat, lambda m: dic[m.group()], text)

templateo='''
'''
new = strMUT(templateo,  { 'SER': str(t)})


#----------------------------------------------------------------------------
def FastoE(fa_path="transcripts_all-simple.fasta"):
  found=0
  with open(fa_path) as f:
      groups = groupby(f, key=lambda x: not x.startswith(">"))
      seqD = {}
      for k,v in groups:
          if not k:
              key, val = list(v)[0].rstrip(), "".join(map(str.rstrip,next(groups)[1],""))
              gene=key[1:].strip()
              genox=gene.split()[0]## isoform 
              found+=1
              seqD[genox] = val.strip() ## sets
  return seqD.keys()

#---------------------------------------------------------------------------

template='''#!/bin/bash
#SBATCH -J stSER
#SBATCH -o stSER.stdout
#SBATCH -e stSER.stderr
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH -p nbi-long    # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-long

source kallisto-0.43.1
source samtools-1.5

## we made 180 20

PREFIX=FULL/TRIMMED/

kallisto quant --threads=16  --index=alyrata384v1.idx  --output-dir=OUT/OUTDIR    --plaintext  -b 100 --single   -l 500 -s 50     ${PREFIX}/trim-FASTQ

'''

BASE=['TBG_L1', 'TBG_H3', 'TBG_H1', 'TBG_L2', 'TBG_C1', 'TBG_C2', 'TBG_L3', 'TBG_H2', 'TBG_C3']
with open('submitter.sh', 'w') as outM:
 for t, fold in enumerate(BASE,1):
  with open('job'+ str(t) + '.py', 'w') as out:
    tag=fold.split('_')[1]
    new = strMUT(template,  { 'SER': fold, 'FASTQ': fold +'.fq', 'OUTDIR' : tag})
    out.write(new)
  outM.write('sbatch ' + 'job'+ str(t) + '.py' + '\n')

