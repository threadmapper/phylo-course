import sys
import itertools
import re
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import pprint 
from itertools import groupby
from dotmap import DotMap

def avg(v):
  return 1.0*sum(v)/len(v)



def strMUT(text, dic):
    """ Replaces keys of dic with values of dic in text. 2005-02 by Emanuel Rumpf """
    pat = "(%s)" % "|".join(map(re.escape, dic.keys()))
    return re.sub(pat, lambda m: dic[m.group()], text)



def getFast(fastaFiler = "93-11_AltReferenceGenome.fasta"):
    return SeqIO.to_dict(SeqIO.parse(fastaFiler, "fasta"))


def getEasy(f):
    for header, group in itertools.groupby(f, key = lambda x: x.startswith(">") ):
      if header:
         line = next(group)
         tag = line[1:].strip().split()
      else:
         sequence = ''.join(line.strip() for line in group)
         yield tag, sequence.strip()

def getSeqD( fastaFiler = "93-11_AltReferenceGenome.fasta" ):
   """
   chinaD = loadWholeSeq("9311_combinedREF.fasta")
   """
   seqD = {}
   tube =[] 
   with open(fastaFiler, 'r' ) as inpFast:
    for h, seq in getEasy(inpFast):
        seqD[tuple(h)]= seq
        tube.append(tuple(h))
   return tube, seqD

def getTubeD( fastaFiler = "93-11_AltReferenceGenome.fasta" ):
   """
   chinaD = loadWholeSeq("9311_combinedREF.fasta")
   """
   seqD = {}
   tube =[] 
   with open(fastaFiler, 'r' ) as inpFast:
    for h, seq in getEasy(inpFast):
        tag = ''.join(h)
        seqD[tag] = seq
        tube.append(tag)
   return tube, seqD

class Viv(dict):
  def __missing__(self, key):
      val = self[key] = type(self)()
      return val  
 