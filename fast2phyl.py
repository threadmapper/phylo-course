import sys
## only keeping first 15 characters in header
## change param keep 

fa_path =  inFile = sys.argv[1]
ph_path = outFile = sys.argv[2]
keep=15
#fa_path = "unaligned.fasta" 
#ph_path = "unaligned.csv" 
#fa_path = "example_2_7.fasta" 
#ph_path = "example_2_7.csv" 

def sanit(rec):
   rec=rec.strip().split()
   st='_'.join(rec)
   st=st.replace(',','_')
   st=st.replace('-','_')
   st=st.replace('/','_')
   st=st.replace(' ','_')
   return st.strip()  

#------------------------------
from itertools import groupby
with open(fa_path) as f:
    groups = groupby(f, key=lambda x: not x.startswith(">"))
    d = {}
    for k,v in groups:
        if not k:
            key, val = list(v)[0].rstrip(), "".join(map(str.rstrip,next(groups)[1],""))
            d[sanit(key[1:])] = val     

tube=[]
maxfound=0
for hedo, seqo in d.items():
    tube.append([hedo,seqo])
    if len(seqo) > maxfound:
        maxfound=len(seqo)

#====== max length seq found
with open(ph_path, 'w') as outf:
   outf.write(" " + str(len(tube)) + " " + str(maxfound) + "\n")
   for h, seqo in tube:
     if len(seqo) > 3 and len(h) > 0:## 
      if len(seqo) < maxfound:# pad  
         pad=maxfound -len(seqo)
         #pad=0
         outf.write(h[:keep] + ' ' + seqo + '-'*(pad) +'\n' )
         #outf.write(h +  '\n' )
      else:
         outf.write(h[:keep]+ ' ' + seqo + '\n' )
         #outf.write(h +  '\n' )
      
      