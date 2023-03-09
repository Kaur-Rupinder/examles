#program for calculating position weights from fasta files consisting aligned sequences

import numpy as np, sys
import matplotlib.pyplot as plt

try:
  
  fasta = sys.argv[1]                                                       # Reading aligned fasta file of sequences from command line argument

except:
  
  print("Usages: python {aligned sequences fasta file.aln}")
  raise SystemExit()
 
seq_name=list()
seq_hash={}
counter=0
with open(fasta,'r') as fopen:
  for line in fopen:
        if line.startswith('>'):
            line=line.split('|')
            seq_name.append(line[3])
            counter=counter+1
        else:
            line=line.rstrip('\n')
            seq_hash.setdefault(seq_name[counter-1],[]).append(line)                                        #creating a hash where identifier is key and sequence is value as list
print("No. of uniq sequences in file=%d" % len(seq_hash))           


dict_pssm={}
for entry in seq_hash:
    seq_hash[entry]=(''.join(line.strip() for line in (seq_hash[entry])))                                   #joining sequence list as a single list

a=(seq_hash.keys()[-1])    
identifier=set()
for i in range(len(seq_hash[a])):                                                                           
                identifier.add(seq_hash[a][i])

if len(identifier)<=5:                                                                 #identifying type of sequence
   residues=['A','C','G','T','-']
else:
    residues=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
for i in range(len(seq_hash[a])):                                                                           #creating position specific hash from all sequences
        for j in seq_hash:
            dict_pssm.setdefault(i,[]).append(seq_hash[j][i])            


fout=open('pssm.txt','w')
fout.write('\t')
fout.write('\t'.join((i for i in residues)))
fout.write("\tSUM\tShannon_entropy\n")
total=float()
shannon_entropy=float()
shannon=list()
for i in (dict_pssm):
    fout.write((str(i)))
    fout.write("\t")
    element={}
    for j in range(len(dict_pssm[i])):
       element.setdefault(dict_pssm[i][j],[]).append(dict_pssm[i][j])
       total=0
       shannon_entropy=0
    for residue in residues:
        if residue in element:
          f=((float(len(element[residue]))/156))
          p=(((float(len((element[residue]))))+(0.25))/(float(156)+1.000))               # correction of frequency matrix by adding pseduo-weight 1.000
          entropy=(-f*(np.log2(f)))
          w=np.log2(p)                                                                               # calculating log of corrected frequency (weight) for each residue at each position
          fout.write("%.3f" % w)
          fout.write("\t")
          
        else:
            p=(float(0.25))/(float(156)+1.000)
            entropy=0
            w=np.log2(p)
            fout.write("%.3f\t" % w)
        shannon_entropy=shannon_entropy+entropy
        total=(total+float(w))                                                                                    #total weight of a column (position)
    shannon.append(shannon_entropy)   
    print >> fout,total,"\t",shannon_entropy
fout.close()
size = np.abs(shannon) * 100
plt.scatter(np.arange(0,len(shannon)),shannon, c='c', s=size)
plt.xlabel('<--  Position  -->',fontsize=40)
plt.ylabel('Entropy/Conservation',fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
