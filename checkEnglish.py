from nltk.corpus import wordnet
import os
import sys

defPath = "/home/tom/projects/UGP2/datafiles/dictionaryCorpus"
fin = open(defPath+"/rr",'r').read()
fout = open(defPath+"/rw",'w')
out=""
spp=[' ','\n','\t']
a=0
b=0
word=""
for c in fin:
	if c in spp:
		if wordnet.synsets(word):
			print a
			a+=1
			print>>fout,word
		word=""
	else:
		word=word+c
		#print word