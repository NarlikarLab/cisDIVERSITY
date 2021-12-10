import sys
import re
import os

def reverseComp(seq):
	toReturn = ""
	for i in range(len(seq)):
		if(seq[i] == 'N' or seq[i] == 'n' ):
			toReturn = 'N' + toReturn
		elif(seq[i] == 'A' or seq[i] == 'a'):
			toReturn = 'T' + toReturn
		elif(seq[i] == 'C' or seq[i] == 'c'):
			toReturn = 'G' + toReturn
		elif(seq[i] == 'G' or seq[i] == 'g'):
			toReturn = 'C' + toReturn
		elif(seq[i] == 'T' or seq[i] == 't'):
			toReturn = 'A' + toReturn
		else:
			print "Something wrong",seq
			sys.exit()
	return(toReturn)

def convertToNumbers(seq):
	out = seq.replace("A","0 ")
	out = out.replace("C","1 ")
	out = out.replace("G","2 ")
	out = out.replace("T","3 ")
	out = out.replace("a","0 ")
	out = out.replace("c","1 ")
	out = out.replace("g","2 ")
	out = out.replace("t","3 ")
	out = out.replace("N","4 ")
	out = out.replace("X","5 ")
	return(out)

fasta = sys.argv[1]
info = sys.argv[2]
pssm = sys.argv[3]

outpath = sys.argv[4]

pssmArray = [0] * 100
signArray = ["+"] * 100
cnt = 0
totalLength = 0
toreverse = []
fin = open(pssm)
for line in fin:
	if not re.search(" total sequences", line): continue
	words = line.split()
	pssmArray[int(words[1][1:-1])-1] = len(words[5])
	totalLength += len(words[5])
	cnt = cnt + 1
fin.close()

fin = open(info)
ffa = open(fasta)
fout1 = open(outpath + "/sitesData.txt","w")

fastal = ffa.readline()
fastal = ffa.readline()
for line in fin:
	line = line.strip()
	words = line.split("\t")
	fasta = ''
        while fastal != '' and fastal[0] != '>':
                fasta += fastal.strip()
	        fastal = ffa.readline()
        fastal = ffa.readline()

	out = ""
	for i in range(3,len(words)):
		c = i - 3
		if pssmArray[c] == 0: continue
			
		if(words[i] == "NA"):
			out = out + "N" *  pssmArray[c]
		else:
#			if signArray[c] == "-":
			if c in toreverse:
				words[i] = -1 * int(words[i])
			if(int(words[i]) < 0):
                                #print fasta
				out = out + reverseComp(fasta[abs(int(words[i])):abs(int(words[i]))+pssmArray[c]])
			else:
				out = out + fasta[abs(int(words[i])):abs(int(words[i]))+pssmArray[c]]
		out = out + "X" 
	fout1.write(convertToNumbers(out[:-2]) + "\n")
	fout1.flush()


fin.close()
ffa.close()
fout1.close()
