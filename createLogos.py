import sys
import re
import os
import commands
weblogo = "weblogo/seqlogo  -F PDF -k 1 -c  -Y -h 4 -M "
files = os.listdir(sys.argv[1])
w = sys.argv[0].split("/")
weblogo = "/".join(w[0:-1]) + "/" + weblogo
base = sys.argv[1]
outbase = sys.argv[2]
for f in files:
	if re.match("sites_", f):
		if f.endswith(".txt"):
			out = outbase + "/" + f[:-4] 
			inp = base + "/" + f
			(a,b) = commands.getstatusoutput("wc  " + inp)
			l = b.split()[0]
			if l== "0": continue
			w = str(int(b.split()[2])/int(l))
			mot = int(f[:-4].split("_")[1])
			cmd = weblogo + "-w " + w +  " -f " + base + "/" + f + " -o " +  out + " -t " + "\"motif " + str(mot) +" ("   + l + ")\"  2>/dev/null"
#			print cmd
			os.system(cmd)
	if re.match("revsites_", f):
		if f.endswith(".txt"):
			out = outbase + "/" + f[:-4] 
			inp = base + "/" + f
			(a,b) = commands.getstatusoutput("wc  " + inp)
			l = b.split()[0]
			if l== "0": continue
			w = str(int(b.split()[2])/int(l))
			mot = int(f[:-4].split("_")[1])  
#			cmd = weblogo + "-w " + w +   " -f " + base + "/" + f + " -o " +  out +" -t " + "\"" + l + " sites\"  2>/dev/null"
			cmd = weblogo + "-w " + w +  " -f " + base + "/" + f + " -o " +  out + " -t " + "\"motif " + str(mot) +" ("   + l + ")\"  2>/dev/null"

#			print cmd
			os.system(cmd)
