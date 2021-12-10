'''
Please uncomment the last two statements to copy the 
cisDiversity.html and md.js files in the current directory
'''

import sys,os,math
import numpy as np
import createHTML as ch
import collections,re
import commands

def reqformat(num,noOfdigits):
    noOfdigits = int(noOfdigits)
    if (num==0):
        return "0"*noOfdigits
    d = math.floor(math.log10(num))+1
    counter = d
    if(d<noOfdigits):
        str1 = str(num)
        while (counter < noOfdigits):
            str1 = "0"+str1
            counter +=1
    else:
        str1 = str(num)
    return str1

def distanceBetweenColors(h1,h2):
    def hexs_to_ints(s):
        return [int(s[i:i+2],16) for i in range(1,7,2)]
    return math.sqrt(sum([abs(a-b)**2 for (a,b) in zip(*map(hexs_to_ints,(h1,h2)))]))

def areDifferentColors(colorarr,i):
    target=colorarr[i]
    flag=True
    dist = 100
    for p in range(i):
        if (distanceBetweenColors(colorarr[p],target)>dist):
            continue
        else:
            flag=False
            break
    return flag

def makejsVars(solutionDir,fastafile,infofile,outdir):
    noOfSeqs = 0
    noOfMotifs = 0
    noOfModules = 0
    allmodules = collections.OrderedDict()
    info = open(infofile,'r')
    count = 0
    for il in info:
        count+=1
    info.seek(0)
    noOfSeqs = count
    seqheaders = np.zeros(count,dtype=object)
    seqlens = collections.OrderedDict()
    il = info.readline()
    ill = filter(None,il.strip().split('\t'))
    noOfMotifs = len(ill)-3
    
    motifSeqs = {}
    
    for i in range(1,noOfMotifs+1):
        motifSeqs[i] = 0
    i = 0
    while il:
        ill = il.strip('\n').split('\t')[:-1]
        i = int(ill[0])-1
        seqheaders[i] = ill[1]
        modno = int(ill[2])
        if modno not in allmodules.keys():
            allmodules[modno] = collections.OrderedDict()
            allmodules[modno][seqheaders[i]] = np.zeros(noOfMotifs,dtype=object)
        else:
            allmodules[modno][seqheaders[i]] = np.zeros(noOfMotifs,dtype=object)
        for j in range(noOfMotifs):
            allmodules[modno][seqheaders[i]][j]=ill[j+3]
            if ill[j+3] != 'NA':
                motifSeqs[j+1] +=1
        il = info.readline()
    
    noOfModules = len(allmodules.keys())
    noOfdigits = math.floor(math.log10(noOfMotifs))+1
    ff = open(fastafile,'r')
    fl = ff.readline()
    sl = 0
    header = ''
    while fl:
        if fl[0] == '>':
            if header != '':
                seqlens[header] = sl
            header = fl[1:].strip()
            sl = 0
        else:
            sl += len(fl.strip())
        fl = ff.readline()
    seqlens[header] = sl
    ff.close()
    
    #hypergeomVals = gethypergeomPvals(allmodules)
    
    motifwidths = np.zeros(noOfMotifs,dtype=int)
    for i in range(noOfMotifs):
        sitesfile = solutionDir + '/sites_%02d.txt'%(i+1)
        sitesf = open(sitesfile,'r')
        motifwidths[i]=len(sitesf.readline().strip('\n'))
        sitesf.close()
        
    moduleMotifPercent = {}
    noOfSeqsPerModule = {}
    modules = allmodules.keys()
    for i in range(noOfModules):
        moduleMotifPercent[modules[i]] = np.zeros(noOfMotifs,dtype=float)
        noOfSeqsPerModule[modules[i]]=len(allmodules[modules[i]].keys())
        for chrom in allmodules[modules[i]].keys(): #seqsInThisModule
            for k in range(noOfMotifs):
                if allmodules[modules[i]][chrom][k]=="NA":
                    continue
                else:
                    moduleMotifPercent[modules[i]][k]+=1
        moduleMotifPercent[modules[i]] = moduleMotifPercent[modules[i]]/noOfSeqsPerModule[modules[i]]

    ### Redefining the colors palette
    moduleMotifColors = {}
    for m in modules:
        moduleMotifColors[m] = np.zeros(noOfMotifs,dtype='S7')
    for i in range(noOfModules):
        for j in range(noOfMotifs):
            moduleMotifColors[modules[i]][j] = "#FF0000"
            
    logoSizes = np.zeros((noOfMotifs,2),dtype=int)
    for i in range(noOfMotifs):
        logo = solutionDir+'/sites_%02d.png'%(i+1)
        if not (os.path.isfile(logo)): continue
        (status,out) = commands.getstatusoutput("file "+logo)
        if status!=0:
            print status,out
            exit()
        w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
        logoSizes[i][0] = w/2.5
        logoSizes[i][1] = h/2

    
    circplot = solutionDir + '/circlePlot.png'
    (status,out)=commands.getstatusoutput("file "+circplot)
    if status!=0:
            print status,out
            exit()
    elif "No such file" in out:
            circplot_size = (0,0)
    else:
            w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
            circplot_size = (w,h)

    fullpartition = solutionDir + '/fullPartition.png'
    (status,out)=commands.getstatusoutput("file "+fullpartition)
    if status!=0:
            print status,out
            exit()
    elif "No such file" in out:
            fpart_size = (0,0)
    else:
            w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
            fpart_size = (w,h)

    allsizes = (logoSizes,circplot_size,fpart_size)
    
    minmotpercent = 0.05;

    jsfile = outdir+'/data.js'
    fastafile = os.path.basename(fastafile)
    files = (solutionDir,jsfile,fastafile)
    ch.writejsVars(noOfSeqs,noOfMotifs,noOfModules,minmotpercent,moduleMotifColors,motifwidths,noOfSeqsPerModule,moduleMotifPercent,seqlens,seqheaders,motifSeqs,allmodules,allsizes,files)
    ch.makeModuleHTML(noOfSeqsPerModule,outdir)
            
def main():
    solutionDir = sys.argv[1]
    fastafile = sys.argv[2]
    basepath = sys.argv[3]
    outdir = solutionDir
    infofile = solutionDir+'/info.txt'

    # check if there are no motifs in the output. Then copy different files
    info = open(infofile)
    il = filter(None,info.readline().strip().split('\t'))
    nm = len(il)-3
    info.close()
    if nm == 0:
        os.system("cp "+basepath+"/cisDiversity_nomotif.html "+outdir+"/cisDiversity.html")
    else:
        os.system("cp " + basepath +"/md.js "+outdir)
        os.system("cp " + basepath +"/cisDiversity.html "+outdir)
        makejsVars(solutionDir,fastafile,infofile,outdir)
main()
