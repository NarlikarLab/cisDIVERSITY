'''
Please uncomment the last two statements to copy the 
moduleDiversity.html and md.js files in the current directory
'''

import sys,os,math
import numpy as np
import createHTML as ch
import collections,re
import commands
#from scipy.stats import hypergeom


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

# def gethypergeomPvals(allmodules):
#     modules = allmodules.keys()
#     # count all the instances of motif present
#     mods= allmodules.keys()
#     seqs=allmodules[mods[0]].keys()
#     noOfmotifs=len(allmodules[mods[0]][seqs[0]])
#     #moduleWiseVals = np.zeros((len(modules),noOfmotifs),dtype=float)
#     moduleWiseVals = {}
#     #totSeqsModule = np.zeros(len(modules),dtype=int)
#     totSeqsModule = {}
    
#     for m in modules:
#         totSeqsModule[m] = 0
#         moduleWiseVals[m] = np.zeros(noOfmotifs,dtype=float)
    
#     totMotOccurrences = np.zeros(noOfmotifs,dtype=int)
#     totSeqs=0
#     for m in modules:
#         for s in allmodules[m]:
#             totSeqsModule[m]+=1
#             for i in range(noOfmotifs):
#                 if (allmodules[m][s][i] !="NA"):
#                     moduleWiseVals[m][i]+=1
#                     totMotOccurrences[i]+=1
#     totSeqs=np.sum(totSeqsModule.values())
#     hypergeomVals = {}
#     for m in modules:
#         hypergeomVals[m] = np.zeros(noOfmotifs,dtype=float)
#     for m in modules:
#         #print "module ",m
#         for i in range(noOfmotifs):
#             M=totSeqs                 # Total number of objects in pop
#             n=totMotOccurrences[i]    # Total number of Motif_i among all the objects in pop
#             N=totSeqsModule[m]        # Total number of objects in sample
#             x=moduleWiseVals[m][i]    # Total number of Motif_i among all the objects in sample
            
#             hypergeomVals[m][i] = 1- hypergeom.cdf(x-1,M,n,N) # prob of getting x or more objects in sample
#             #print x,M,n,N,hypergeomVals[m][i]
        
#     del totMotOccurrences,totSeqsModule,moduleWiseVals
#     return hypergeomVals

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
    ill = il.strip('\n').split('\t')[:-1]
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
    flag = 0
    while fl:
        if fl[0]=='>':
            header = fl[1:-1]
            flag =1
            fl = ff.readline()
        elif flag==1:
            sl = len(fl.strip('\n'))
            seqlens[header]= sl
            flag=0
            fl = ff.readline()
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
            # if (hypergeomVals[modules[i]][j] < 0.00001):
            #     moduleMotifColors[modules[i]][j] = "#FF0000"
            # else:
            #     moduleMotifColors[modules[i]][j] = "#008000"
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
    # for i in range(noOfMotifs):
    #     h = 75
    #     w = motifwidths[i]*16.4
    #     logoSizes[i][0] = w
    #     logoSizes[i][1] = h

    
    circplot = solutionDir + '/circlePlot.png'
    (status,out)=commands.getstatusoutput("file "+circplot)
    if status!=0:
            print status,out
            exit()
    w,h = map(int,re.findall('(\d+)\s*x\s*(\d+)',out)[-1])
    circplot_size = (w,h)

    fullpartition = solutionDir + '/fullPartition.png'
    (status,out)=commands.getstatusoutput("file "+fullpartition)
    if status!=0:
            print status,out
            exit()
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
    #outdir = sys.argv[3]
    outdir = solutionDir
    infofile = solutionDir+'/info.txt'
    os.system("cp " + basepath +"/md.js "+outdir)
    os.system("cp " + basepath +"/cisDiversity.html "+outdir)
    makejsVars(solutionDir,fastafile,infofile,outdir)
main()
