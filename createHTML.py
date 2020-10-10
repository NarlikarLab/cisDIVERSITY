from HTMLStrings import *
import numpy as np
def writeAllMotifs(seq,fp):
    k=0
    l=len(seq)
    while (k<l and seq[k]=="NA"):
        k=k+1
    if k==l:
        return
    else:
        motifpos = abs(int(seq[k]))
        if int(seq[k])>0:
            direc="R"
        else:
            direc="L"
        motifid = k+1
        fp.write("\t"+str(motifid)+":{\"pos\":"+str(motifpos)+",\"dir\":\""+direc+"\"}")
        k=k+1
        for j in range(k,l):
            if (seq[j]=="NA"):
                pass
            else:
                motifpos = abs(int(seq[j]))
                if int(seq[j])>0:
                    direc="R"
                else:
                    direc="L"
                motifid = j+1
                fp.write(", "+str(motifid)+":{\"pos\":"+str(motifpos)+",\"dir\":\""+direc+"\"}")

                
def writejsVars(noOfSeqs,noOfMotifs,noOfModules,minmotpercent,moduleMotifColors,motifwidths,noOfSeqsPerModule,moduleMotifPercent,seqlens,seqheaders,motifSeqs,allmodules,allsizes,files):
    solutionDir = files[0]
    jsfile = files[1]
    fastafile = files[2]
    jsf = open(jsfile,'w')
    
    #wholedata
    jsf.write("var solutionDir=\""+solutionDir+"\";\n")
    jsf.write("var wholedata={\n\"totalseqs\":"+str(noOfSeqs)+",\n\"motifs\":"+str(noOfMotifs)+",\n\"modules\":"+str(noOfModules)+"\n};")
    #inputfasta
    jsf.write("\nvar inputfasta=\""+fastafile+"\";\n")
    
    #maxlength
    jsf.write("\nvar maxlength="+str(max(seqlens.values()))+";")

    #maxlengthSeqTgas
    jsf.write("\nvar maxlengthSeqTags="+str(max(map(len,seqheaders)))+";")

    #Minimum motif percent (threshold)
    jsf.write("\nvar minmotpercent="+str(minmotpercent)+";")

    #logoSizes
    logoSizes = allsizes[0]
    jsf.write("\nvar logoSizes={ ")
    for i in range(noOfMotifs-1):
        jsf.write(str(i+1)+": ["+str(logoSizes[i][0])+" ,"+str(logoSizes[i][1])+"],\n")
    i=i+1
    jsf.write(str(i+1)+": ["+str(logoSizes[i][0])+" ,"+str(logoSizes[i][1])+"]\n};")

    #cplot
    cplot_size = allsizes[1]
    jsf.write("\nvar cplot= ["+str(cplot_size[0])+","+str(cplot_size[1])+"];")

    #fullpart
    fpart_size = allsizes[2]
    jsf.write("\nvar fullpart= ["+str(fpart_size[0])+","+str(fpart_size[1])+"];\n")

    #moduleMotifColors
    modules = allmodules.keys()
    jsf.write("\n\nvar moduleMotifColors={ ")
    for i in range(noOfModules):
        jsf.write(str(modules[i])+": [")
        for j in range(noOfMotifs-1):
            jsf.write("\""+moduleMotifColors[modules[i]][j]+"\", ")
        if i == (noOfModules-1):
            jsf.write("\""+moduleMotifColors[modules[i]][noOfMotifs-1]+"\"]\n")
        else:
            jsf.write("\""+moduleMotifColors[modules[i]][noOfMotifs-1]+"\"],\n")
    jsf.write("};")


    #motifwidths
    jsf.write("\n\nvar motifwidths={\n")
    for k in range(noOfMotifs):
        if (k == noOfMotifs-1):
            jsf.write(str(k+1)+":"+str(motifwidths[k])+"\n};")
        else:
            jsf.write(str(k+1)+":"+str(motifwidths[k])+",\n")

    #modseqs
    jsf.write("\n\nvar modseqs={\n")
    modules = sorted(noOfSeqsPerModule.keys())
    for i in range(len(modules)-1):
        jsf.write(str(modules[i])+": "+str(noOfSeqsPerModule[modules[i]])+",\n")
    i = len(modules)-1
    jsf.write(str(modules[i])+": "+str(noOfSeqsPerModule[modules[i]])+"\n};")
    
    #motifSeqs
    jsf.write("\n\nvar motifSeqs={\n")
    motifs = sorted(motifSeqs.keys())
    for i in range(len(motifs)-1):
        k = motifs[i]
        jsf.write(str(k)+": "+str(motifSeqs[k])+",\n")
    i =len(motifs)-1
    k = motifs[i]
    jsf.write(str(k)+": "+str(motifSeqs[k])+"\n};")
    
    #seqdetails seqheader:[seqno,length]
    jsf.write("\n\nvar seqdetails={\n")
    for k in range(len(seqheaders)-1):
        jsf.write("\""+seqheaders[k]+"\": ["+str(k)+", "+str(seqlens[seqheaders[k]])+"],\n" )
    k = len(seqheaders)-1
    jsf.write("\""+seqheaders[k]+"\": ["+str(k)+", "+str(seqlens[seqheaders[k]])+"]\n};")
    
    #allmodules
    jsf.write("\nvar allmodules={\n")
    mods = sorted(allmodules.keys())
    for k1 in range(len(mods)): 
        k = mods[k1] #k is the module
        jsf.write(str(k)+":{\n")
        seqs =allmodules[k].keys()
        if len(seqs)>1:
            for i in range(len(seqs)-1):
                header = seqs[i]
                motifseq = allmodules[k][header]
                jsf.write("\t\""+header+"\":{\n")
                writeAllMotifs(motifseq,jsf)
                jsf.write("},\n")
        else:
            i=-1
        header = seqs[-1]
        jsf.write("\t\""+header+"\":{\n")
        motifseq = allmodules[k][header]
        writeAllMotifs(motifseq,jsf)
        if k1 == len(mods)-1:
            jsf.write("}\n\t}\n") #closes all sequence headers and this module
        else:
            jsf.write("}\n\t},\n")
    jsf.write("};\n") # closes allmodules
    
    #moduleMotifPercent
    jsf.write("\nvar moduleMotifPercent={\n")
    modules = moduleMotifPercent.keys()
    for i in range(len(modules)):
        jsf.write(str(modules[i])+ ":[ ")
        for j in range(noOfMotifs-1):
            jsf.write(str(moduleMotifPercent[modules[i]][j])+", ")
        j = noOfMotifs-1
        if i==len(modules)-1:
            jsf.write(str(moduleMotifPercent[modules[i]][j])+"]\n")
        else:
            jsf.write(str(moduleMotifPercent[modules[i]][j])+"],\n")
    jsf.write("};\n")
    
    jsf.close()
    
def makeModuleHTML(noOfSeqsPerModule,outdir):            
   for mod in noOfSeqsPerModule.keys():
       modhtml = outdir+"/Module_"+str(mod)+".html"
       mh = open(modhtml,'w')
       mh.write(htmlContent%{'MODNO':str(mod),'MODID':'Module_'+str(mod)})

       #Now write the table
       mh.write(endContent%{'MODID':'Module_'+str(mod),'MODNO':str(mod)})
       mh.close()
    
            
