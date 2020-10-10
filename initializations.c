#include "model.h"
#include <getopt.h>


void setAllConst(structForAllConstants *allCons)
{
    allCons -> pseudo = PSEUDO;
    allCons -> pseudoGamma = PSEUDOGAMMA;
    allCons -> pseudoP = PSEUDOP;
    allCons -> pseudoA = PSEUDOA;
    allCons -> iter = ITER;
    allCons -> improve = IMPROVE;
    allCons -> zprob = ZPROB;
    allCons -> trials = TRIALS;
    allCons -> minw = MINW;
    allCons -> maxw = MAXW;
    allCons -> initw = INITW;
    allCons -> minsites = MINSITES;
    allCons -> minseqs = MINSEQS;
    allCons -> rev = 1;
    allCons -> backorder = BACKORDER;
    allCons -> r = REG;
    allCons -> m = MOTIFS;
    allCons -> verbose = 1;
    allCons -> seed = SEED;
    allCons -> entropy = 0;
    allCons -> threads = NUM_THREADS;
}

void computePartZ(modelStruct *eModel,int mot, structForAllConstants *allCons)
{
    int i,j;
    double sum,t;
    eModel -> partZ[mot] = 0;
    if (fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return;
    for(i = 0; i< eModel ->motifs[mot]->w; i++)
    {
	sum =0;
	for(j = 0; j< 4; j ++)
	{
	    t = eModel -> motifs[mot] -> counts[j][i] ;
	    if(t < 0) t = 0;
	    sum = sum + t;
	    eModel ->partZ[mot] +=  (t - allCons -> pseudo) * logWithCheck(t,1);
	}
	eModel ->partZ[mot] -=  (sum - allCons -> pseudo * 4) * logWithCheck(sum,1);
    }
} 


double **getBackground(FASTA *X, structForAllConstants *allCons)
{
    double **backs;
    int i,j, start,k;
    long int ind;
    MMStruct *M;
    wordsStruct *nWords;
    nWords = ALLOC(sizeof(wordsStruct));
    nWords->order = allCons -> backorder;
    nWords->wordCounts = ALLOC(sizeof(long int*) * (allCons -> backorder+1));
    for(i=0; i<= allCons -> backorder; i++)
    {
	nWords->wordCounts[i] = ALLOC(sizeof(long int) * pow(4,i+1));
    }
    M = ALLOC(sizeof(MMStruct));
    M->order = allCons -> backorder;
    M->probs = ALLOC(sizeof(double*) * (M->order+1));
    for(i=0; i<=M->order; i++)
    {
	M->probs[i] = ALLOC(sizeof(double) * pow(4,i+1));
    }
    backs = ALLOC(sizeof(double *) * X -> num);
    for(j=0;j<X -> num; j++)
    {
	backs[j] = ALLOC(sizeof(double) * X -> length[j]);
	for(i=0; i<=nWords->order; i++)
	{
	    memset(nWords->wordCounts[i], 0 , sizeof(long int) * pow(4,i+1));
	    memset(M -> probs[i], 0 , sizeof(long int) * pow(4,i+1));
	}
	for(i=0; i<X -> length[j] ; i++)
	{
	    if(X -> seq[j][i] == 4)
		continue;
	    for(k=0;k<=nWords->order ;k++)
	    {
		if( k + i >= X -> length[j]) continue;
		ind = getIndexInt(X -> seq[j], i,k);
		if(ind == -1) continue;
		nWords->wordCounts[k][ind]++;
		M -> probs[k][ind]++;
	    }
	}
	normalizeMM(M);
	for(i=0; i<X -> length[j] ; i++)
	{
	    if(X -> seq[j][i] == 4)
	    {
		backs[j][i] = -1;
		continue;
	    }
	    start = maxPieceWithoutN(X -> seq[j], i,M->order );
	    
	    ind = getIndexInt(X -> seq[j], start, i-start );
	    if(ind == -1)
		backs[j][i] = -1;
	    else
		backs[j][i] = pow(M->probs[i-start][ind],1);
	} 
    }
    freeWords(nWords);
    freeMM(M);
    return(backs);
}



FASTA *readFASTA(char *filename,structForAllConstants *allCons )
{
    FASTA *sequences;
    FILE *fp;
    int maxlength = 0, length = 0, i, nonDNAChars = 0, back,front,j;
    char c;
    char buffer[HEADERLIMIT];
    if( NULL == (fp = fopen(filename, "r")))
    {
	fileError(filename, " is not readable");
    }

    sequences = ALLOC(sizeof(FASTA));
    sequences -> num = 0;
    sequences -> revFlag = allCons -> rev;
    while(!feof(fp))
    {
	fscanf(fp,"%c",&c);
	if(c == '>')
	{
	    if(!(fgets(buffer,HEADERLIMIT,fp)))    
	    {
		fileError(filename, " has too long a header");
	    }
	    if(maxlength < length)			   
			maxlength = length;
		
	    length = 0;
	    sequences -> num++;
	}
	else
	{
	    if(c == '\n')
		continue;
	    length++;
	}
    }
    if(maxlength < length)
		maxlength = length;


    if( sequences -> revFlag == 1)
    {
	maxlength *= 2;
    }
    fprintf(stderr, "Read %s: %d sequences, maximum sequence length %d\n", filename, sequences -> num, maxlength);
			
    maxlength+=3;
    sequences -> seq = ALLOC(sizeof(int *) * sequences -> num);
    sequences -> noN = ALLOC(sizeof(int *) * sequences -> num);
    sequences -> header = ALLOC(sizeof(char *) * sequences -> num);
    sequences -> length = ALLOC(sizeof(int ) * sequences -> num);
    sequences -> nonRepLength = ALLOC(sizeof(int ) * sequences -> num);
    sequences -> tableOfBackground = ALLOC(sizeof(double **)  * sequences -> num);
    
    for(i = 0; i <  sequences -> num; i++)
    {
	sequences -> seq[i] = ALLOC(sizeof(int) * maxlength);
	sequences -> header[i] = ALLOC(sizeof(char) * HEADERLIMIT);
		
    }

    rewind(fp);
    i = -1;
    while(!feof(fp))
    {
	fscanf(fp, "%c", &c);

	switch (c)
	{
	case '>':
	    if(!(fgets(buffer,HEADERLIMIT,fp)))    
	    {
		fileError(filename, " has too long a header");
	    }
	    i++;
	    strtok(buffer,"\n");
	    strcpy(sequences -> header[i], buffer);
	    
	    sequences -> length[i] = 0;
	    sequences -> nonRepLength[i] = 0;
	    break;
	case 'a':
	case 'A':
	    sequences -> seq[i][sequences -> length[i]] = 0;
	    sequences -> length[i]++;
	    sequences -> nonRepLength[i]++;
	    break;
        case 'c':
	case 'C':
	    sequences -> seq[i][sequences -> length[i]] = 1;
	    sequences -> length[i]++;
	    sequences -> nonRepLength[i]++;
	    break;
	case 'g':
	case 'G':
	    sequences -> seq[i][sequences -> length[i]] = 2;
	    sequences -> length[i]++;
	    sequences -> nonRepLength[i]++;
	    break;        
        case 't':
	case 'T':
	    sequences -> seq[i][sequences -> length[i]] = 3;
	    sequences -> length[i]++;
	    sequences -> nonRepLength[i]++;
	    break;
	case '\n':
	    break;		
	default:
	    sequences -> seq[i][sequences -> length[i]] = 4;
	    sequences -> length[i]++;
	    nonDNAChars++;
	    break;
	}
    }
    if(nonDNAChars > 0)
	fprintf(stderr, "Found %d non ACGT characters in %s\n", nonDNAChars, filename);
    fclose(fp);

    if( sequences -> revFlag == 1)
    {
	for(i = 0; i <  sequences -> num; i++)
	{
	    sequences->seq[i][sequences -> length[i]] = 4;
	    back = sequences->length[i] * 2;
	    front = 0;
	    while(back > sequences->length[i])
	    {
		sequences -> seq[i][back] = 3 - sequences -> seq[i][front];
		if(sequences -> seq[i][back] < 0)
		    sequences -> seq[i][back] = 4;
		back--;
		front++;
	    }
	    sequences->length[i] = sequences->length[i] * 2 + 1;
	    sequences -> nonRepLength[i] =  sequences -> nonRepLength[i] * 2 + 1;

	}
    }
    
    for(i = 0; i <  sequences -> num; i++)
    {
	sequences -> noN[i] = ALLOC(sizeof(int ) * sequences->length[i]);
	memset(sequences -> noN[i] ,0,sizeof(int ) * sequences->length[i]);
	setMaxSequenceWithoutN(sequences, i);
	sequences -> tableOfBackground[i] = ALLOC(sizeof(double *)  * sequences -> length[i]);
	for(j=0 ;j < sequences -> length[i]; j++)
	{
	    sequences -> tableOfBackground[i][j] = ALLOC(sizeof(double) * allCons -> maxw);
	    memset(sequences -> tableOfBackground[i][j] ,0,sizeof(double) * allCons -> maxw);
	}
    }
    
    return(sequences);
}

void setBackTable(FASTA *X, double **backs, int maxW)
{
    int i,j,k;
    double ans;
    for(i =0;i< X-> num;i++)
    {
	for(j=0; j< X -> length[i]; j++)
	{
	    ans = 1;
	    for(k=0;k< maxW; k++)
	    {
		if( j + k >= X -> length[i]) break;
		if(backs[i][j+k] == -1 || backs[i][j+k] == 0 )
		{
		    break;
		}
		ans = ans * backs[i][j+k];
		X -> tableOfBackground[i][j][k] = ans;
	    }
	}
    }
}



pssms *readMotifs(char *filename, structForAllConstants *allCons)
{
    int i,j;
    FILE *fp;
    char c;
    pssms *motifs;
    double temp;
    int k,w,tot;
    
    if( NULL == (fp = fopen(filename, "r")))
    {
	fileError(filename, " is not readable");
    }
    fscanf(fp,"%c",&c);
    fscanf(fp,"Number of sequences: %d\n",&i);
    fscanf(fp,"Number of modes: %d\n",&tot);
    motifs = ALLOC(sizeof(pssms ));
    motifs -> number = tot;
    motifs -> motif = ALLOC(sizeof (phiStruct *) * tot);
    fscanf(fp,"Likelihood: %lf\n\n",&temp);
    for(i=0;i<tot;i++)
    {
	fscanf(fp,"\nMOTIF: %d (%d sequences, width %d)\n\n",&j, &k, &w);
	if(i != j)
	{
	    printf("error in reading motif file\n");
	    exit(0);
	}
	motifs -> motif[i] = ALLOC(sizeof (phiStruct));
	motifs -> motif[i] -> w = w;
	motifs -> motif[i] -> phi = ALLOC(sizeof(double *) * 4);
	for(j = 0; j < 4; j++)
	{
	   motifs -> motif[i] -> phi[j]  = ALLOC(sizeof(double) * w);
	}
	for(j=0;j<4;j++)
	{
	    fscanf(fp,"%c [",&c);
	    for(w = 0; w< motifs -> motif[i] -> w-1 ;w ++)
	    {
		fscanf(fp,"%lf\t",&motifs -> motif[i] -> phi[j][w]);
		motifs -> motif[i] -> phi[j][w] = (motifs -> motif[i] -> phi[j][w] * k + allCons -> pseudo)/(k+allCons -> pseudo);
	    }
	    fscanf(fp,"%lf]\n",&motifs -> motif[i] -> phi[j][w]);
	    motifs -> motif[i] -> phi[j][w] = (motifs -> motif[i] -> phi[j][w] * k + allCons -> pseudo)/(k+allCons -> pseudo);
	    
	}
	fscanf(fp,"\n");
    }
    return(motifs); 
}

void justCopyModel(modelStruct *from, modelStruct *to)
{
    to -> totDenom = from -> totDenom ;
    to -> minW = from -> minW ;
    to -> maxW = from -> maxW ;

    copyI(from,to);
    copyz(from,to);
    copyCounts(from,to);

}



modelStruct  *makeACopyOfModel(modelStruct *from)
{
    int i,j;
    modelStruct *to;
    to = ALLOC(sizeof(modelStruct));
    to -> numOfSequences = from -> numOfSequences;
    to -> numOfRegProgs = from -> numOfRegProgs ;
    to -> numOfMotifs = from -> numOfMotifs;
    to -> minW = from -> minW ;
    to -> maxW = from -> maxW ;
    to -> totMotifs = ALLOC(sizeof(int) * to -> numOfSequences);
  
    to -> z = ALLOC(sizeof(int *) * to -> numOfSequences);
    to -> zDenom = ALLOC(sizeof(double *) * to -> numOfSequences);
    for(i = 0; i <  to -> numOfSequences; i++)
    {
	to -> z[i] = ALLOC(sizeof(int) * to -> numOfMotifs);
	to -> zDenom[i] = ALLOC(sizeof(double) * to -> numOfMotifs);
    }
    to -> totDenom = from -> totDenom ;
    
    to -> I = ALLOC(sizeof(int) * to -> numOfSequences);
    memset(to -> I ,0,sizeof(int ) * to -> numOfSequences);

    to -> motifs = ALLOC(sizeof(countStruct *) * to -> numOfMotifs);

     
    for(i = 0; i <  to -> numOfMotifs; i++)
    {
	to -> motifs[i] = ALLOC(sizeof(countStruct));
	to -> motifs[i] -> w = from -> motifs[i] -> w;
	to -> motifs[i] -> counts = ALLOC(sizeof(double *) * 4);
	to -> motifs[i] -> total = 0;
	for(j = 0; j < 4; j++)
	{
	    to -> motifs[i] -> counts[j]  = ALLOC(sizeof(double) * to -> maxW );
	}
		 
    }
    
    to -> gammaCounts = ALLOC(sizeof(double) * to -> numOfRegProgs);
    to -> partZ = ALLOC(sizeof(double) * to -> numOfMotifs);

    to -> VAcounts = ALLOC(sizeof(double *) * to -> numOfRegProgs);
    to -> VPcounts = ALLOC(sizeof(double *) * to -> numOfRegProgs);
    for(j =0 ;j< to -> numOfRegProgs; j++)
    {
	to -> VAcounts[j] = ALLOC(sizeof(double) * to -> numOfMotifs);
	to -> VPcounts[j] = ALLOC(sizeof(double) * to -> numOfMotifs);
	memset( to -> VAcounts[j], 0, sizeof( double) * to -> numOfMotifs);
	memset( to -> VPcounts[j], 0, sizeof( double) * to -> numOfMotifs);
    }
    copyI(from,to);
    copyz(from,to);
    copyCounts(from,to);
    return(to);
}

void copyCounts(modelStruct *from,modelStruct *to)
{
    int i,j,k;
    for(i = 0; i <  to -> numOfMotifs; i++)
    {
	to -> motifs[i] -> w = from -> motifs[i] -> w;
	for(j = 0; j < 4; j++)
	    for(k = 0;k< to -> motifs[i] -> w ;k++)
		to -> motifs[i] -> counts[j][k]  = from -> motifs[i] -> counts[j][k] ;
	to -> motifs[i] -> total = from -> motifs[i] -> total;
	to -> partZ[i] = from -> partZ[i];
    }
}




void copyz(modelStruct *from,modelStruct *to)
{
    int i,k;
    for (i = 0; i< to -> numOfSequences;i++)
    {
	for(k =0 ; k< to -> numOfMotifs; k++)
	{
	    to -> z[i][k] = from ->  z[i][k];
	    to -> zDenom[i][k] =  from -> zDenom[i][k];
	}
	to -> totMotifs[i] = from -> totMotifs[i];
    }
    for (i = 0; i< to -> numOfRegProgs;i++)
    {
	for(k =0 ; k< to -> numOfMotifs; k++)
	{
	    to -> VAcounts[i][k] = from -> VAcounts[i][k]   ;
	    to -> VPcounts[i][k] = from -> VPcounts[i][k]   ;

	}
    }
}

void copyI(modelStruct *from,modelStruct *to)
{
    int i;
    for (i = 0; i< to -> numOfRegProgs;i++)
    {
	to -> gammaCounts[i] = 	from -> gammaCounts[i];
    }
    for (i = 0; i< from -> numOfSequences;i++)
    {
	to -> I[i] = from -> I[i] ;
    }
}

modelStruct *cleanUpModel(modelStruct *eModel, structForAllConstants *allCons)
{
    modelStruct *cleaned;
    int **regI, **motI;
    int i,j,k;
    cleaned = ALLOC(sizeof(modelStruct));
    cleaned  -> numOfSequences  = eModel  -> numOfSequences ;
    cleaned -> numOfRegProgs = eModel -> numOfRegProgs;
    cleaned -> numOfMotifs = eModel  -> numOfMotifs;
    cleaned -> minW = eModel -> minW;
    cleaned -> maxW = eModel -> maxW;
    regI = ALLOC(sizeof(int *) *  eModel -> numOfRegProgs);
    motI = ALLOC(sizeof(int *) *  eModel -> numOfMotifs);
    for(i =0;i< cleaned -> numOfRegProgs;i++)
    {
	regI[i] = ALLOC(sizeof(int) * 3);
	regI[i][1] = i;
	regI[i][0] = eModel -> gammaCounts[i] - allCons -> pseudoGamma;
	regI[i][2] = i;
    }
    for(i =0;i< cleaned -> numOfMotifs;i++)
    {
	motI[i] = ALLOC(sizeof(int) * 3);
	motI[i][1] = i;
	motI[i][2] = i;
	motI[i][0] = eModel -> motifs[i] -> total - 4*  allCons-> pseudo ;
    }
   
    selectionSort(regI, eModel -> numOfRegProgs);
    selectionSort(motI, eModel -> numOfMotifs);
    
    cleaned -> z = ALLOC(sizeof(int *) * eModel -> numOfSequences);
    cleaned -> zDenom = ALLOC(sizeof(double *) * eModel -> numOfSequences);
   
    cleaned -> totDenom = eModel -> totDenom ;
    
    cleaned -> I = ALLOC(sizeof(int) * cleaned -> numOfSequences);
    memset(cleaned -> I ,0,sizeof(int ) * cleaned -> numOfSequences);
 
    cleaned -> gammaCounts = ALLOC(sizeof(double) * cleaned -> numOfRegProgs);
    cleaned -> partZ = ALLOC(sizeof(double) * cleaned -> numOfMotifs);

    for (i = 0; i< cleaned -> numOfRegProgs;i++)
    {
	cleaned -> gammaCounts[i] =  eModel -> gammaCounts[regI[i][1]];
    }
    for (i = 0; i< cleaned -> numOfSequences;i++)
    {
	cleaned -> I[i] = regI[eModel -> I[i]][2] ;
    }
    cleaned -> totMotifs = ALLOC(sizeof(int) * cleaned -> numOfSequences);

    for(i = 0; i <  cleaned -> numOfSequences; i++)
    {
	cleaned -> z[i] = ALLOC(sizeof(int) * cleaned -> numOfMotifs);
	cleaned -> zDenom[i] = ALLOC(sizeof(double) * cleaned -> numOfMotifs);
	for(k =0 ; k< cleaned -> numOfMotifs; k++)
	{
	    cleaned -> z[i][k] = eModel ->  z[i][motI[k][1]];
	    cleaned -> zDenom[i][k] =  eModel -> zDenom[i][motI[k][1]];
	}
	cleaned -> totMotifs[i] = eModel -> totMotifs[i];
    }
    
    cleaned -> VAcounts = ALLOC(sizeof(double *) * cleaned -> numOfRegProgs);
    cleaned -> VPcounts = ALLOC(sizeof(double *) * cleaned -> numOfRegProgs);
    for(j =0 ;j< cleaned -> numOfRegProgs; j++)
    {
	cleaned -> VAcounts[j] = ALLOC(sizeof(double) * cleaned -> numOfMotifs);
	cleaned -> VPcounts[j] = ALLOC(sizeof(double) * cleaned -> numOfMotifs);
	for(k =0 ; k< cleaned -> numOfMotifs; k++)
	{
	    cleaned -> VAcounts[j][k] = eModel -> VAcounts[regI[j][1]][motI[k][1]]   ;
	    cleaned -> VPcounts[j][k] = eModel -> VPcounts[regI[j][1]][motI[k][1]]   ;

	}
    }
    cleaned -> motifs = ALLOC(sizeof(countStruct *) * cleaned -> numOfMotifs);

    for(i = 0; i <  cleaned -> numOfMotifs; i++)
    {
	cleaned -> motifs[i] = ALLOC(sizeof(countStruct));
	cleaned -> motifs[i] -> w = eModel -> motifs[motI[i][1]] -> w;
	cleaned -> motifs[i] -> counts = ALLOC(sizeof(double *) * 4);
	for(j = 0; j < 4; j++)
	{
	    cleaned -> motifs[i] -> counts[j]  = ALLOC(sizeof(double) * cleaned -> maxW );
	    for(k = 0;k< cleaned -> motifs[i] -> w ;k++)
		cleaned -> motifs[i] -> counts[j][k]  = eModel -> motifs[motI[i][1]] -> counts[j][k] ;
	}
	cleaned -> motifs[i] -> total = eModel -> motifs[motI[i][1]] -> total;
	cleaned -> partZ[i] = eModel -> partZ[motI[i][1]];
    		 
    }

    for(i =0;i< cleaned -> numOfRegProgs;i++)
    {
	free(regI[i]);
    }
    for(i =0;i< cleaned -> numOfMotifs;i++)
    {
	free(motI[i]);
    }
    free(regI);
    free(motI);
    return(cleaned);

    
}


modelStruct *initializeRandomModel(FASTA *X,structForAllConstants *allCons,unsigned int *seeds)
{
    modelStruct *eModel;
    int i,j;

    eModel = ALLOC(sizeof(modelStruct));

    eModel -> numOfSequences = X -> num;
    eModel -> numOfRegProgs = allCons -> r;
    eModel -> numOfMotifs = allCons -> m;
    eModel -> minW = allCons -> minw;
    eModel -> maxW = allCons -> maxw;
    eModel -> totMotifs = ALLOC(sizeof(int ) * eModel -> numOfSequences);
    memset(eModel -> totMotifs ,0,sizeof(int ) * eModel -> numOfSequences);
    
    if(eModel -> numOfRegProgs > pow(2,eModel -> numOfMotifs))
    	printError("Number of regulatory programs has to be less than 2^numberOfMotifs\n");
    

    eModel -> z = ALLOC(sizeof(int *) * eModel -> numOfSequences);
    eModel -> zDenom = ALLOC(sizeof(double *) * eModel -> numOfSequences);

    for(i = 0; i <  eModel -> numOfSequences; i++)
    {
	eModel -> z[i] = ALLOC(sizeof(int) * eModel -> numOfMotifs);
	memset(	eModel -> z[i] ,0,sizeof(int ) * eModel -> numOfMotifs);
	eModel -> zDenom[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(	eModel -> zDenom[i] ,0,sizeof(double ) * eModel -> numOfMotifs);
    }

    eModel -> totDenom = 0;
    
    eModel -> I = ALLOC(sizeof(int) * eModel -> numOfSequences);
    memset(eModel -> I ,0,sizeof(int ) * eModel -> numOfSequences);

    eModel -> VAcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    eModel -> VPcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    for(i = 0; i <  eModel -> numOfRegProgs; i++)
    {
	eModel -> VAcounts[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(eModel -> VAcounts[i], 0, eModel -> numOfMotifs);
	eModel -> VPcounts[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(eModel -> VPcounts[i], 0, eModel -> numOfMotifs);
    }
    eModel -> motifs = ALLOC(sizeof(countStruct *) * eModel -> numOfMotifs);

    for(i = 0; i <  eModel -> numOfMotifs; i++)
    {
	eModel -> motifs[i] = ALLOC(sizeof(countStruct));
	eModel -> motifs[i] -> w = allCons -> initw;
	eModel -> motifs[i] -> counts = ALLOC(sizeof(double *) * 4);
	eModel -> motifs[i] -> total = 0;
	for(j = 0; j < 4; j++)
	{
	    eModel -> motifs[i] -> counts[j]  = ALLOC(sizeof(double) * eModel -> maxW );
	    memset( eModel -> motifs[i] -> counts[j], 0, sizeof(double) * eModel -> maxW);
	}
	
    }
    
    
    eModel -> gammaCounts = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset( eModel -> gammaCounts, 0 , sizeof( double ) * eModel -> numOfRegProgs);
    eModel -> partZ = ALLOC(sizeof(double) * eModel -> numOfMotifs);
    memset( eModel -> partZ, 0 , sizeof( double ) * eModel -> numOfMotifs);

    initializeRandomI(eModel,allCons,seeds);
    initializeRandomz(X,eModel,allCons,seeds);
    computeCounts(X,eModel,allCons);
    return(eModel);
}




modelStruct *initializeNullModel(FASTA *X,structForAllConstants *allCons)
{
    modelStruct *eModel;
    int i,j,k;

    eModel = ALLOC(sizeof(modelStruct));

    eModel -> numOfSequences = X -> num;
    eModel -> numOfRegProgs = allCons -> r;
    eModel -> numOfMotifs = allCons -> m;
    eModel -> minW = allCons -> minw;
    eModel -> maxW = allCons -> maxw;
    eModel -> totMotifs = ALLOC(sizeof(int ) * eModel -> numOfSequences);
    memset(eModel -> totMotifs ,0,sizeof(int ) * eModel -> numOfSequences);
    
    if(eModel -> numOfRegProgs > pow(2,eModel -> numOfMotifs))
    	printError("Number of regulatory programs has to be less than 2^numberOfMotifs\n");
    

    eModel -> z = ALLOC(sizeof(int *) * eModel -> numOfSequences);
    eModel -> zDenom = ALLOC(sizeof(double *) * eModel -> numOfSequences);

    for(i = 0; i <  eModel -> numOfSequences; i++)
    {
	eModel -> z[i] = ALLOC(sizeof(int) * eModel -> numOfMotifs);
	memset(	eModel -> z[i] ,0,sizeof(int ) * eModel -> numOfMotifs);
	eModel -> zDenom[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(	eModel -> zDenom[i] ,0,sizeof(double ) * eModel -> numOfMotifs);
    }

    eModel -> totDenom = 0;
    
    eModel -> I = ALLOC(sizeof(int) * eModel -> numOfSequences);
    memset(eModel -> I ,0,sizeof(int ) * eModel -> numOfSequences);

    eModel -> VAcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    eModel -> VPcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    for(i = 0; i <  eModel -> numOfRegProgs; i++)
    {
	eModel -> VAcounts[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(eModel -> VAcounts[i], 0, eModel -> numOfMotifs);
	eModel -> VPcounts[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(eModel -> VPcounts[i], 0, eModel -> numOfMotifs);
    }
    eModel -> motifs = ALLOC(sizeof(countStruct *) * eModel -> numOfMotifs);

    for(i = 0; i <  eModel -> numOfMotifs; i++)
    {
	eModel -> motifs[i] = ALLOC(sizeof(countStruct));
	eModel -> motifs[i] -> w = allCons -> initw;
	eModel -> motifs[i] -> counts = ALLOC(sizeof(double *) * 4);
	eModel -> motifs[i] -> total = 0;
	for(j = 0; j < 4; j++)
	{
	    eModel -> motifs[i] -> counts[j]  = ALLOC(sizeof(double) * eModel -> maxW );
	    memset( eModel -> motifs[i] -> counts[j], 0, sizeof(double) * eModel -> maxW);
	}
	
    }
    
    
    eModel -> gammaCounts = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset( eModel -> gammaCounts, 0 , sizeof( double ) * eModel -> numOfRegProgs);
    eModel -> partZ = ALLOC(sizeof(double) * eModel -> numOfMotifs);
    memset( eModel -> partZ, 0 , sizeof( double ) * eModel -> numOfMotifs);
    //empty I
    for (i = 0; i< eModel -> numOfRegProgs;i++)
    {
	eModel -> gammaCounts[i] =  allCons -> pseudoGamma;
    }
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	eModel -> I[i] =  0;
	eModel -> gammaCounts[eModel -> I[i] ]++;
    }
    //empty z
    for(i =1; i< eModel -> numOfRegProgs; i++)
	for(k=0; k< eModel -> numOfMotifs; k++)
	{
	    eModel -> VAcounts[i][k] = allCons -> pseudoA;
	    eModel -> VPcounts[i][k] = allCons -> pseudoP;
	}
    for(k=0; k< eModel -> numOfMotifs; k++)
    {
	eModel -> VAcounts[0][k] = allCons -> pseudoA + eModel -> numOfSequences;
	eModel -> VPcounts[0][k] = allCons -> pseudoP;
    }
    
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	for(k =0 ; k< eModel -> numOfMotifs; k++)
	{
	    eModel -> z[i][k] = -1;
	}
    }
    
    computeCounts(X,eModel,allCons);
    return(eModel);
}


void computeCounts(FASTA *X, modelStruct *eModel, structForAllConstants *allCons)
{
    int i,j,k;
    double t;
    memset( eModel -> partZ, 0 , sizeof( double ) * eModel -> numOfMotifs);
    
    for(k = 0; k < eModel->numOfMotifs; k++)
    {
	for(i=0;i < eModel -> motifs[k] -> w;i++)
	{
	    for(j=0;j<4;j++)
		eModel -> motifs[k] -> counts[j][i] = allCons -> pseudo;
	}
	eModel -> motifs[k] -> total =  allCons -> pseudo * 4;

    }
    for(i=0;i < eModel -> numOfSequences; i++)
	for(k=0;k<eModel ->numOfMotifs;k++)
	{
	    if (eModel -> z[i][k] == -1) continue;
	    for(j=0; j < eModel ->motifs[k] -> w; j++)
		eModel -> motifs[k] -> counts[X -> seq[i][ eModel -> z[i][k] + j  ]][j] ++;
	    eModel -> motifs[k] -> total ++;
	}
    for(k = 0; k < eModel -> numOfMotifs; k++)
    {
	if (fabs(eModel -> motifs[k] -> total - allCons -> pseudo * 4) < EPS) continue;
	for(i = 0; i< eModel ->motifs[k]->w; i++)
	{
	    for(j = 0; j< 4; j ++)
	    {
		t = eModel -> motifs[k] -> counts[j][i];
		eModel ->partZ[k] +=  (t - allCons -> pseudo) * logWithCheck(t,1);
	    }
	    eModel ->partZ[k] -=  (eModel -> motifs[k] -> total - allCons -> pseudo * 4) * logWithCheck(eModel -> motifs[k]->total  ,1);
	}
	
    }
}

void computeVAVP(modelStruct *eModel,  structForAllConstants *allCons)
{
    int i,k;
    for(i =0; i< eModel -> numOfRegProgs; i++)
	for(k=0; k< eModel -> numOfMotifs; k++)
	{
	    eModel -> VAcounts[i][k] = allCons -> pseudoA;
	    eModel -> VPcounts[i][k] = allCons -> pseudoP;
	}
    for(i=0;i< eModel->numOfSequences;i++)
    	for(k =0;k <eModel -> numOfMotifs;k++)
	{
	    if(eModel -> z[i][k] == -1)
		eModel -> VAcounts[eModel -> I[i]][k]++;
	    else
	    {
		eModel -> VPcounts[eModel -> I[i]][k]++;
	    }
	}
}
void computeCountsMotif(FASTA *X, modelStruct *eModel, int mot, structForAllConstants *allCons)
{
    int i,j;
    double t;
    for(i=0;i < eModel -> motifs[mot] -> w;i++)
    {
	for(j=0;j<4;j++)
	    eModel -> motifs[mot] -> counts[j][i] = allCons -> pseudo;
    }
    eModel -> motifs[mot] -> total = allCons -> pseudo * 4;

    
    for(i=0;i < eModel -> numOfSequences; i++)
    {
	if (eModel -> z[i][mot] == -1) continue;
	for(j=0; j < eModel ->motifs[mot] -> w; j++)
	    eModel -> motifs[mot] -> counts[X -> seq[i][ eModel -> z[i][mot] + j  ]][j] ++;
	eModel -> motifs[mot] -> total ++;
    }
    eModel -> partZ[mot] = 0;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return;

    for(i = 0; i< eModel ->motifs[mot]->w; i++)
    {
	for(j = 0; j< 4; j ++)
	{
	    t = eModel -> motifs[mot] -> counts[j][i];
	    eModel ->partZ[mot] +=  (t -  allCons -> pseudo) * logWithCheck(t,1);
	}
	eModel ->partZ[mot] -=  (eModel -> motifs[mot]->total - allCons -> pseudo*4) * logWithCheck(eModel -> motifs[mot]->total  ,1);
    }
	
    
}

void initializeRandomIspecial(modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds, int number)
{
    
    int i;
    for (i = 0; i< eModel -> numOfRegProgs;i++)
    {
	eModel -> gammaCounts[i] =  allCons -> pseudoGamma;
    }
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	eModel -> I[i] =  rand_r(seeds) % number;
	eModel -> gammaCounts[eModel -> I[i] ]++;
    }
}

void initializeRandomI(modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds)
{
    int i;
    for (i = 0; i< eModel -> numOfRegProgs;i++)
    {
	eModel -> gammaCounts[i] =  allCons -> pseudoGamma;
    }
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	eModel -> I[i] =  rand_r(seeds) % eModel -> numOfRegProgs;
	eModel -> gammaCounts[eModel -> I[i] ]++;
    }
}



void initializeRandomz(FASTA *X, modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds)
{
    int i,k;
    int mot;
    double R;
    int *givenMotifs;
    int newk;
    for(i =0; i< eModel -> numOfRegProgs; i++)
	for(k=0; k< eModel -> numOfMotifs; k++)
	{
	    eModel -> VAcounts[i][k] = allCons -> pseudoA;
	    eModel -> VPcounts[i][k] = allCons -> pseudoP;
	}
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	for(k =0 ; k< eModel -> numOfMotifs; k++)
	{
	    eModel -> z[i][k] = -1;
	}
    }
    givenMotifs = ALLOC(sizeof (int) * eModel -> numOfMotifs);
    for(k = 0; k < eModel -> numOfMotifs;k++)
    {
	givenMotifs[k] = k;
    }

    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	randomPermute(givenMotifs, eModel -> numOfMotifs,seeds);

	for(newk =0 ; newk< eModel -> numOfMotifs; newk++)
	{
	    k = givenMotifs[newk];
	    R = 0;
	    while(R == 0)
	    {
		R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
	    }
	    if(R <= allCons -> zprob) mot = 1;
	    else mot = 0;
	    if(mot == 0)
	    {
		eModel -> z[i][k] = -1;
	    }
	    else
	    {
		pickAz(X,eModel,i,k,seeds);
		if(eModel -> z[i][k] != -1)
		{
		    eModel -> zDenom[i][k] = scorePieceByTable(eModel,k,i,eModel->z[i][k],X);
		    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][k],1);
		    if(eModel -> totDenom >= 0)
			printf("debug470");
		    
		}
		
	    }
	    if(eModel -> z[i][k] == -1)
		eModel -> VAcounts[eModel -> I[i]][k]++;
	    else
	    {
		eModel -> VPcounts[eModel -> I[i]][k]++;
		eModel -> totMotifs[i]++;
	    }
	}
    }
    free(givenMotifs);
    
}


void pickAz(FASTA *X, modelStruct *eModel, int seq, int mot,unsigned int *seeds)
{
    int flag;
    flag  = 1;
    while(flag < 20)
    {
	eModel -> z[seq][mot] = rand_r(seeds) % (X -> length[seq] - eModel -> motifs[mot] -> w + 1);
	if(pieceContainsNorOverlaps(X,eModel,seq,mot)) flag ++;
	else break;
    }
    if(flag == 20)
    {
	eModel -> z[seq][mot] = -1;
    }
}



int checkConsistency(FASTA *X, modelStruct *eModel, structForAllConstants *allCons)
{
    int i,j,k;
    countStruct **motifs;
    double *gammaCounts;
    int toReturn = 1;
    double **VAcounts;
    double **VPcounts;
    double **zDenom;
    double totDenom = 0;
    int totMotifs = 0;
    double partZ, sum ,t;
    
    motifs = ALLOC(sizeof(countStruct *) * eModel -> numOfMotifs);
    zDenom = ALLOC(sizeof(double *) * eModel -> numOfSequences);
    
    for(i = 0; i <  eModel -> numOfMotifs; i++)
    {
	motifs[i] = ALLOC(sizeof(countStruct));
	motifs[i] -> w = eModel -> motifs[i] -> w;
	if(motifs[i] -> w < eModel -> minW || motifs[i] -> w > eModel -> maxW)
	{
	    printf("motif size %d for motif %d is not within limits.\n", motifs[i] -> w, i );
	    toReturn = 0;
	}
	motifs[i] -> counts = ALLOC(sizeof(double *) * 4);
	motifs[i] -> total = 0;
	for(j = 0; j < 4; j++)
	{
	    motifs[i] -> counts[j]  = ALLOC(sizeof(double) * eModel -> maxW );
	    memset(  motifs[i] -> counts[j], 0, sizeof(double) * eModel -> maxW);
	}	 
    }
    for(i=0;i< eModel -> numOfSequences;i++)
    {
	zDenom[i] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(zDenom[i] ,0,sizeof(double ) * eModel -> numOfMotifs);
	totMotifs = 0;
	for(k=0;k<eModel ->numOfMotifs;k++)
	{
	    if (eModel -> z[i][k] == -1) continue;
	    for(j=0; j < eModel ->motifs[k] -> w; j++)
		motifs[k] -> counts[X -> seq[i][ eModel -> z[i][k] + j  ]][j] ++;
	    motifs[k] -> total ++;
	    totMotifs ++;
	    zDenom[i][k] =  scorePieceByTable(eModel,k,i,eModel->z[i][k],X);
	    if(fabs(zDenom[i][k] - eModel -> zDenom[i][k]) > EPS)
	    {
		printf("zDenom for i %d, k %d is %lf, but recorded as %lf.\n", i, k,zDenom[i][k] , eModel -> zDenom[i][k] );
		toReturn = 0;
	    }
	    totDenom += logWithCheck(eModel -> zDenom[i][k],1);
	}
	if(fabs(totMotifs - eModel -> totMotifs[i]) > EPS)
	{
	    printf("totalMotifs for i %d is %d, but recorded as %d.\n", i, totMotifs,eModel -> totMotifs[i] );
	    toReturn = 0;
		
	}
    }
    if(fabs(totDenom - eModel -> totDenom) > EPS)
    {
	printf("totDenom is %lf, but recorded as %lf, diff is %lf.\n",  eModel -> totDenom,totDenom,totDenom- eModel -> totDenom);
	toReturn = 0;
    }

    
    for(k=0;k<eModel ->numOfMotifs;k++)
    {
	partZ = 0;
	for(j=0; j < eModel ->motifs[k] -> w; j++)
	{
	    sum = 0;
	    for(i=0; i < 4; i++)
	    {
		if(fabs(motifs[k] -> counts[i][j] - (eModel -> motifs[k] -> counts[i][j] - allCons -> pseudo)) > EPS)
		{
		    printf("%d motif counts are not consistent. zs say it should be %lf, but is recorded as %lf.\n", k,eModel -> motifs[k] -> counts[i][j] -  allCons -> pseudo, motifs[k] -> counts[i][j]  );
		    toReturn = 0;
		}
		t = motifs[k] -> counts[i][j]+  allCons -> pseudo;
		if(t < 0) t = 0;
		sum = sum + t;
		partZ +=  motifs[k] -> counts[i][j] * logWithCheck(t,1);
	    }
	    partZ -=  (sum - allCons -> pseudo * 4) * logWithCheck(sum  ,1);
	}
	if(fabs(partZ - eModel -> partZ[k]  ) > EPS)
        {
            printf("%d motif part Z is not consistent. zs say it should be %lf, but is recorded as %lf.\n", k,eModel -> partZ[k], partZ);
            toReturn = 0;
        }
    }
    
    gammaCounts = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset( gammaCounts,0, sizeof(double) * eModel -> numOfRegProgs);
    for(i=0;i< eModel -> numOfSequences;i++)
    {
	gammaCounts[eModel -> I[i]]++;
    }
	
    for(j=0;j<  eModel -> numOfRegProgs; j++)
	if(fabs(gammaCounts[j] - (eModel -> gammaCounts[j] - allCons -> pseudoGamma)) > EPS)
	{
	    printf("%d gamma counts are not consistent. Is say it should be %lf, but is recorded as %lf.\n", j,eModel -> gammaCounts[j]-  allCons -> pseudoGamma  , gammaCounts[j]  );
	    toReturn = 0;
		    
	}
    VAcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    VPcounts = ALLOC(sizeof(double *) * eModel -> numOfRegProgs);
    for(j =0 ;j< eModel -> numOfRegProgs; j++)
    {
	VAcounts[j] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	VPcounts[j] = ALLOC(sizeof(double) * eModel -> numOfMotifs);
	memset(VPcounts[j], 0, sizeof( double) * eModel -> numOfMotifs);
	memset(VAcounts[j], 0, sizeof( double) * eModel -> numOfMotifs);
    }
    for (i = 0; i< eModel -> numOfSequences;i++)
    {
	for(k =0 ; k< eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[i][k] == -1)
		VAcounts[eModel->I[i]][k]++;
	    else 
		VPcounts[eModel->I[i]][k]++;
	}
    }
    for(j =0 ;j< eModel -> numOfRegProgs; j++)
    {
	for(k =0 ; k< eModel -> numOfMotifs; k++)
	{
	    if(fabs(VAcounts[j][k]  - (eModel ->VAcounts[j][k] -  allCons -> pseudoA)) > EPS)
	    {
		printf("%d %d VAcounts are not consistent. Is and Zs say it should be %lf, but is recorded as %lf.\n", j,k,eModel ->VAcounts[j][k] - allCons -> pseudoA  , VAcounts[j][k]   );
		toReturn = 0;

	    }
	    if(fabs(VPcounts[j][k]  - (eModel ->VPcounts[j][k]- allCons -> pseudoP)) > EPS)
	    {
		printf("%d %d VPcounts are not consistent. Is and Zs say it should be %lf, but is recorded as %lf.\n", j,k,eModel ->VPcounts[j][k]- allCons -> pseudoP   , VPcounts[j][k]   );
		toReturn = 0;

	    }


	}
	
    }
    for(k = 0; k< eModel -> numOfMotifs; k++)
	freecountStruct(motifs[k]);
    free(motifs);
    for(j =0 ;j< eModel -> numOfRegProgs; j++)
    {
	free(VAcounts[j]);
	free(VPcounts[j]);
    }
    free(VAcounts);
    free(VPcounts);
    free(gammaCounts);
    for(i=0;i< eModel -> numOfSequences;i++)
    {
	free(zDenom[i]);
	if(pieceContainsNorOverlaps(X,eModel,i,0)==2)
	{
	    printf("seq %d has motifs that overlap\n",i);
	    for(k=0;k< eModel->numOfMotifs;k++)
	    {
		printf("%d at %d;", k, eModel -> z[i][k]); 
	    }
	    toReturn = 0;
	}
    }
    free(zDenom);
    return(toReturn);
}


void cleanUp(modelStruct *eModel, int initW,  structForAllConstants *allCons)
{
    int i,j,k;
    for(k=0; k< eModel -> numOfMotifs; k++)
    {
	if(eModel -> motifs[k] -> total >= allCons -> minsites + allCons -> pseudo * 4) continue;
	if(fabs(eModel -> motifs[k] -> total - allCons -> pseudo * 4) < EPS) continue;
	for(i=0;i<eModel->numOfSequences;i++)
	{
	    if(eModel -> z[i][k] == -1) continue;
	    eModel -> z[i][k] = -1;
	    eModel -> totMotifs[i]--;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][k],1);
	    eModel -> zDenom[i][k] = 0;
	    eModel -> VAcounts[eModel -> I[i]][k]++;
	    eModel -> VPcounts[eModel -> I[i]][k]--;

	}
	eModel -> motifs[k] -> w = initW;
	for(i=0;i < eModel -> motifs[k] -> w;i++)
	{
	    for(j=0;j<4;j++)
		eModel -> motifs[k] -> counts[j][i] = allCons -> pseudo;
	}
	eModel -> motifs[k] -> total = allCons -> pseudo * 4;
	computePartZ(eModel, k, allCons);
    }
}

int removeSmallProgs(modelStruct *eModel, int s,structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j,k;
    double *empty;
    double R,sum = 0,emptyTot = 0;
    int current, new;
    empty = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset(empty, 0, sizeof(double ) * eModel -> numOfRegProgs);
    
    for(j=0;j< eModel -> numOfRegProgs; j++)
    {
	if(eModel -> gammaCounts[j] > s +  allCons -> pseudoGamma)
	{
	    empty[j] = 1;
	    sum++;
	}
	else
	{
	    if(eModel -> gammaCounts[j] ==  allCons -> pseudoGamma)
		emptyTot++;
	}
    }
    if(emptyTot + sum ==  eModel -> numOfRegProgs) return(0);
    if(sum == 0)
    {
	empty[0] = 1;
	sum = 1;
    }
    empty[0] = empty[0] / sum;
    for(j =0; j < eModel -> numOfRegProgs -1 ;j++)
    {
	empty[j+1] = empty[j+1]/sum + empty[j];
    }
    for(i=0;i< eModel -> numOfSequences; i++)
    {
	current = eModel -> I[i];
	if(eModel -> gammaCounts[current] > s +  allCons -> pseudoGamma) continue;
	R = 0;
	while(R == 0)
	{
	    R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
	}
	new = binarySearchSample(R,empty,0,eModel -> numOfRegProgs );
	eModel -> I[i] = new;
	eModel ->gammaCounts[current]--;
	eModel ->gammaCounts[new]++;
	for(k=0;k<eModel->numOfMotifs;k++)
	{
	    if( eModel -> z[i][k] == -1)
	    {
		eModel -> VAcounts[current][k]--;
		eModel -> VAcounts[new][k]++;
	    }
	    else
	    {
		eModel -> VPcounts[current][k]--;
		eModel -> VPcounts[new][k]++;
	    }
	}
			      
    }
    free(empty);
    return(1);
}

void checkConsistencyAndPrint(FASTA *X, modelStruct *eModel, char *errorString, structForAllConstants *allCons)
{
    int consistency = checkConsistency(X, eModel,allCons);
    if(consistency == 0)
    {
	fprintf(stderr,  "%s ", errorString);
	exit(0);
    }
}




int removeSmallProgsCompletely(modelStruct *eModel, int s, structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j,k;
    double *empty, *nonempty;
    double emptyOne = -1;
    int current, new = -1;
    int total = 0;
    double sum = 0;
    double R;
    empty = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    nonempty = ALLOC(sizeof(double) * eModel -> numOfRegProgs);

    memset(empty, 0, sizeof(double ) * eModel -> numOfRegProgs);
    memset(nonempty, 0, sizeof(double ) * eModel -> numOfRegProgs);
    
    for(j=0;j< eModel -> numOfRegProgs; j++)
    {
	if(eModel -> gammaCounts[j] < s +  allCons -> pseudoGamma)
	{
	    empty[j] = 1;
	    new = j;
	    total = total + eModel -> gammaCounts[j] -  allCons -> pseudoGamma;
	}
	else
	{
	    nonempty[j] = 1;
	    sum++;
	}
	if(eModel -> gammaCounts[j] ==  allCons -> pseudoGamma)
	    emptyOne = j;
	    
    }
    
    if(new == -1)
    {
	free(empty);
	return(0);
    }
    if(emptyOne == -1)
	emptyOne = new;
    if(total >= allCons -> minseqs)
    {
	new = emptyOne;
	for(i=0;i< eModel -> numOfSequences; i++)
	{
	    current = eModel -> I[i];
	    if(empty[current] == 0) continue;
	    eModel -> I[i] = new;
	    eModel ->gammaCounts[current]--;
	    eModel ->gammaCounts[new]++;
	    for(k=0;k<eModel->numOfMotifs;k++)
	    {
		if( eModel -> z[i][k] == -1)
		{
		    eModel -> VAcounts[current][k]--;
		    eModel -> VAcounts[new][k]++;
		}
		else
		{
		    eModel -> VPcounts[current][k]--;
		    eModel -> VPcounts[new][k]++;
		}
	    }
	}		      
    }
    else
    {
	if(sum == 0)
	{
	    nonempty[0] = 1;
	    sum = 1;
	}
	nonempty[0] = nonempty[0] / sum;
	for(j =0; j < eModel -> numOfRegProgs -1 ;j++)
	{
	    nonempty[j+1] = nonempty[j+1]/sum + nonempty[j];
	}
	for(i=0;i< eModel -> numOfSequences; i++)
	{
	    current = eModel -> I[i];
	    if(empty[current] == 0) continue;
	    R = 0;
	    while(R == 0)
	    {
		R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
	    }
	    new = binarySearchSample(R,nonempty,0,eModel -> numOfRegProgs );
	    eModel -> I[i] = new;
	    eModel ->gammaCounts[current]--;
	    eModel ->gammaCounts[new]++;
	    for(k=0;k<eModel->numOfMotifs;k++)
	    {
		if( eModel -> z[i][k] == -1)
		{
		    eModel -> VAcounts[current][k]--;
		    eModel -> VAcounts[new][k]++;
		}
		else
		{
		    eModel -> VPcounts[current][k]--;
		    eModel -> VPcounts[new][k]++;
		}
	    }
	}
	
    }

    free(empty);
    free(nonempty);
    return(1);
}
