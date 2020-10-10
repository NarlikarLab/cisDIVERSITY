

#include "model.h"
#include <getopt.h>



double scorePieceByMotifAO(FASTA *X, modelStruct *eModel, int mot,int seq)
{
    double ans = 1;
    int i;
    int w = eModel -> motifs[mot] -> w;
    int start = eModel->z[seq][mot];
    if(start < 0 ) return 0;
    if(X -> noN[seq][start] < w) return 0;
    if(start +  w > X -> length[seq]) return 0;
    
    for(i=start; i< start +w ; i++)
    {
	ans = ans * eModel -> motifs[mot] -> counts[X -> seq[seq][i]][i-start];
    }
    return(ans);
}


double scorePieceByMotif(FASTA *X, modelStruct *eModel, int mot,int seq)
{
    double ans = 1;
    int i;
    int w = eModel -> motifs[mot] -> w;
    int start = eModel->z[seq][mot];
    double total = eModel -> motifs[mot] -> total;
    if(pieceContainsNorOverlaps(X,eModel,seq, mot)) return 0;
    if(start < 0 ) return 0;
    if(start +  w > X -> length[seq]) return 0;
    
    for(i=start; i< start +w ; i++)
    {
	ans = ans * eModel -> motifs[mot] -> counts[X -> seq[seq][i]][i-start] / total;
    }
    return(ans);
}

double scorePieceByInitializationMotif(FASTA *X, modelStruct *eModel, int mot, phiStruct *motif, int seq, int pos)
{
    double ans = 1;
    int i;
    int w = motif -> w;
    int start = pos;
    int old;
    old = eModel -> z[seq][mot];
    eModel -> z[seq][mot] = start;
    if(pieceContainsNorOverlaps(X,eModel,seq, mot)) return 0;
    if(start < 0 ) return 0;
    if(start +  w > X -> length[seq]) return 0;
    
    for(i=start; i< start +w ; i++)
    {
	ans = ans * motif -> phi[X -> seq[seq][i]][i-start];

    }
    eModel -> z[seq][mot]  = old;
    return(ans);
}


void removeMotifCounts(FASTA *X,modelStruct *eModel, int mot, int seq, structForAllConstants *allCons)
{
    int i, start,end,j;
    double t;
    start = eModel -> z[seq][mot];
    if (start == -1) return;
    end = eModel -> z[seq][mot] + eModel -> motifs[mot] -> w;
    for(i =start ; i<end ;i++)
    {
	eModel -> motifs[mot] -> counts[X->seq[seq][i]][i-start]--;
    }
    eModel -> motifs[mot] -> total--;
    eModel -> partZ[mot] = 0;
    if (fabs(eModel -> motifs[mot] -> total - allCons -> pseudoP * 4) < EPS) return;
    for(i = 0; i< eModel ->motifs[mot]->w; i++)
    {
	for(j = 0; j< 4; j ++)
	{
	    t = eModel -> motifs[mot] -> counts[j][i];
	    eModel ->partZ[mot] +=  (t - allCons -> pseudoP) * logWithCheck(t,1);
	}
	eModel ->partZ[mot] -=  (eModel -> motifs[mot]->total - allCons -> pseudoP * 4)   * logWithCheck(eModel -> motifs[mot]->total  ,1);
    }
    
}

void addMotifCounts(FASTA *X,modelStruct *eModel, int mot, int seq,structForAllConstants *allCons)
{
    int i, start,end,j;
    double t;
    start = eModel -> z[seq][mot];
    end = eModel -> z[seq][mot] + eModel -> motifs[mot] -> w;
    for(i =start ; i<end ;i++)
    {
	eModel -> motifs[mot] -> counts[X->seq[seq][i]][i-start]++;
    }
    eModel -> motifs[mot] -> total++;
    eModel -> partZ[mot] = 0;
    for(i = 0; i< eModel ->motifs[mot]->w; i++)
    {
	for(j = 0; j< 4; j ++)
	{
	    t = eModel -> motifs[mot] -> counts[j][i];
	    eModel ->partZ[mot] +=  (t - allCons -> pseudoP) * logWithCheck(t,1);
	}
	eModel ->partZ[mot] -=  (eModel -> motifs[mot]->total - allCons -> pseudoP * 4)   * logWithCheck(eModel -> motifs[mot]->total  ,1);
    }
    
}

int sampleI(FASTA *X, modelStruct *eModel,  int seq,  double temp,unsigned int *seeds)
{
    int j,k;
    double sum = 0;
    double *weights ;
    double R;
    int new;
    int change;
    change = eModel -> I[seq];
    eModel -> gammaCounts[ eModel -> I[seq]]--;
    eModel -> I[seq] = -1;
    weights = ALLOC(sizeof(double) * eModel -> numOfRegProgs);

    for(k=0;k < eModel -> numOfMotifs; k++)
    {
	if(eModel -> z[seq][k] > -1)
	    eModel -> VPcounts[change][k]--;
	else
	    eModel -> VAcounts[change][k]--;
    } 
    for(j = 0; j< eModel -> numOfRegProgs; j++)
    {
	weights[j] = eModel -> gammaCounts[j];
	for(k=0;k < eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[seq][k] > -1)
		weights[j] = weights[j] * eModel -> VPcounts[j][k];
	    else
		weights[j] = weights[j] *eModel ->  VAcounts[j][k];
	    weights[j] = weights[j]/(eModel -> VPcounts[j][k] + eModel ->  VAcounts[j][k]);

	}
	sum = sum + weights[j];
	if(weights[j] == 0)
	{
	    printf("I weight becoming 0\n");
	    exit(0);
	}
    }
    weights[0] = weights[0]/sum;
    for(j =0; j < eModel -> numOfRegProgs -1 ;j++)
    {
	weights[j+1] = weights[j+1]/sum + weights[j];
    }
    
    R = 0;
    while(R == 0)
    {
	R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
    }
    new = binarySearchSample(R,weights,0,eModel -> numOfRegProgs );
    eModel -> I[seq] = new;
    free(weights);
    eModel -> gammaCounts[ new]++;
    for(k=0;k < eModel -> numOfMotifs; k++)
    {
	if(eModel -> z[seq][k] > -1)
	    eModel -> VPcounts[new][k]++;
	else
	    eModel -> VAcounts[new][k]++;
    }
    if(new == change) return 0;
    else return 1;
}

int sampleIlimitedProgs(FASTA *X, modelStruct *eModel,  int seq,  int regProgs,unsigned int *seeds)
{
    int j,k;
    double sum = 0;
    double *weights ;
    double R;
    int new;
    int change;
    change = eModel -> I[seq];
    eModel -> gammaCounts[ eModel -> I[seq]]--;
    eModel -> I[seq] = -1;
    weights = ALLOC(sizeof(double) * regProgs);

    for(k=0;k < eModel -> numOfMotifs; k++)
    {
	if(eModel -> z[seq][k] > -1)
	    eModel -> VPcounts[change][k]--;
	else
	    eModel -> VAcounts[change][k]--;
    } 
    for(j = 0; j< regProgs; j++)
    {
	weights[j] = eModel -> gammaCounts[j];
	for(k=0;k < eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[seq][k] > -1)
		weights[j] = weights[j] * eModel -> VPcounts[j][k];
	    else
		weights[j] = weights[j] *eModel ->  VAcounts[j][k];
	    weights[j] = weights[j]/(eModel -> VPcounts[j][k] + eModel ->  VAcounts[j][k]);

	}
	sum = sum + weights[j];
	if(weights[j] == 0)
	{
	    printf("I weight becoming 0\n");
	    exit(0);
	}
    }
    weights[0] = weights[0]/sum;
    for(j =0; j < regProgs -1 ;j++)
    {
	weights[j+1] = weights[j+1]/sum + weights[j];
    }
    
    R = 0;
    while(R == 0)
    {
	R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
    }
    new = binarySearchSample(R,weights,0,regProgs );
    eModel -> I[seq] = new;
    free(weights);
    eModel -> gammaCounts[ new]++;
    for(k=0;k < eModel -> numOfMotifs; k++)
    {
	if(eModel -> z[seq][k] > -1)
	    eModel -> VPcounts[new][k]++;
	else
	    eModel -> VAcounts[new][k]++;
    }
    
    if(new == change) return 0;
    else return 1;
}





int pickBestI(FASTA *X, modelStruct *eModel,  int seq,structForAllConstants *allCons)
{
    int j,k;
    double weight ;
    double best ;
    int new;
    int change;
    change = eModel -> I[seq];
  
    for(j = 0; j< eModel -> numOfRegProgs; j++)
    {
	weight = logWithCheck(eModel -> gammaCounts[j]/(allCons -> pseudoGamma * eModel -> numOfRegProgs + eModel -> numOfSequences - 1),1);
	for(k=0;k < eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[seq][k] > -1)
		weight = weight + logWithCheck(eModel -> VPcounts[j][k]/X->nonRepLength[seq],1);
	    else
		weight = weight + logWithCheck(eModel ->  VAcounts[j][k],1);
	    weight = weight - logWithCheck(eModel -> VPcounts[j][k] + eModel ->  VAcounts[j][k],1);

	   
	}
	if(weight == 0)
	{
	    printf("I weight becoming 0\n");
	    exit(0);
	}
	if(j == 0)
	{
	    best = weight;
	    new = j;
	}
	    
	if(weight > best)
	{
	    best = weight;
	    new = j;
	}
    }
    if(change == new) return(0);
    eModel -> I[seq] = new;
    eModel -> gammaCounts[ new]++;
    eModel -> gammaCounts[change]--;
    for(k=0;k < eModel -> numOfMotifs; k++)
    {
	if(eModel -> z[seq][k] > -1)
	{
	    eModel -> VPcounts[new][k]++;
	    eModel -> VPcounts[change][k]--;
	}
	else
	{
	    eModel -> VAcounts[new][k]++;
	    eModel -> VAcounts[change][k]--;
	}
    }
    return 1;
}



int sampleLeftEnd(FASTA *X, modelStruct *eModel, int mot, structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j;
    double sum = 0,t;
    double weights[3];
    double counts[4];
    double backL, backC;
    double R;
    int new;
    double probY, probN;
    double totalt;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return 0; 
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < allCons -> minsites) return 0; 
    backL=0;
    backC =0;
    memset(weights, 0, sizeof(double) * 3 );
    memset(counts, 0, sizeof(double) * 4);
    
    for(j = 0; j< 4; j ++)
    {
	t = eModel -> motifs[mot] -> counts[j][0];
	weights[1] = weights[1] + t * logWithCheck(t,1);
    }
    weights[1] = weights[1] -  eModel -> motifs[mot] -> total * logWithCheck(eModel -> motifs[mot] -> total,1);
    for(i =0 ; i< eModel -> numOfSequences; i++)
    {
	if( eModel->z[i][mot] < 0) continue;
	eModel->z[i][mot]--;
	if(pieceContainsNorOverlaps(X, eModel, i, mot))
	{
	    eModel->z[i][mot]++;
	    probY = eModel -> VPcounts[eModel -> I[i]][mot] /X->nonRepLength[i] ;
	    probN = eModel -> VAcounts[eModel -> I[i]][mot];
	    probY = probY / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    probN = probN / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    
	    weights[0] +=  logWithCheck(probN,1) - logWithCheck(probY,1);
	    
	}
	else
	{
	    eModel->z[i][mot]++;
	    counts[X -> seq[i][ eModel -> z[i][mot] - 1]] ++;
	    backL = backL + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] - 1][0],1);
	}
	backC = backC + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot]][0],1);
    }
    weights[1] = weights[1] - backC;
    weights[0] = weights[1] + weights[0];
    totalt = 0;
    for(j = 0; j< 4; j ++)
    {
	t = counts[j] + allCons -> pseudo;
	weights[0] = weights[0] + t * logWithCheck(t,1);
	totalt += t;
    }
    weights[0] = weights[0] - totalt * logWithCheck(totalt,1);
    weights[0] = weights[0] - backL;
    
    sum = normalizeWeights(weights,3);
    if(eModel -> motifs[mot] -> w == eModel -> maxW) weights[0] = 0;
    if(eModel -> motifs[mot] -> w <= eModel -> minW) weights[2] = 0;
    sum = 0;
    for(i = 0; i<3;i++)
    {
	sum = sum + weights[i];
    }
    if(sum == 0)
    {
	printf("sum of weights in motif width is 0");
	exit(0);
    }
    R = 0;
    while(R == 0)
    {
	R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
    }
    weights[0] = weights[0]/sum;
    for(i =0;i < 2 ;i++)
    {
	weights[i+1] = weights[i+1]/sum + weights[i];
    }
    new = binarySearchSample(R,weights,0,3 );
    if(new== 1) return 0;
    if(new == 0)
    {
	eModel -> motifs[mot] -> w++;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel->z[i][mot]--;
	    if(pieceContainsNorOverlaps(X, eModel, i, mot))
	    {
		eModel->z[i][mot] = -1;
		eModel ->totMotifs[i]--;
		eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
		eModel -> zDenom[i][mot] = 0;
		eModel -> VAcounts[eModel -> I[i]][mot] ++;
		eModel -> VPcounts[eModel -> I[i]][mot] --;
		continue;
	    }
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    if(new== 2)
    {
	if(eModel -> motifs[mot] -> w <= eModel -> minW)
	{
	    printf("widthr is small\n");
	}
	eModel -> motifs[mot] -> w--;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel->z[i][mot]++;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    return(1);
}




int trimLeftEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons)
{
    int i,j;
    double maxi,t;
    double weights[3];
    double counts[4];
    double backL, backC;
    int start, end;
    int new;
    double probY, probN;
    double totalt;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return 0; 
    backL=0;
    backC =0;
    memset(weights, 0, sizeof(double) * 3 );
    memset(counts, 0, sizeof(double) * 4);
    
    for(j = 0; j< 4; j ++)
    {
	t = eModel -> motifs[mot] -> counts[j][0];
	weights[1] = weights[1] + t * logWithCheck(t,1);
    }
    weights[1] = weights[1] -  eModel -> motifs[mot] -> total * logWithCheck(eModel -> motifs[mot] -> total,1);
    for(i =0 ; i< eModel -> numOfSequences; i++)
    {
	if( eModel->z[i][mot] < 0) continue;
	eModel->z[i][mot]--;
	if(pieceContainsNorOverlaps(X, eModel, i, mot))
	{
	    eModel->z[i][mot]++;
	    probY = eModel -> VPcounts[eModel -> I[i]][mot] /X->nonRepLength[i] ;

	    probN = eModel -> VAcounts[eModel -> I[i]][mot];
	    probY = probY / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    probN = probN / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    
	    weights[0] +=  logWithCheck(probN,1) - logWithCheck(probY,1);
	    
	}
	else
	{
	    eModel->z[i][mot]++;
	    counts[X -> seq[i][ eModel -> z[i][mot] - 1]] ++;
	    backL = backL + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] - 1][0],1);
	}
	backC = backC + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot]][0],1);
    }
    weights[1] = weights[1] - backC;
    weights[0] = weights[1] + weights[0];
    totalt = 0;
    for(j = 0; j< 4; j ++)
    {
	t = counts[j] + allCons -> pseudo;
	totalt += t;
	weights[0] = weights[0] + t * logWithCheck(t,1);
    }
    weights[0] = weights[0] - totalt * logWithCheck(totalt,1);
    weights[0] = weights[0] - backL;
    start = 0;
    end = 3;
    if(eModel -> motifs[mot] -> w == eModel -> maxW) start = 1;
    if(eModel -> motifs[mot] -> w <= eModel -> minW) end = 2;
    maxi = weights[start];
    new = start;
    for(i = start+1; i<end;i++)
    {
	if(maxi < weights[i])
	{
	    new = i;
	    maxi = weights[i];
	}
    }
    if(new== 1) return 0;
    if(new == 0)
    {
	eModel -> motifs[mot] -> w++;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel->z[i][mot]--;
	    if(pieceContainsNorOverlaps(X, eModel, i, mot))
	    {
		eModel->z[i][mot] = -1;
		eModel -> totMotifs[i]--;
		eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
		eModel -> zDenom[i][mot] = 0;
		eModel -> VAcounts[eModel -> I[i]][mot] ++;
		eModel -> VPcounts[eModel -> I[i]][mot] --;
		continue;
	    }
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    if(new== 2)
    {
	if(eModel -> motifs[mot] -> w <= eModel -> minW)
	{
	    printf("widthr is small\n");
	}
	eModel -> motifs[mot] -> w--;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel->z[i][mot]++;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    return(1);
}





int sampleRightEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j;
    double sum = 0,t;
    double weights[3];
    double counts[4];
    double backR, backC;
    double R;
    int new;
    double probY, probN;
    int currW;
    double totalt;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return 0;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < allCons -> minsites) return 0; 

    backR = 0;
    backC = 0;
    memset(weights, 0, sizeof(double) * 3 );
    memset(counts, 0, sizeof(double) * 4);
    currW = eModel -> motifs[mot] -> w;
    for(j = 0; j< 4; j ++)
    {
	t = eModel -> motifs[mot] -> counts[j][currW-1];
	weights[1] = weights[1] + t * logWithCheck(t,1);
    }
    weights[1] = weights[1] -  eModel -> motifs[mot] -> total * logWithCheck(eModel -> motifs[mot] -> total,1);
    eModel -> motifs[mot] -> w = currW + 1;
    for(i =0 ; i< eModel -> numOfSequences; i++)
    {
	if( eModel->z[i][mot] < 0) continue;
	if(pieceContainsNorOverlaps(X, eModel, i, mot))
	{
	    probY = eModel -> VPcounts[eModel -> I[i]][mot]/X->nonRepLength[i] ;
	    probN = eModel -> VAcounts[eModel -> I[i]][mot];
	    probY = probY / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    probN = probN / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    
	    weights[2] +=  logWithCheck(probN,1) - logWithCheck(probY,1);
	    
	}
	else
	{
	    counts[X -> seq[i][ eModel -> z[i][mot] + currW]] ++;
	    backR = backR + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] + currW][0],1);
	}
	backC = backC + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] + currW - 1][0],1);
    }
    eModel -> motifs[mot] -> w = currW;
    weights[1] = weights[1] - backC;
    weights[2] = weights[1] + weights[2];
    totalt = 0;
    for(j = 0; j< 4; j ++)
    {
	t = counts[j] + allCons -> pseudo;
	weights[2] = weights[2] + t * logWithCheck(t,1);
	totalt += t;
    }
    weights[2] = weights[2] - totalt * logWithCheck(totalt,1);
    weights[2] = weights[2] - backR;
    sum = normalizeWeights(weights,3);

    sum = 0;
    for(i = 0; i<3;i++)
    {
	sum = sum + weights[i];
    }
    if(sum == 0)
    {
	printf("sum of weights in motif width is 0");
	exit(0);
    }
    weights[0] = weights[0]/sum;
    for(i =0;i < 2 ;i++)
    {
	weights[i+1] = weights[i+1]/sum + weights[i];
    }
    R = 0;
    while(R == 0)
    {
	R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
    }
    if(eModel -> motifs[mot] -> w == eModel -> maxW) weights[2] = 0;
    if(eModel -> motifs[mot] -> w <= eModel -> minW) weights[0] = 0;

    new = binarySearchSample(R,weights,0,3 );
    if(new== 1) return 0;
    if(new == 2)
    {
	eModel -> motifs[mot] -> w = currW + 1;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    if(pieceContainsNorOverlaps(X, eModel, i, mot))
	    {
		eModel->z[i][mot] = -1;
		eModel -> totMotifs[i]--;

		eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
		eModel -> zDenom[i][mot] = 0;
		eModel -> VAcounts[eModel -> I[i]][mot] ++;
		eModel -> VPcounts[eModel -> I[i]][mot] --;
		continue;
	    }
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    if(new== 0)
    {
	if(eModel -> motifs[mot] -> w <= eModel -> minW)
	{
	    printf("widthr is small: this should not have happened\n");
	    exit(0);
	}
	eModel -> motifs[mot] -> w = currW - 1;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    return(1);
}



int trimRightEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons)
{
    int i,j;
    double maxi = 0,t;
    double weights[3];
    double counts[4];
    double backR, backC;
    int new, start,end;
    double probY, probN;
    int currW;
    double totalt;
    if(fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) < EPS) return 0; 
    backR = 0;
    backC = 0;
    memset(weights, 0, sizeof(double) * 3 );
    memset(counts, 0, sizeof(double) * 4);
    currW = eModel -> motifs[mot] -> w;
    for(j = 0; j< 4; j ++)
    {
	t = eModel -> motifs[mot] -> counts[j][currW-1];
	weights[1] = weights[1] + t * logWithCheck(t,1);
    }
    weights[1] = weights[1] -  eModel -> motifs[mot] -> total * logWithCheck(eModel -> motifs[mot] -> total,1);
    eModel -> motifs[mot] -> w = currW + 1;
    for(i =0 ; i< eModel -> numOfSequences; i++)
    {
	if( eModel->z[i][mot] < 0) continue;
	if(pieceContainsNorOverlaps(X, eModel, i, mot))
	{
	    probY = eModel -> VPcounts[eModel -> I[i]][mot]/X->nonRepLength[i] ;
	    probN = eModel -> VAcounts[eModel -> I[i]][mot];
	    probY = probY / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    probN = probN / (eModel -> VPcounts[eModel -> I[i]][mot] +eModel -> VAcounts[eModel -> I[i]][mot] );
	    
	    weights[2] +=  logWithCheck(probN,1) - logWithCheck(probY,1);
	    
	}
	else
	{
	    counts[X -> seq[i][ eModel -> z[i][mot] + currW]] ++;
	    backR = backR + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] + currW][0],1);
	}
	backC = backC + logWithCheck(X -> tableOfBackground[i][eModel -> z[i][mot] + currW - 1][0],1);
    }
    eModel -> motifs[mot] -> w = currW;
    weights[1] = weights[1] - backC;
    weights[2] = weights[1] + weights[2];
    totalt = 0;
    for(j = 0; j< 4; j ++)
    {
	t = counts[j];
	weights[2] = weights[2] + t * logWithCheck(t,1);
	totalt += t;
    }
    weights[2] = weights[2] - totalt * logWithCheck(totalt,1);
    weights[2] = weights[2] - backR;
    start = 0;
    end = 3;
    if(eModel -> motifs[mot] -> w == eModel -> maxW) end = 2;
    if(eModel -> motifs[mot] -> w <= eModel -> minW) start = 1;
    maxi = weights[start];
    new =start;
    for(i=start+1; i< end; i++)
    {
	if(weights[i] > maxi)
	{
	    maxi = weights[i];
	    new = i;
	}
    }
    if(new== 1) return 0;
    if(new == 2)
    {
	eModel -> motifs[mot] -> w = currW + 1;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    if(pieceContainsNorOverlaps(X, eModel, i, mot))
	    {
		eModel->z[i][mot] = -1;
		eModel -> totMotifs[i]--;
		eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
		eModel -> zDenom[i][mot] = 0;
		eModel -> VAcounts[eModel -> I[i]][mot] ++;
		eModel -> VPcounts[eModel -> I[i]][mot] --;
		continue;
	    }
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    if(new== 0)
    {
	if(eModel -> motifs[mot] -> w <= eModel -> minW)
	{
	    printf("widthr is small: this should not have happened\n");
	    exit(0);
	}
	eModel -> motifs[mot] -> w = currW - 1;
	for(i =0 ; i< eModel -> numOfSequences; i++)
	{
	    if( eModel->z[i][mot] < 0) continue;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][mot],1);
	    eModel -> zDenom[i][mot] =  scorePieceByTable(eModel,mot,i,eModel->z[i][mot],X);
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[i][mot],1);
	}
	computeCountsMotif(X, eModel, mot,allCons);
    }
    return(1);
}





int pickBestz(FASTA *X, modelStruct *eModel, int mot, int seq,structForAllConstants *allCons)
{
    int i;
    double weight ;
    double probY, probN ,curWeight, thisBack;
    int new;
    int change;
    change  = eModel ->z[seq][mot] ;
    new = -1;
    
    probY = eModel -> VPcounts[eModel -> I[seq]][mot]/X->nonRepLength[seq] ;
    probN = eModel -> VAcounts[eModel -> I[seq]][mot] ;
    probN = probN/probY;
    weight = probN;
    
    for(i=0;i < X->length[seq];i++)
    {
	eModel ->  z[seq][mot] = i;
  
	curWeight = scorePieceByMotif(X, eModel, mot,seq);
	if (curWeight == 0)
	{
	    eModel ->  z[seq][mot] = -1;
	    continue;
	}
	thisBack = scorePieceByTable(eModel,mot,seq,i,X);
	curWeight/= thisBack;
	
	eModel ->  z[seq][mot] = -1;
	if(curWeight > weight)
	{
	    weight = curWeight;
	    new = i;
	}
    }
    if(change == new)
    {
	eModel ->  z[seq][mot] = change;
	return 0;
    }
    if(new == -1)
    {
	eModel -> z[seq][mot] = change;
	removeMotifCounts(X,eModel, mot, seq,allCons);  
    	eModel -> z[seq][mot] = -1;
    	eModel -> totDenom -= logWithCheck(eModel -> zDenom[seq][mot],1);
    	if(eModel -> totDenom > 0)
    	    printf("debug 1028");
    	eModel -> zDenom[seq][mot] = 0;
	eModel -> z[seq][mot] = new;

    }
    else
    {
	if(change == -1)
	{
	    eModel -> z[seq][mot] = new;
	    addMotifCounts(X,eModel, mot, seq,allCons);
	    eModel -> zDenom[seq][mot] = scorePieceByTable(eModel,mot,seq,eModel -> z[seq][mot],X)  ;
	    
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[seq][mot],1) ;
	    if(eModel -> totDenom > 0)
		printf("debug 314");
	    if(new + eModel -> motifs[mot] -> w > X -> length[seq])
		printf("debug 316");

	}
	else
	{
	    eModel -> z[seq][mot] = change;
	    removeMotifCounts(X,eModel, mot, seq,allCons);  
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[seq][mot],1);
	    if(eModel -> totDenom > EPS)
		printf("debug 1054");
	    eModel -> zDenom[seq][mot] = 0;
	    eModel -> z[seq][mot] = new;
	    addMotifCounts(X,eModel, mot, seq, allCons);
	    eModel -> zDenom[seq][mot] = scorePieceByTable(eModel,mot,seq,eModel -> z[seq][mot],X)  ;
	    
	    eModel -> totDenom += logWithCheck(eModel -> zDenom[seq][mot],1) ;
	    if(eModel -> totDenom > 0)
		printf("debug 314");
	    if(new + eModel -> motifs[mot] -> w > X -> length[seq])
		printf("debug 316");
	}

    }
    
    if(new >= 0 && change == -1)
    {
	eModel -> VPcounts[eModel -> I[seq]][mot]++;
	eModel -> VAcounts[eModel -> I[seq]][mot]--;
	eModel -> totMotifs[seq]++;


    }

    if(new == -1)
    {
	eModel -> VAcounts[eModel -> I[seq]][mot]++;
	eModel -> VPcounts[eModel -> I[seq]][mot]--;
	eModel -> totMotifs[seq]--;
    }
    return 1;
}


int binarySearchSample(double value, double *weights, int start, int end)
{
    int mid;
    if(start >= end)
    {
	return(start);
    }
    else
    {
	mid = (start+end)/2;
	if(weights[mid] < value) return(binarySearchSample(value, weights, mid+1,end));
	else
	{
	    if(mid == 0) return(0);
	    else
	    {
		if(value > weights[mid-1]) return(mid);
		else return(binarySearchSample(value, weights, start,mid));
	    }
	}
    }
}

double fullPosterior(FASTA *X, modelStruct *eModel,structForAllConstants *allCons )
{
    int  j,k;
    double partI, partZ, partBack;
    double sum = 0;
    partBack = -eModel -> totDenom;
    partI = 0;
    partZ = 0;
    for(k = 0; k < eModel -> numOfMotifs; k++)
    {
	partZ = partZ + eModel -> partZ[k];
    }
    for(j =0; j < eModel -> numOfRegProgs; j++)
    {
	partI = partI + (eModel -> gammaCounts[j] - allCons -> pseudoGamma) * logWithCheck( eModel ->gammaCounts[j],1);
	for(k =0 ;k<eModel -> numOfMotifs; k++)
	{
	    partI = partI + (eModel -> VAcounts[j][k] - allCons -> pseudoA) * logWithCheck(eModel -> VAcounts[j][k],1);
	    partI = partI + (eModel -> VPcounts[j][k] - allCons -> pseudoP)* logWithCheck(eModel -> VPcounts[j][k],1);
	    sum = eModel -> VAcounts[j][k]  + eModel -> VPcounts[j][k] ;
	    partI = partI - (sum - allCons -> pseudoA - allCons -> pseudoP) * logWithCheck(sum,1);
	}	
    }
    partI = partI - eModel -> numOfSequences * logWithCheck(eModel -> numOfSequences , 1);
    return(partBack + partI + partZ );
}


int clubNonMotifs(FASTA *X, modelStruct *eModel,structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j,k,current,new;
    double *empty;
    double sum, R;
    sum = 0;
    empty = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset(empty, 0, sizeof(double ) * eModel -> numOfRegProgs);
    for(i =0;i< eModel->numOfRegProgs; i++)
    {
	if(fabs(eModel ->gammaCounts[i] - allCons -> pseudoGamma) < EPS)
	{
	    empty[i] = 1;
	    sum ++;
	}
    }
    if(sum == 0 )
    {
	free(empty);
	return 0;
    }
    empty[0] = empty[0] / sum;
    for(j =0; j < eModel -> numOfRegProgs -1 ;j++)
    {
	empty[j+1] = empty[j+1]/sum + empty[j];
    }
    for(i=0;i< eModel -> numOfSequences; i++)
    {
	if(eModel -> totMotifs[i] == 0)
	{
	    current = eModel -> I[i];
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
		eModel -> VAcounts[current][k]--;
		eModel -> VAcounts[new][k]++;
	    }
	}
    }
    free(empty);
    return(1);
}



int clubNonMotifsSeq(FASTA *X, modelStruct *eModel,int seq,structForAllConstants *allCons,unsigned int *seeds)
{
    int i,j,k,current,new;
    double *empty;
    double sum, R;
    sum = 0;
    empty = ALLOC(sizeof(double) * eModel -> numOfRegProgs);
    memset(empty, 0, sizeof(double ) * eModel -> numOfRegProgs);
    for(i =0;i< eModel->numOfRegProgs; i++)
    {
	if(fabs(eModel ->gammaCounts[i] - allCons -> pseudoGamma) < EPS)
	{
	    empty[i] = 1;
	    sum ++;
	}
    }
    if(sum == 0 )
    {
	free(empty);
	return 0;
    }
    empty[0] = empty[0] / sum;
    for(j =0; j < eModel -> numOfRegProgs -1 ;j++)
    {
	empty[j+1] = empty[j+1]/sum + empty[j];
    }
    i = seq;
    current = eModel -> I[i];
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
	eModel -> VAcounts[current][k]--;
	eModel -> VAcounts[new][k]++;	
    }
    free(empty);
    return(1);
}

int clubNonMotifsTogether(FASTA *X, modelStruct *eModel,structForAllConstants *allCons)
{
    int i,k,current;
    int empty;
    double sum;
    int toReturn = 0;
    sum = 0;
    for(i =0;i< eModel->numOfRegProgs; i++)
    {
	if(fabs(eModel ->gammaCounts[i] - allCons -> pseudoGamma) < EPS)
	{
	    empty= i;
	    sum ++;
	    break;
	}
    }
    if(sum == 0 )
    {
	return toReturn;
    }
    for(i=0;i< eModel -> numOfSequences; i++)
    {
	if(eModel -> totMotifs[i] == 0)
	{
	    current = eModel -> I[i];
	    eModel -> I[i] = empty;
	    eModel ->gammaCounts[current]--;
	    eModel ->gammaCounts[empty]++;
	    for(k=0;k<eModel->numOfMotifs;k++)
	    {
		eModel -> VAcounts[current][k]--;
		eModel -> VAcounts[empty][k]++;
	    }
	    toReturn = 1;
	}
    }
    return(toReturn);
}
