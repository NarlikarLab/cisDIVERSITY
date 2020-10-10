
#include "model.h"

int samplezAO(FASTA *X, modelStruct *eModel, int mot, int seq,double temp,structForAllConstants *allCons,unsigned int *seeds)
{
    int i;
    double sum = 0;
    double *weights ;
    double probY, probN ,R, thisBack;
    int new;
    int change;
    int n, *overlaps;
    int who = RUSAGE_SELF; 
    struct rusage usage; 
    int ret; 


    if (X->nonRepLength[seq] == 0) X->nonRepLength[seq] = 1;
    overlaps = ALLOC(sizeof(int)*eModel ->numOfMotifs);
    change  = eModel ->z[seq][mot] ;
    weights = ALLOC(sizeof(double) * (X->length[seq] + 1));

    
    if(eModel -> z[seq][mot] > -1)
    {
	removeMotifCounts(X,eModel, mot, seq,allCons); 
	eModel -> z[seq][mot] = -1;
	eModel -> totDenom -= logWithCheck(eModel -> zDenom[seq][mot],1);
	if(eModel -> totDenom > EPS)
	    printf("debug 869");
	eModel -> zDenom[seq][mot] = 0;
	probY = (eModel -> VPcounts[eModel -> I[seq]][mot]-1)/X->nonRepLength[seq] ;
	probN = eModel -> VAcounts[eModel -> I[seq]][mot]/temp;
	eModel -> totMotifs[seq]--;
    }
    else
    {
	probY = eModel -> VPcounts[eModel -> I[seq]][mot]/X->nonRepLength[seq] ;
	probN = (eModel -> VAcounts[eModel -> I[seq]][mot]-1)/temp;
    }
    probN = probN * pow(eModel -> motifs[mot] -> total,eModel -> motifs[mot] ->w) ;
    probN = probN/probY;

    for(i=0;i < X->length[seq];i++)
    {
	eModel ->  z[seq][mot] = i;
  
	weights[i] = scorePieceByMotifAO(X, eModel, mot,seq);
  
	if (weights[i] == 0)
	{
	    eModel ->  z[seq][mot] = -1;
	    continue;
	}
	thisBack = scorePieceByTable(eModel,mot,seq,i,X);
	if(thisBack == 0)
	{
	    fprintf(stderr,"Something went wrong in computing the background score\n");
	    printf("seq %d, pos %d,  back %f", seq,i, X -> tableOfBackground[seq][i][eModel -> motifs[mot] ->w -1 ]);
	    exit(0);
	}
	weights[i] /= thisBack;
	sum += weights[i];
	eModel ->  z[seq][mot] = -1;
    }
   
    sum += probN;
    if(sum == 0)
    {
	printf("DEBUG.... sum is 0:sampling function.c; probN is %lf", probN);
	printf("eModel -> VAcounts[eModel -> I[seq]][mot] is %lf", eModel -> VAcounts[eModel -> I[seq]][mot]);
	exit(0);
	    
    }
    weights[X->length[seq]] = probN;
    weights[0] = weights[0]/sum;
    for(i =0;i < X->length[seq] ;i++)
    {
	weights[i+1] = weights[i+1]/sum + weights[i];
    }
    R = 0;
    while(R == 0)
    {
	R =  (0.0 + rand_r(seeds))/(float)RAND_MAX;
    }
    new = binarySearchSample(R,weights,0,X->length[seq] + 1 );
    if(new < X->length[seq])
    {
	eModel -> z[seq][mot] = new;
	addMotifCounts(X,eModel, mot, seq,allCons);
	eModel -> totMotifs[seq]++;
 
	n = pieceOverlapsWith(X,eModel,seq,mot,overlaps);
	while(n >0)
	{
	    n--;
	    removeMotifCounts(X,eModel, overlaps[n], seq,allCons);
	    eModel -> z[seq][overlaps[n]] = -1;
	    eModel -> totDenom -= logWithCheck(eModel -> zDenom[seq][ overlaps[n]],1) ;
	    eModel -> zDenom[seq][ overlaps[n]] = 0;
	    eModel -> VAcounts[eModel -> I[seq]][overlaps[n]]++;
	    eModel -> VPcounts[eModel -> I[seq]][overlaps[n]]--;
	    eModel -> totMotifs[seq]--;
	}
	eModel -> zDenom[seq][mot] = scorePieceByTable(eModel,mot,seq,eModel -> z[seq][mot],X)  ;

	eModel -> totDenom += logWithCheck(eModel -> zDenom[seq][mot],1) ;
	if(eModel -> totDenom > 0)
	    printf("debug 314");
	if(new + eModel -> motifs[mot] -> w > X -> length[seq])
	    printf("debug 316");


    }
    else
    {
	new = -1;
	eModel -> z[seq][mot] = new;
    }
    free(weights);
    free(overlaps);

    if(new >= 0 && change == -1)
    {
	eModel -> VPcounts[eModel -> I[seq]][mot]++;
	eModel -> VAcounts[eModel -> I[seq]][mot]--;
    }
    if(new == -1 && change > -1)
    {
	eModel -> VAcounts[eModel -> I[seq]][mot]++;
	eModel -> VPcounts[eModel -> I[seq]][mot]--;

    }
    
    if(change == new) return 0;

    return 1;

    
    
}

