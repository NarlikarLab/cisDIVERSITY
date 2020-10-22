#include "model.h"
#include <getopt.h>

modelStruct *learnModel(FASTA *X, modelStruct *eModel, FILE *fp1, FILE *fp2,structForAllConstants *allCons,unsigned int *seeds)
{
    int i,k;
    modelStruct *best, *best1;
    int iter;
    double current;
    double bestSoFar,bestSoFar1;
    double tempForPrinting;
    int change;
    int noimprovement;
    double temp = 1;
    int newi;
    int initW;
    int inits ; //20;

    int *givenMotifs;
    givenMotifs = ALLOC(sizeof (int) * eModel -> numOfMotifs);
    for(k = 0; k < eModel -> numOfMotifs;k++)
    {
	givenMotifs[k] = k;
    }
    initW = eModel -> motifs[0] -> w;
    best = makeACopyOfModel(eModel);
    best1 = makeACopyOfModel(eModel);
    if(allCons -> verbose >= 1)
    {
	printv(eModel, fp2,allCons);
	printAllMotifs(best, fp2,allCons);
	printAllConsensus(best, fp2,allCons);
    }
    bestSoFar =  fullPosterior(X, best,allCons);
    tempForPrinting = bestSoFar;
    current = bestSoFar;
    if(fp1 != NULL)    fprintf(fp1, "-1 %lf\n", bestSoFar);
    change = 0;
    iter = 0;
    noimprovement = 0;

    while(iter < allCons -> iter && noimprovement <= allCons -> improve)
    {
	if(fp1 != NULL) fprintf(fp1, "%d %lf\n", iter, current);

	if(iter == allCons -> improve) temp = 1;
	if(temp < 1) temp = 1;
	if(iter % 10 == 0 )
	{
	    for(k =0; k< eModel -> numOfMotifs; k++)
	    {
		
		change = sampleLeftEnd(X, eModel, k,allCons,seeds);
		if(change == 1)
		{
		    current = fullPosterior(X, eModel,allCons);
		    if(current > bestSoFar)
		    {
			justCopyModel(eModel,best);
			bestSoFar = current;
			noimprovement= 0;
		    }
		}
		resetEmptyMotifWidth(eModel, k, allCons);
		
		change = sampleRightEnd(X, eModel, k,allCons,seeds);
		
		if(change == 1)
		{
		    current = fullPosterior(X, eModel,allCons);
		    if(current > bestSoFar)
		    {
			justCopyModel(eModel,best);
			bestSoFar = current;
			noimprovement= 0;
			
		    }
		}
		resetEmptyMotifWidth(eModel, k, allCons);
	    }
	    clubNonMotifs(X,eModel,allCons,seeds);
	}
	
	for(i =0;i< eModel -> numOfSequences;i++)
	{
	    for(k = 0;k < eModel -> numOfMotifs; k++)
	    {
		change = samplezAO(X, eModel, givenMotifs[k],i,temp,allCons,seeds);	
		resetEmptyMotifWidth(eModel, givenMotifs[k], allCons);

		if(change == 1)  
		{
		    if(eModel -> totMotifs[i] == 0)
			clubNonMotifsSeq(X,eModel,i,allCons,seeds); 
		    current = fullPosterior(X, eModel,allCons);
		    if(current > bestSoFar)
		    {
			justCopyModel(eModel,best);

			bestSoFar = current;
			noimprovement= 0;
		    }
		}
	    }
	
	

	    change = sampleI(X, eModel, i, 1,seeds);
		
	    if(change == 1)
	    {
		    
		current = fullPosterior(X, eModel,allCons);
		if(current > bestSoFar)
		{
		    justCopyModel(eModel,best);
		    
		    bestSoFar = current;
		    noimprovement= 0;
		}
	    }
	    
	}	    
	iter ++;
	if(tempForPrinting < bestSoFar)
	{
	    if(allCons -> verbose >= 1)
	    {
		fprintf(fp2,"\nIter %d: Current posterior %lf, best %lf\n",iter, current,bestSoFar);
		fprintf(fp2,"\n");
		printv(best, fp2,allCons);
		printAllConsensus(best, fp2,allCons);
		fflush(fp2);
	    }
	    tempForPrinting = bestSoFar;
	    noimprovement = 0;
	}
	else
	{
	    noimprovement++;
	}
	if(allCons -> verbose >= 1)
	{
	    fprintf(fp2, ".");

	}

    }
    cleanUp(best, initW,allCons);
    clubNonMotifsTogether(X,best,allCons) ;
    removeSmallProgsCompletely(best, allCons -> minseqs,allCons,seeds);
    if(allCons -> entropy == 1)
	removeMotifsWithLowEntropy(best,allCons,fp2);
    
    iter = 0;
    justCopyModel(best, eModel);
    current = fullPosterior(X, eModel,allCons);
    bestSoFar = current;
    while(iter < allCons -> improve)
    {
	for(i =0;i< eModel -> numOfSequences;i++)
	{
	    for(k = 0;k < eModel -> numOfMotifs; k++)
	    {
		change = samplezAO(X, eModel, givenMotifs[k],i,temp,allCons,seeds);	
		resetEmptyMotifWidth(eModel, givenMotifs[k], allCons);
		
		if(change == 1)  
		{
		    if(eModel -> totMotifs[i] == 0)
			clubNonMotifsSeq(X,eModel,i,allCons,seeds); 
		    current = fullPosterior(X, eModel,allCons);
		    if(current > bestSoFar)
		    {
			justCopyModel(eModel,best);
			bestSoFar = current;			
		    }
		}
	    }
	}
	iter ++;
    }

    iter = 0;
    justCopyModel(best, eModel);
    current = fullPosterior(X, eModel,allCons);
    bestSoFar = current;
    initializeRandomI(eModel,allCons,seeds);
    computeVAVP(eModel,allCons);
    current = fullPosterior(X, eModel,allCons);
    fflush(fp2);
    iter = 0;
    while(iter < allCons -> iter * 2)
    {
	if(fp1 != NULL) fprintf(fp1, "%d %lf\n", iter, current);
	for(i =0;i< eModel -> numOfSequences;i++)
	{
	    change = sampleI(X, eModel, i, 1,seeds);
	    if(change == 1)
	    {
		current = fullPosterior(X, eModel,allCons);
		if(current > bestSoFar)
		{
		    justCopyModel(eModel,best);
		    bestSoFar = current;
		}
	    }
	    
	}
	iter++;

    }
    if(allCons -> entropy == 1)
	removeMotifsWithLowEntropy(best,allCons,fp2);

    bestSoFar =  fullPosterior(X, best,allCons);
    justCopyModel(best, eModel);

    noimprovement = 1;
    if(allCons -> verbose >=1) fprintf(fp2, "\nDoing hill climbing\n");
    while(noimprovement == 1)
    {
	noimprovement = 0;
	for(i =0;i< best -> numOfSequences;i++)
	{
	    for(k = 0;k < best -> numOfMotifs; k++)
	    {
		change = pickBestz(X, eModel, k,i,allCons);
		if(change == 1)
		{
		    clubNonMotifsTogether(X,eModel,allCons);
		    current = fullPosterior(X, eModel,allCons);
		    if(current <= bestSoFar)
		    {
			noimprovement= 0;
			justCopyModel( best, eModel);
			
			continue;
		    }
		    bestSoFar = current;
		    justCopyModel( eModel,best);
		    if(allCons -> verbose >=1) fprintf(fp2, ".");
		    noimprovement = 1;
		}
		
	    }
	    change =  pickBestI(X, eModel,  i,allCons);
	    if(change == 1)
	    {
		clubNonMotifsTogether(X,eModel,allCons);
		current = fullPosterior(X, eModel,allCons);
		if(current <= bestSoFar)
		{
		    noimprovement= 0;
		    justCopyModel( best, eModel);

		    continue;		    
		}
		bestSoFar = current;
		justCopyModel( eModel,best);
		if(allCons -> verbose >=1) fprintf(fp2, ".");
		noimprovement = 1;
	    }
	}
	
	for(k=0;k < best -> numOfMotifs; k++)
	{
	    change = trimLeftEnd(X, eModel, k,allCons);
	    if(change == 1)
	    {
		clubNonMotifsTogether(X,eModel,allCons);

		current = fullPosterior(X, eModel,allCons);
		if(current <= bestSoFar)
		{
		    noimprovement= 0;
		    justCopyModel( best, eModel);

		    continue;
		    
		}
		bestSoFar = current;
		justCopyModel( eModel,best);
		if(allCons -> verbose >=1) fprintf(fp2, ".");
		noimprovement = 1;

	    }
	    
	    change = trimRightEnd(X, eModel, k,allCons);
	    if(change == 1)
	    {
		clubNonMotifsTogether(X,eModel,allCons);

		current = fullPosterior(X, eModel,allCons);
		if(current <= bestSoFar)
		{
		    noimprovement= 0;
		    justCopyModel( best, eModel);
		    continue;
		}
		bestSoFar = current;
		justCopyModel( eModel,best);
		if(allCons -> verbose >=1) fprintf(fp2, ".");
		noimprovement = 1;
	    }
	    
	}
	if(fp1 != NULL) fprintf(fp1, "0 %lf\n", current);
    }
    if(allCons -> verbose >=1) fprintf(fp2,"\nBest before cleaning: %lf\n", fullPosterior(X, best,allCons));

    
    cleanUp(best, initW,allCons);
    clubNonMotifsTogether(X,best,allCons) ;
    removeSmallProgsCompletely(best, allCons -> minseqs,allCons,seeds);
    if(allCons -> entropy == 1)
	removeMotifsWithLowEntropy(best,allCons,fp2);
    
    justCopyModel(best,eModel);
    current = fullPosterior(X, eModel,allCons);
    bestSoFar = current;

    inits = eModel -> numOfRegProgs;
    for(newi = 1; newi < inits + 1; newi++)
    {
	justCopyModel(best,eModel);
	justCopyModel(best,best1);
	current = fullPosterior(X, eModel,allCons);
	bestSoFar1 = current;

	initializeRandomIspecial(eModel,allCons,seeds,newi);
	computeVAVP(eModel,allCons);
	current = fullPosterior(X, eModel,allCons);

	if(allCons -> verbose >=1)
	{
	    fprintf(fp2,"posterior after randomizing for %d programs: %lf\n",newi ,current);
	    fflush(fp2);
	}
	if(newi == 1) continue;
	iter = 0;
	while(iter < allCons -> iter * 2)
	{
	    if(fp1 != NULL) fprintf(fp1, "%d %lf\n", iter, current);
	    for(i =0;i< eModel -> numOfSequences;i++)
	    {
		change = sampleIlimitedProgs(X, eModel, i, newi,seeds);
		if(change == 1)
		{
		    current = fullPosterior(X, eModel,allCons);
		    if(current > bestSoFar1)
		    {
			justCopyModel(eModel,best1);
			bestSoFar1 = current;

		    }
		}
		
	    }
	    iter++;

	}
	removeSmallProgsCompletely(best1, allCons -> minseqs,allCons,seeds);
	clubNonMotifsTogether(X,best1,allCons);
	current = fullPosterior(X, best1,allCons);
	if(current > bestSoFar)
	{
	    justCopyModel(best1,best);
	    bestSoFar = current;
		    
	}
	
    }



    free(givenMotifs);
    return(best);
}
 
