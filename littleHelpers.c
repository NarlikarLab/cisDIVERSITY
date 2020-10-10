
#include "model.h"
#include <getopt.h>

void resetEmptyMotifWidth(modelStruct *eModel, int mot, structForAllConstants *allCons)
{

    if( fabs(eModel -> motifs[mot] -> total - allCons -> pseudo * 4) <= EPS)
    {
	if(eModel -> motifs[mot] -> w >  allCons -> initw)
	    eModel -> motifs[mot] -> w = allCons -> initw;
	computePartZ(eModel,mot,allCons);
    }
}


int numberOfCols(modelStruct *eModel,int mot, double cutoff)
{
    int i, j, tot = 0;
    double ans = 2.0;
    for(i = 0; i< eModel ->motifs[mot]->w; i++)
    {
        ans = 2.0;
        for(j = 0; j< 4; j ++)
        {
            ans = ans + eModel -> motifs[mot] -> counts[j][i]/eModel -> motifs[mot] -> total * logWithCheck(eModel -> motifs[mot] -> counts[j][i]/eModel -> motifs[mot] -> total,2 );
        }
        if(ans > cutoff) tot++;
    }
    return tot;
}

void removeMotifsWithLowEntropy(modelStruct *eModel, structForAllConstants *allCons,FILE *fp2)
{
    int i,j,k;
    int cols;
    for(k =0 ; k< eModel -> numOfMotifs; k++)
    {
        if( fabs(eModel -> motifs[k] -> total - allCons -> pseudo * 4) <= EPS &&  fp2 != NULL)
        {
            fprintf(fp2,"motif %d is empty\n", k);

            continue;
        }

        cols = numberOfCols(eModel,k,1.0);
	if(fp2 != NULL)
	    fprintf(fp2,"Motif %d has %d columns\n", k,cols);
        if(cols <= 3)
        {
	    if(fp2 != NULL)
	    {
		fprintf(fp2,"Removing motif %d\n", k);
		printConsensus(eModel -> motifs[k], fp2);
		fprintf(fp2,"\n");
		fflush(fp2);
	    }
            for(i = 0; i< eModel -> numOfSequences;i++)
            {
                if (eModel -> z[i][k] == -1) continue;
                eModel -> z[i][k] = -1;
                eModel -> totDenom -= logWithCheck(eModel -> zDenom[i][k],1);
                eModel -> zDenom[i][k] = 0;
                eModel -> totMotifs[i]--;
            }
            eModel -> motifs[k] -> total = allCons -> pseudo * 4;

            for(i =0; i< eModel -> numOfRegProgs; i++)
            {
                eModel -> VAcounts[i][k] = allCons -> pseudoA + (eModel -> gammaCounts[i]-  allCons -> pseudoGamma  );
                eModel -> VPcounts[i][k] = allCons -> pseudoP;
            }
            eModel -> motifs[k] -> w = allCons -> initw;
            eModel -> partZ[k] = 0;
            for(i = 0; i< eModel ->motifs[k]->w; i++)
            {
                for(j = 0; j< 4; j ++)
                {
                    eModel -> motifs[k] -> counts[j][i] = allCons -> pseudo;
                }
            }

        }
    }

}

void selectionSort(int **array, int n)
{
    int i,j,cnt;
    int *temp;
    int max, maxid;
    temp = ALLOC(sizeof(int) * n);
    for(i=0;i<n;i++)
    {
	temp[i] = array[i][0];
    }
    i = 0;
    cnt =0;
    while(i<n)
    {
	if(temp[i] == -1)
	{
	    i = i+1;
	    continue;
	}
	max = temp[i];
	maxid = i;
	for(j=0;j<n;j++)
	{
	    if(temp[j] == -1) continue;
	    if(temp[j] > max)
	    {
		max = temp[j];
		maxid = j;
	    }
	}
	array[cnt][0] = max;
	array[cnt][1] = maxid;
	array[maxid][2] = cnt;
	cnt ++;
	temp[maxid] = -1;
    }
    free(temp);
}

void *ALLOC(size_t size) 
{
    void *foo;
    if(!(foo = malloc(size))) 
    {
	fprintf(stderr, "malloc failed.\n");
	exit(0);
    }
    return foo;
}

void fileError(char *filename, char *errorString)
{
    fprintf(stderr, "%s %s.\n", filename, errorString);
    exit(0);
}
void freeBacks(double **b, FASTA *X)
{
    int j;
    for(j=0;j<X -> num; j++)
    {
	free(b[j] );
    }
    free(b);
}
void freecountStruct(countStruct *motif)
{
    int i;
    for(i = 0; i< 4; i++)
    {
	free(motif -> counts[i]);
    }
    free(motif->counts);
    free(motif);
}

void freePhi(phiStruct *motif)
{
    int i;
    for(i = 0; i< 4; i++)
    {
	free(motif -> phi[i]);
    }
    free(motif -> phi);
    free(motif);
}

void freepssms(pssms *motifs)
{
    int i;
    for(i = 0; i< motifs -> number; i++)
    {
	freePhi(motifs -> motif[i]);
    }
    free(motifs -> motif);
    free(motifs);
}

void freeMM(MMStruct *M)
{
    int i;
    for(i=0; i<= M->order; i++)
    {
	free(M->probs[i]);
    }
    free(M->probs);
    free(M);
}

void freeFASTA(FASTA *X)
{
    int i,j;
    for(i =0;i< X -> num; i++)
    {
	free(X -> seq[i]);
	free(X -> header[i]);
	free(X -> noN[i]);
	for(j =0;j< X -> length[i]; j++)
	    free(X -> tableOfBackground[i][j]);
	free(X ->tableOfBackground[i]);
    }
    free(X -> length);
    free(X -> seq);
    free(X -> header);
    free(X -> noN);
    free(X ->tableOfBackground);
    free(X);
}


void freeWords(wordsStruct *N)
{
    int i;
    for(i=0; i<= N->order; i++)
    {
	free(N->wordCounts[i]);
    }
    free(N->wordCounts);
    free(N);
}

void freeModel(modelStruct *eModel)
{
    int i,k;
    for(k = 0;k < eModel -> numOfMotifs; k++)
    {
	freecountStruct(eModel -> motifs[k]);
    }
    free(eModel -> motifs);
    free( eModel -> totMotifs );
    for(i = 0;i < eModel -> numOfSequences; i++)
    {
	free(eModel -> z[i]);
	free(eModel -> zDenom[i]);
    }
    free(eModel -> z);
    free(eModel -> zDenom);
    free(eModel -> I);
    free(eModel -> gammaCounts);
    free(eModel -> partZ);
    for(i = 0; i <  eModel -> numOfRegProgs; i++)
    {
	free(eModel -> VAcounts[i]);
	free(eModel -> VPcounts[i]);
    }
    free(eModel -> VAcounts);
    free(eModel -> VPcounts);
    free(eModel);
}


long int getIndexInt(int *seq, int start, int length)
{
    long int index = 0;
    int alpha, k;
    for(k=0;k<length+1; k++)
    {
	switch (seq[start+k])
	{
	case 0:
	    alpha = 0;
	    break;
	case 1:
	    alpha = 1;
	    break;
	case 2:
	    alpha = 2;
	    break;
	case 3:
	    alpha = 3;
	    break;
	default:
//			printf("\n%c %d %d\n", seq[start+k], start,k);
	    //exit(1);
	    return(-1);
	}
	index = index + pow(4, length-k) * alpha;
    }
    return(index);
    
}


double logWithCheck(double x, float base)
{
    if(x > 0) 
    {
	if(base == 1) return log(x);
	if(base < 0) return (-1000);
	return(log(x)/log(base));
    }
    else 
    {
	return(0);
    }
} 



double min(double a, double b)
{
    if(a < b )
	return a;
    return b;
}
double max(double a, double b)
{
    if(a < b )
	return b;
    return a;
}

void setMaxSequenceWithoutN(FASTA *X, int seq)
{
    int i,j;
    int prev = -1;
    for(i=0;i< X ->length[seq];i++)
    {
	if(X -> seq[seq][i] == 4)
	{
	    for(j = prev + 1;j < i;j++)
	    {
		X -> noN[seq][j] = i - j; 
	    }
	    prev = i;
	}
    }
    for(j = prev + 1;j < i;j++)
    {
	X -> noN[seq][j] = i - j; 
    }

}


int maxPieceWithoutN(int *seq, int pos, int order )
{
    int i,start;
    start = max(0, pos-order);
    for(i = pos; i >= start; i--)
    {
	if(seq[i] == 4)
	    return(i+1);
    }
    return(start);
}


void normalizeMM(MMStruct *model)
{
    int i,j,k;
    double sum;
    for(i =0 ; i< model->order+1; i++)
    {
	for(j=0; j< pow(4,i); j++)
	{
	    sum = 0.0;
	    for( k =0;k<4;k++)
	    {
//		model->probs[i][j*4+k] = model->probs[i][j*4+k] + 0.1;
//		model->probs[i][j*4+k] = pow(model->probs[i][j*4+k],0.9);
		sum += model->probs[i][j*4+k];
	    }
	    if(sum == 0)
	    {
		sum = 1;
	    }
	    for( k =0;k<4;k++)
	    {
		model->probs[i][j*4+k] /= sum;
//		model->probs[i][j*4+k] = pow(model->probs[i][j*4+k],0.99);
	    }
	}
    }
} 



void normalizeMMwithPseudoCount(MMStruct *model)
{
    int i,j,k;
    double sum;
    for(i =0 ; i< model->order+1; i++)
    {
	for(j=0; j< pow(4,i); j++)
	{
	    sum = 0.0;
	    for( k =0;k<4;k++)
	    {
//		model->probs[i][j*4+k] = model->probs[i][j*4+k] + 0.1;
//		model->probs[i][j*4+k] = pow(model->probs[i][j*4+k],0.9);
		sum += model->probs[i][j*4+k];
	    }
	    
	    for( k =0;k<4;k++)
	    {
		model->probs[i][j*4+k] = (model->probs[i][j*4+k] + 1) /(sum + 4);
//		model->probs[i][j*4+k] = pow(model->probs[i][j*4+k],0.99);
	    }
	}
    }
} 


double normalizeWeights(double *weights, int n)
{
    int i;
    double ans,diff;
    double min = DBL_MAX;
    double max = -DBL_MAX;
	
    for(i=0;i<n;i++)
    {
	if(weights[i] < min)
	    min = weights[i];
	if(weights[i] > max)
	    max = weights[i];
    }
    diff = (max + min)/2;
    ans = 0;
    for(i=0;i<n;i++)
    {
	weights[i] = exp(weights[i] - diff);
	ans = ans + weights[i];
    }
    return(ans);
}

int pieceContainsN(int *seq, int start, int width)
{
    int i;
    for(i=start; i < start+width; i++)
    {
	if(seq[i] == 4)
	    return(1);
    } 
    return(0);
}

int pieceContainsNquick(FASTA *X, int start, int width, int seq)
{
    if(X -> noN[seq][start] < width) return(0);
    return(1);
}


int pieceOverlapsWith(FASTA *X,modelStruct *eModel,int seq,int mot,  int *overlaps)
{
    int k;
    int startMot, endMot,startMot1,endMot1;
    int *done;
    int n;
    startMot = eModel -> z[seq][mot];
    endMot = eModel -> z[seq][mot] +  eModel -> motifs[mot] -> w -1;
    n = 0;
    done = ALLOC(sizeof(int)*eModel ->numOfMotifs);
    memset(done,0, sizeof(int ) * eModel -> numOfMotifs);
    memset(overlaps,0, sizeof(int ) * eModel -> numOfMotifs);
    
    for(k = 0; k< eModel -> numOfMotifs;k++)
    {
	if( k == mot) continue;
	if( eModel -> z[seq][k] == -1) continue;
	startMot1 = eModel -> z[seq][k];
	endMot1 = eModel -> z[seq][k] +  eModel -> motifs[k] -> w - 1;
	if(startMot <= endMot1 && startMot >= startMot1)
	{
	    overlaps[n] = k;
	    done[k] =1;
	    n++;
	    continue;
	}
	if(endMot <= endMot1 && endMot >= startMot1)
	{
	    overlaps[n] = k;
	    done[k] =1;
	    n++;
	    continue;
	}
	if(startMot1 <= endMot && startMot1 >= startMot)
	{
	    overlaps[n] = k;
	    done[k] =1;
	    n++;
	    continue;
	}
	if(endMot1 <= endMot && endMot1 >= startMot) 
	{
	    overlaps[n] = k;
	    done[k] =1;
	    n++;
	    continue;
	}

    }
    if(X -> revFlag == 1)
    {
	for(k = 0; k< eModel -> numOfMotifs;k++)
	{
	    if( k == mot) continue;
	    if(done[k] == 1) continue;
	    if( eModel -> z[seq][k] == -1) continue;
	    endMot1 = reverseLocation(X, seq,eModel -> z[seq][k]);
	    startMot1 = endMot1 - eModel -> motifs[k] -> w;
	    endMot1 --;
	    
	    if(startMot <= endMot1 && startMot >= startMot1)
	    {
		overlaps[n] = k;
		n++;
		continue;
	    }
	    if(endMot <= endMot1 && endMot >= startMot1) 
	    {
		overlaps[n] = k;
		n++;
		continue;
	    }
	    if(startMot1 <= endMot && startMot1 >= startMot)
	    {
		overlaps[n] = k;
		n++;
		continue;
	    }
	    if(endMot1 <= endMot && endMot1 >= startMot) 
	    {
		overlaps[n] = k;
		n++;
		continue;
	    }

	}
    }
    free(done);
    return(n);
}

int pieceContainsNorOverlaps(FASTA *X,modelStruct *eModel,int seq,int mot)
{
    int k;
    int startMot, endMot,startMot1,endMot1;
    
    startMot = eModel -> z[seq][mot];
    if(startMot < 0) return(1);
    endMot = eModel -> z[seq][mot] +  eModel -> motifs[mot] -> w - 1;
    if(endMot > X -> length[seq]) return(1);
    if(X -> noN[seq][startMot] < eModel -> motifs[mot] -> w) return(1);
    for(k = 0; k< eModel -> numOfMotifs;k++)
    {
	if( k == mot) continue;
	if( eModel -> z[seq][k] == -1) continue;
	startMot1 = eModel -> z[seq][k];
	endMot1 = eModel -> z[seq][k] +  eModel -> motifs[k] -> w - 1;
	if(startMot <= endMot1 && startMot >= startMot1) return(2);
	if(endMot <= endMot1 && endMot >= startMot1) return(2);
	if(startMot1 <= endMot && startMot1 >= startMot) return(2);
	if(endMot1 <= endMot && endMot1 >= startMot) return(2);

    }
    if(X -> revFlag == 1)
    {
	for(k = 0; k< eModel -> numOfMotifs;k++)
	{
	    if( k == mot) continue;
	    if( eModel -> z[seq][k] == -1) continue;
	    endMot1 = reverseLocation(X, seq,eModel -> z[seq][k]);
	    startMot1 = endMot1 - eModel -> motifs[k] -> w;
	    endMot1 --;
	    if(startMot <= endMot1 && startMot >= startMot1) return(2);
	    if(endMot <= endMot1 && endMot >= startMot1) return(2);
	    if(startMot1 <= endMot && startMot1 >= startMot) return(2);
	    if(endMot1 <= endMot && endMot1 >= startMot) return(2);

	}
    }
    return(0);
}

int reverseLocation(FASTA *X, int seq, int pos)
{
    return(X -> length[seq] - pos);
}




void randomPermute(int *array,int n,unsigned int *seeds)
{
    int i;
    int j,temp;
    for(i = n-1; i>=0; --i)
    {
	j = rand_r(seeds) % (i+1);
	temp = array[i];
	array[i] = array[j];
	array[j] = temp;
    }
}


double scorePieceByTable(modelStruct *eModel, int mot, int seq, int pos, FASTA *X)
{
    double ans;
    ans = X -> tableOfBackground[seq][pos][eModel -> motifs[mot] -> w-1];
    if(ans == 0)
    {
	fprintf(stderr,"Something went wrong in computing the background score\n");
	printf("seq %d, pos %d,  back %f", seq,pos,ans);
	
	exit(0);
	    
	
    }
    
    return(ans);
}
