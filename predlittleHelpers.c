
#include "predmodel.h"
#include <getopt.h>


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

void freeOriginalModel(originalModel *eModel)
{
    int i,k;
    for(k = 0;k < eModel -> numOfMotifs; k++)
    {
        freePhi(eModel -> motifs[k]);
    }
    free(eModel -> motifs);
    free(eModel -> gammaFraction);
    for(i = 0; i <  eModel -> numOfRegProgs; i++)
    {
        free(eModel -> VAfraction[i]);
        free(eModel -> VPfraction[i]);
    }
    free(eModel -> VAfraction);
    free(eModel -> VPfraction);
    free(eModel);
}


void freePredModel(predSeqStruct *toReturn,originalModel *eModel )
{
    int i,k;
    for(k = 0;k < toReturn -> numOfSequences; k++)
    {
        for(i = 0;i < eModel -> numOfRegProgs; i++)
            free(toReturn -> Z[k][i]);
        free(toReturn -> Z[k]);
        free(toReturn -> likelihoods[k]);
        free(toReturn -> avglikelihoods[k]);
    }
    free(toReturn -> Z);
    free(toReturn -> likelihoods);
    free(toReturn -> avglikelihoods);
    free(toReturn -> backs);
    free(toReturn -> maxLikelihood);
    free(toReturn -> totMotifs);
    free(toReturn);
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


int piecePredOverlapsWith(FASTA *X,originalModel *eModel,predSeqStruct *toReturn,int seq,int mot, int prog, int *overlaps)
{
    int k;
    int startMot, endMot,startMot1,endMot1;
    int *done;
    int n;
    startMot = toReturn -> Z[seq][prog][mot];
    endMot = toReturn -> Z[seq][prog][mot] +  eModel -> motifs[mot] -> w;
    n = 0;
    done = ALLOC(sizeof(int)*eModel ->numOfMotifs);
    memset(done,0, sizeof(int ) * eModel -> numOfMotifs);
    memset(overlaps,0, sizeof(int ) * eModel -> numOfMotifs);

    for(k = 0; k< eModel -> numOfMotifs;k++)
    {
        if( k == mot) continue;
        if( toReturn -> Z[seq][prog][k] == -1) continue;
        startMot1 = toReturn -> Z[seq][prog][k];
        endMot1 = toReturn -> Z[seq][prog][k] +  eModel -> motifs[k] -> w;
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
            if( toReturn -> Z[seq][prog][k] == -1) continue;
            endMot1 = reverseLocation(X, seq,toReturn -> Z[seq][prog][k]);
            startMot1 = endMot1 - eModel -> motifs[k] -> w;
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


double scorePieceByTable(originalModel *eModel, int mot, int seq, int pos, FASTA *X)
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
