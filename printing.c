
#include "model.h"
#include <getopt.h>



char DNA(int i)
{
    switch (i)
    {
    case 0:
	return('A');
    case 1:
	return('C');
    case 2:
	return('G');
    case 3:
	return('T');
    default:
	return('N');
    }

}

void printError(char *str)
{
    fprintf(stderr, "%s  \n",str);
    fprintf(stderr, "Exiting...  \n");
    exit(0);

}

void printRevSites(FASTA *X, modelStruct *eModel, FILE **fp)
{
    int i,j,k;
    char dna[4];
    dna[0] = 'T';
    dna[1] = 'G';
    dna[2] = 'C';
    dna[3] = 'A';
    
    for(i =0; i<eModel ->numOfSequences ; i++)
    {
	for(k =0; k< eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[i][k] == -1) continue;
	    for(j= eModel -> motifs[k] -> w -1; j>=0; j--)
	    {
		fprintf(fp[k], "%c", dna[X -> seq[i][eModel -> z[i][k]+j]]);
	    }
	    fprintf(fp[k],"\n");
	}
		
    }
}

void printSites(FASTA *X, modelStruct *eModel, FILE **fp)
{
    int i,j,k;
    char dna[4];
    dna[0] = 'A';
    dna[1] = 'C';
    dna[2] = 'G';
    dna[3] = 'T';
    
    for(i =0; i<eModel ->numOfSequences ; i++)
    {
	for(k =0; k< eModel -> numOfMotifs; k++)
	{
	    if(eModel -> z[i][k] == -1) continue;
	    for(j=0; j< eModel -> motifs[k] -> w; j++)
	    {
		fprintf(fp[k], "%c", dna[X -> seq[i][eModel -> z[i][k]+j]]);
	    }
	    fprintf(fp[k],"\n");
	}
		
    }
}

void printConsensus(countStruct *motif, FILE *fp)
{
    int j;
    double A, C, G,T;
    for(j=0;j<motif -> w;j++)
    {
	A = (double) motif -> counts[0][j]/motif -> total;
	C = (double) motif -> counts[1][j]/motif -> total;
	G = (double) motif -> counts[2][j]/motif -> total;
	T = (double) motif -> counts[3][j]/motif -> total;
	if(A > 0.8)
	{
	    fprintf(fp, "A");
	    continue;
	}
	if(C > 0.8)
	{
	    fprintf(fp, "C");
	    continue;
	}
	if(G > 0.8)
	{
	    fprintf(fp, "G");
	    continue;
	}
	if(T > 0.8)
	{
	    fprintf(fp, "T");
	    continue;
	}
	if(A > 0.5)
	{
	    fprintf(fp, "a");
	    continue;
	}
	if(C > 0.5)
	{
	    fprintf(fp, "c");
	    continue;
	}
	if(G > 0.5)
	{
	    fprintf(fp, "g");
	    continue;
	}
	if(T > 0.5)
	{
	    fprintf(fp, "t");
	    continue;
	}
	if(A + C > 0.9)
	{
	    fprintf(fp, "M");
	    continue;
	}
	if(A + G > 0.9)
	{
	    fprintf(fp, "R");
	    continue;
	}
	if(A + T > 0.9)
	{
	    fprintf(fp, "W");
	    continue;
	}
	if(C + G > 0.9)
	{
	    fprintf(fp, "S");
	    continue;
	}
	if(C + T > 0.9)
	{
	    fprintf(fp, "Y");
	    continue;
	}
	if(G + T > 0.9)
	{
	    fprintf(fp, "K");
	    continue;
	}
	fprintf(fp, "N");
	
    }
    //   fprintf(fp, "\n");
}



void printRevConsensus(countStruct *motif, FILE *fp)
{
    int j;
    double A, C, G,T;
    for(j=motif -> w-1; j >= 0 ;j--)
    {
	A = (double) motif -> counts[3][j]/motif -> total;
	C = (double) motif -> counts[2][j]/motif -> total;
	G = (double) motif -> counts[1][j]/motif -> total;
	T = (double) motif -> counts[0][j]/motif -> total;
	if(A > 0.8)
	{
	    fprintf(fp, "A");
	    continue;
	}
	if(C > 0.8)
	{
	    fprintf(fp, "C");
	    continue;
	}
	if(G > 0.8)
	{
	    fprintf(fp, "G");
	    continue;
	}
	if(T > 0.8)
	{
	    fprintf(fp, "T");
	    continue;
	}
	if(A > 0.5)
	{
	    fprintf(fp, "a");
	    continue;
	}
	if(C > 0.5)
	{
	    fprintf(fp, "c");
	    continue;
	}
	if(G > 0.5)
	{
	    fprintf(fp, "g");
	    continue;
	}
	if(T > 0.5)
	{
	    fprintf(fp, "t");
	    continue;
	}
	if(A + C > 0.9)
	{
	    fprintf(fp, "M");
	    continue;
	}
	if(A + G > 0.9)
	{
	    fprintf(fp, "R");
	    continue;
	}
	if(A + T > 0.9)
	{
	    fprintf(fp, "W");
	    continue;
	}
	if(C + G > 0.9)
	{
	    fprintf(fp, "S");
	    continue;
	}
	if(C + T > 0.9)
	{
	    fprintf(fp, "Y");
	    continue;
	}
	if(G + T > 0.9)
	{
	    fprintf(fp, "K");
	    continue;
	}
	fprintf(fp, "N");
	
    }
//    fprintf(fp, "\n");
}

void printMotif(countStruct *motif, FILE *fp,structForAllConstants *allCons)
{
    int i, j;
    for(i=0;i<4;i++)
    { 
	fprintf(fp, "%c\t[ ", DNA(i));
	for(j=0;j<motif -> w;j++)
	{
	    fprintf(fp,"%0.1f\t", motif -> counts[i][j] - allCons -> pseudo);
	}
	fprintf(fp, "]\n");
    }
}


void printAllMotifs(modelStruct *eModel, FILE *fp,structForAllConstants *allCons)
{
    int k;
    for(k =0;k< eModel -> numOfMotifs;k++)
    {
	if(fabs(eModel -> motifs[k] -> total - allCons -> pseudo *4) < EPS) continue;
	fprintf(fp,"Motif #%d\n", k+1);
	fprintf(fp,"total sequences %0.0f\n", eModel -> motifs[k] -> total -  allCons -> pseudo *4);
	printMotif(eModel -> motifs[k],fp,allCons);
    }
}

void printAllConsensus(modelStruct *eModel, FILE *fp,structForAllConstants *allCons)
{
    int k;
    for(k=0;k < eModel -> numOfMotifs;k++)
    {
	if(fabs(eModel -> motifs[k] -> total - allCons -> pseudo *4) <= EPS) continue;
	fprintf(fp,"Motif #%d; ", k+1);
	fprintf(fp,"total sequences %0.0f; ", eModel -> motifs[k] -> total -allCons -> pseudo *4);
	printConsensus(eModel-> motifs[k],fp);
	fprintf(fp," / " );
	printRevConsensus(eModel-> motifs[k],fp);
	fprintf(fp,"\n");
	
    }
}

void printinfo(FASTA *X,modelStruct *eModel, FILE *fp,structForAllConstants *allCons)
{
    int i,k;
    for(i =0; i<eModel ->numOfSequences ; i++)
    {
	fprintf(fp,"%d\t%s\t%d\t",i+1, X -> header[i], eModel -> I[i]+1);
	for(k =0; k< eModel -> numOfMotifs; k++)
	{
	    if(fabs(eModel -> motifs[k] -> total - allCons -> pseudo *4) <= EPS) continue;
	    if(eModel -> z[i][k] == -1)
	    {
		fprintf(fp,"NA\t");
		continue;
	    }
	    if( eModel -> z[i][k] >= X -> length[i]/2 && X -> revFlag ==1)
	    {
		fprintf(fp, "-%d\t", reverseLocation(X,i,eModel -> z[i][k]) - eModel -> motifs[k] -> w );
	    }
	    else	
		fprintf(fp, "%d\t", eModel -> z[i][k]);
	}
	fprintf(fp,"\n");
		
    }
}

void printv(modelStruct *eModel, FILE *fp,structForAllConstants *allCons )
{
    int j,k;
    int r,m;
    r = 0;
    m = 0;
    fprintf(fp, "Number of sequences: %d\n", eModel -> numOfSequences);

    for(j=0;j< eModel->numOfRegProgs; j++)
    {
	if(fabs(eModel -> gammaCounts[j] - allCons -> pseudoGamma) <= EPS) continue;
	r++;
    }
    for(k=0;k< eModel ->numOfMotifs; k++)
    {
	if(fabs(eModel -> motifs[k] -> total  - allCons -> pseudo *4) <= EPS) continue;
	m++;
    }
    
    fprintf(fp, "Number of programs: %d\n", r);
    fprintf(fp, "Number of motifs: %d\n", m);
    fprintf(fp, "Motif widths: ");
    for(k=0;k< eModel ->numOfMotifs; k++)
    {
	if(fabs(eModel -> motifs[k] -> total  - allCons -> pseudo *4) <= EPS) continue;
	fprintf(fp, "%d ", eModel -> motifs[k] -> w);
    }
    fprintf(fp, "\n");
    for(j=0;j< eModel->numOfRegProgs; j++)
    {
	if(fabs(eModel -> gammaCounts[j] - allCons -> pseudoGamma) <= EPS) continue;
	fprintf(fp,"Program %d: total %0.0f\n",j+1 ,eModel -> gammaCounts[j]  - allCons -> pseudoGamma);			    
	for(k=0;k< eModel ->numOfMotifs; k++)
	{
	    if(fabs(eModel -> VPcounts[j][k] - allCons -> pseudoP) <= EPS) continue;
	    fprintf(fp, " motif #%d: %0.0f of %0.0f (%f posterior)\n", k+1, eModel -> VPcounts[j][k] - allCons -> pseudoP, ( eModel ->VPcounts[j][k] +  eModel ->VAcounts[j][k] -allCons -> pseudoP - allCons -> pseudoA ), eModel -> VPcounts[j][k]/ ( eModel ->VPcounts[j][k] +  eModel ->VAcounts[j][k] )); 
	}
	fprintf(fp, "\n"); 
    }
}
