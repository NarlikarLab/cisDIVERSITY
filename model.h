	#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
 
#define PSEUDO 0.1  //0.5
#define PSEUDOGAMMA 1.0
#define PSEUDOP  0.1 //0.99 //0.9
#define PSEUDOA 0.1 //0.5
#define BACKORDER 2
#define HEADERLIMIT 10000
#define ZPROB 1.0 //0.3 //1 //0.5
#define REG 15
#define MOTIFS 30
#define ITER 5000
#define TRIALS 10
#define MINW 6
#define MAXW 200
#define INITW 12
#define EPS 0.00001
#define MINSITES 20
#define MINSEQS 20
#define SEED 1
#define IMPROVE 500
#define NUM_THREADS 2
typedef struct
{
    int order;
    double **probs;
} MMStruct;

typedef struct
{
    double pseudo;
    double pseudoGamma;
    double pseudoP;
    double pseudoA;
    int iter;
    int improve;
    int trials;
    int minw;
    int maxw;
    int initw;
    int minsites;
    int minseqs;
    int rev;
    int backorder;
    int r;
    int m;
    int verbose;
    int seed;
    int threads;
    char *basepath ;
    int entropy;
    double zprob;
} structForAllConstants;

typedef struct
{
    int order;
    long int **wordCounts;
} wordsStruct;



typedef struct
{
    int w;
    double **phi;
} phiStruct;

typedef struct
{
    phiStruct **motif;
    int number;
} pssms;


typedef struct
{
    int w;
    double **counts;
    double total;
} countStruct;



typedef struct
{
    int num;
    int *length;
    int *nonRepLength;
    int **seq;
    char **header;
    int revFlag;
    int **noN;
    double ***tableOfBackground;
} FASTA;


typedef struct
{
    int numOfSequences; //n
    int numOfRegProgs; //r
    int numOfMotifs; //m
    int minW;
    int maxW;
    int **z;
    double **zDenom;
    double totDenom;
    int *I;
    countStruct **motifs;
    int *totMotifs;
    double **VAcounts;
    double **VPcounts;
    double *gammaCounts;
    double *partZ;
} modelStruct;


 
typedef struct {
    int tStart;
    int tEnd;
    char *prefix;
    FASTA *X;
    structForAllConstants *allCons;
    double best;
    int bestRun;
    unsigned int randseed;
} threadStruct ;


//main


//initializations
int checkConsistency(FASTA *X, modelStruct *eModel, structForAllConstants *allCons);
void checkConsistencyAndPrint(FASTA *X, modelStruct *eModel, char *errorString, structForAllConstants *allCons);

void cleanUp(modelStruct *eModel, int initW,  structForAllConstants *allCons);
modelStruct *cleanUpModel(modelStruct *eModel,structForAllConstants *allCons);

void computeCounts(FASTA *X, modelStruct *eModel, structForAllConstants *allCons);
void computeCountsMotif(FASTA *X, modelStruct *eModel, int mot, structForAllConstants *allCons);

void computeVAVP(modelStruct *eModel,  structForAllConstants *allCons);
void computeCounts(FASTA *X, modelStruct *eModel, structForAllConstants *allCons);
void computePartZ(modelStruct *eModel,int mot, structForAllConstants *allCons);

void copyCounts(modelStruct *from,modelStruct *to);

void copyI(modelStruct *from,modelStruct *to);
void copyz(modelStruct *from,modelStruct *to);

double **getBackground(FASTA *X, structForAllConstants *allCons);
modelStruct *initializeNullModel(FASTA *X,structForAllConstants *allCons);

modelStruct *initializeRandomModel(FASTA *X, structForAllConstants *allCons,unsigned int *seeds);
void initializeRandomI(modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds);
void initializeRandomIspecial(modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds, int number);

void initializeRandomz(FASTA *X, modelStruct *eModel, structForAllConstants *allCons,unsigned int *seeds);

void justCopyModel(modelStruct *from, modelStruct *to);

modelStruct *makeACopyOfModel(modelStruct *from);

FASTA *readFASTA(char *filename,structForAllConstants *allCons);
pssms *readMotifs(char *filename, structForAllConstants *allCons);

void pickAz(FASTA *X, modelStruct *eModel, int seq, int mot,unsigned int *seeds);
int removeSmallProgs(modelStruct *eModel, int s,structForAllConstants *allCons,unsigned int *seeds);
int removeSmallProgsCompletely(modelStruct *eModel, int s, structForAllConstants *allCons,unsigned int *seeds);

void setAllConst(structForAllConstants *allCons);

void setBackTable(FASTA *X, double **backs, int maxW);


//littleHelpers
void *ALLOC(size_t size) ;
int compare ( const void *pa, const void *pb ) ;
void freeBacks(double **b, FASTA *X);
void freecountStruct(countStruct *motif);
void freeFASTA(FASTA *X);

void fileError(char *filename, char *errorString);
void freeMM(MMStruct *M);
void freeModel(modelStruct *eModel);

void freePhi(phiStruct *motif);
void freepssms(pssms *motifs);

void freeWords(wordsStruct *N);
long int getIndexInt(int *seq, int start, int length);
double logWithCheck(double x, float base);
double max(double a, double b);
int maxPieceWithoutN(int *seq, int pos, int order );
double min(double a, double b);
void normalizeMM(MMStruct *model);
void normalizeMMwithPseudoCount(MMStruct *model);
double normalizeWeights(double *weights, int n);
int numberOfCols(modelStruct *eModel,int mot, double cutoff);
int pieceContainsN(int *seq, int start, int width);
int pieceContainsNorOverlaps(FASTA *X,modelStruct *eModel,int seq,int mot);
int pieceContainsNquick(FASTA *X, int start, int width, int seq);
int pieceOverlapsWith(FASTA *X,modelStruct *eModel,int seq,int mot,  int *overlaps);

void randomPermute(int *array,int n,unsigned int *seeds);
void removeMotifsWithLowEntropy(modelStruct *eModel, structForAllConstants *allCons, FILE *fp);

void resetEmptyMotifWidth(modelStruct *eModel, int mot, structForAllConstants *allCons);

int reverseLocation(FASTA *X, int seq, int pos);
void selectionSort(int **array, int n);

void setMaxSequenceWithoutN(FASTA *X, int seq);
double scorePieceByTable(modelStruct *eModel, int mot, int seq, int pos, FASTA *X);


//probabilityFunctions
void addMotifCounts(FASTA *X,modelStruct *eModel, int mot, int seq,structForAllConstants *allCons);
int binarySearchSample(double value, double *weights, int start, int end);
int clubNonMotifs(FASTA *X, modelStruct *eModel,structForAllConstants *allCons,unsigned int *seeds);

int clubNonMotifsSeq(FASTA *X, modelStruct *eModel,int seq,structForAllConstants *allCons,unsigned int *seeds);
int clubNonMotifsTogether(FASTA *X, modelStruct *eModel,structForAllConstants *allCons);

double fullPosterior(FASTA *X, modelStruct *eModel,structForAllConstants *allCons);

int pickBestI(FASTA *X, modelStruct *eModel,  int seq,structForAllConstants *allCons);
int pickBestz(FASTA *X, modelStruct *eModel, int mot, int seq,structForAllConstants *allCons);

void removeMotifCounts(FASTA *X,modelStruct *eModel, int mot, int seq,structForAllConstants *allCons);
int sampleI(FASTA *X, modelStruct *eModel, int seq, double temp,unsigned int *seeds);
int sampleIlimitedProgs(FASTA *X, modelStruct *eModel,  int seq,  int regProg,unsigned int *seeds);

int sampleIOriginal(FASTA *X, modelStruct *eModel,  int seq, double **backs, double temps,structForAllConstants *allCons);

int sampleLeftEnd(FASTA *X, modelStruct *eModel, int mot, structForAllConstants *allCons,unsigned int *seeds);
int sampleRightEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons,unsigned int *seeds);


int samplez(FASTA *X, modelStruct *eModel, int mot, int seq, double temp);
double scorePieceByTable(modelStruct *eModel, int mot, int seq, int pos, FASTA *X);
double scorePieceByMotif(FASTA *X, modelStruct *eModel, int mot,int seq);
double scorePieceByMotifAO(FASTA *X, modelStruct *eModel, int mot,int seq);

double scorePieceByInitializationMotif(FASTA *X, modelStruct *eModel, int mot, phiStruct *motif, int seq, int pos);
int trimLeftEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons);
int trimRightEnd(FASTA *X, modelStruct *eModel, int mot,structForAllConstants *allCons);


//model
modelStruct *learnModel(FASTA *X, modelStruct *eModel, FILE *fp1,FILE *fp2,structForAllConstants *allCons,unsigned int *seeds);


//printing
char DNA(int i);
void printAllConsensus(modelStruct *eModel, FILE *fp,structForAllConstants *allCons);

void printAllMotifs(modelStruct *eModel, FILE *fp,structForAllConstants *allCons);
void printConsensus(countStruct *motif, FILE *fp);

void printError(char *str);
void printinfo(FASTA *X,modelStruct *eModel, FILE *fp,structForAllConstants *allCons);

void printMotif(countStruct *motif, FILE *fp,structForAllConstants *allCons);
void printRevSites(FASTA *X, modelStruct *eModel, FILE **fp);
void printSites(FASTA *X, modelStruct *eModel, FILE **fp);

void printv(modelStruct *eModel, FILE *fp,structForAllConstants *allCons);
 

//sampling
int samplezAO(FASTA *X, modelStruct *eModel, int mot, int seq, double temp,structForAllConstants *allCons,unsigned int *seeds);


