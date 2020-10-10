#include "model.h"
#include <getopt.h>
#include <pthread.h>
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>

char fastafile[1000] = {'\0'};
double globalBest = -DBL_MAX;
int globalBestRun = -100;
void printUsage(char *a)
{
     fprintf(stderr, "usage: %s  [options]  \n",a);
     fprintf(stderr, "       -f inputFastaFile (compulsory)\n");
     fprintf(stderr, "       -o outputDir (will create it, compulsory)\n");
     fprintf(stderr, "       -h (to print full set of options)\n");

     exit(0);
}

void printFullUsage(char *a)
{
     fprintf(stderr, "usage: %s  [options]  \n",a);
     fprintf(stderr, "       -f input fasta file (compulsory)\n");
     fprintf(stderr, "       -o output directory (will create it, compulsory)\n\n");
     fprintf(stderr, "Change these if you think you have more than default modules/motifs in your data:\n");
     fprintf(stderr, "       -r maximum number of regulatory programs/modules (default: %d)\n",REG);
     fprintf(stderr, "       -m maximum number of motifs (default: %d)\n\n",MOTIFS);
     fprintf(stderr, "Change these if you want a quick (& dirty) or slow (& better) solution:\n");
     fprintf(stderr, "       -i maximum number of iterations per run (default: %d)\n",ITER);
     fprintf(stderr, "       -t number of distinct runs (default: %d)\n",TRIALS);
     fprintf(stderr, "       -n number of processors (cores) to use (default: %d)\n",NUM_THREADS);
     fprintf(stderr, "       -v 0,1,2 (verbosity: 0 to print little, 1 for progress, which is the default, and 2 for that and log likelihoods at each iteration -- can create large files)\n\n");

     fprintf(stderr, "Change these if you want to change the initialization/hyperparameters of the model:\n");
     fprintf(stderr, "       -a random seed (default: %d)\n",SEED);
     fprintf(stderr, "       -s starting width of each motif (default: %d)\n",INITW);
     fprintf(stderr, "       -p pseudocount for the PWMs (default: %lf)\n",PSEUDO);
     fprintf(stderr, "       -g pseudocount for the modules (default: %lf)\n",PSEUDOGAMMA);
     fprintf(stderr, "       -Y pseudocount for a modules to have a PWM (default: %lf)\n",PSEUDOP);
     fprintf(stderr, "       -N pseudocount for a modules to NOT have a PWM (default: %lf)\n\n",PSEUDOA);

     fprintf(stderr, "Change these only if you read the paper carefully and you know what you are doing:\n");

//     fprintf(stderr, "       -c number of iterations for no improvement (default: %d)\n",IMPROVE);
     fprintf(stderr, "       -w minimum width of a motif (default: %d)\n",MINW);
     fprintf(stderr, "       -W maximum width of a motif (default: %d)\n",MAXW);
     fprintf(stderr, "       -S minimum number of sites to create a motif (default: %d)\n",MINSITES);
     fprintf(stderr, "       -C minimum number of sequences in a regulatory program/module (default: %d)\n",MINSEQS);
     fprintf(stderr, "       -e (if low information sequence signatures to be removed; default is to keep all)\n");
     fprintf(stderr, "       -R (if motifs should be searched for on a same strand; default: allow reverse complements)\n");    
     fprintf(stderr, "       -b order of the Markov model for the background sequences (default: %d)\n",BACKORDER);
//     fprintf(stderr, "       -z initial probability of presence of a motif in a sequence (default: %lf)\n",ZPROB);
}

void initialChecks(structForAllConstants *allCons)
{
    if(allCons -> minw > allCons -> maxw)
    {
	if(allCons -> maxw == MAXW)
	{
	    fprintf(stderr,"minimum width of motif %d is too large, making it %d\n", allCons -> minw,MAXW);
	    allCons -> minw = MAXW;
	}
	else
	{
	    if(allCons -> minw == MINW)
	    {
		fprintf(stderr,"maximum width of motif %d is too small, making it %d\n", allCons -> maxw,MINW);
		allCons -> maxw = MINW;
	    }
	    else
	    {
		fprintf(stderr,"minimum and maximum width of motifs (%d, %d) are weird, will reset to default: %d, %d\n", allCons -> minw, allCons -> maxw,MINW,MAXW);
		allCons -> minw = MINW;
		allCons -> maxw = MAXW;
	    }
	}
    }
    if(allCons -> trials < allCons -> threads)
	allCons -> threads =  allCons -> trials;

    if(allCons -> maxw < allCons -> initw || allCons -> initw < allCons -> minw)
	allCons -> initw = (allCons -> maxw + allCons -> minw)/2;

}
void printSettings(FILE *fp, char *filename,structForAllConstants *allCons)
{
    fprintf(fp, "Name of input file: %s\n", filename);
    fprintf(fp, "Maximum number of regulatory programs: %d\n", allCons->r);
    fprintf(fp, "Maximum number of motifs: %d\n", allCons -> m);
    fprintf(fp, "Maximum number of iterations per run: %d\n", allCons -> iter);
    fprintf(fp, "Number of distinct runs: %d\n", allCons -> trials);
    fprintf(fp, "Starting width of each motif : %d\n", allCons -> initw);
    fprintf(fp, "Minimum width of a motif: %d\n", allCons -> minw);
    fprintf(fp, "Maximum width of a motif: %d\n", allCons -> maxw);
    fprintf(fp, "Initial probability of a motif in a sequence: %lf\n", allCons -> zprob);
    if(allCons -> verbose == 2)
	fprintf(fp, "Will print all log posteriors in a file called posteriors in each run directory  (caution: may create huge files)\n");
    if(allCons -> verbose == 0)
	fprintf(fp, "Will not print much\n");
    if(allCons -> verbose >= 1)
	fprintf(fp, "Will print progress in a file called progress in each run directory \n");
    if(allCons -> rev == 1)
	fprintf(fp, "Will look on both strands\n");
    if(allCons -> rev == 0)
	fprintf(fp, "Will look only on one strand\n");
    fprintf(fp, "Minimum number of sites to create a motif: %d\n", allCons -> minsites);
    fprintf(fp, "Minimum number of sequences to create a regulatory program: %d\n", allCons -> minseqs);
    fprintf(fp, "Pseudocount for the PWMs: %lf\n",allCons -> pseudo);
    fprintf(fp, "Pseudocount for the regulatory programs: %lf\n",allCons -> pseudoGamma);
    fprintf(fp, "Pseudocount for module to have a PWM: %lf\n",allCons -> pseudoP);
    fprintf(fp, "Pseudocount for module to NOT have a PWM: %lf\n",allCons -> pseudoA);
    fprintf(fp, "Order of the background Markov model: %d\n",allCons -> backorder);
    fprintf(fp, "Likelihood must not improve for %d iterations (convergence criterion)\n",allCons -> improve);
    fprintf(fp, "Will use %d cores\n",allCons -> threads);
    fprintf(fp, "Will use %d as the primary seed for the random number generator\n",allCons -> seed);
    if(allCons -> entropy == 1)
	fprintf(fp, "Will remove motifs with low information content\n");
    else
	fprintf(fp, "Will NOT remove motifs with low information content\n");
}

void parseOpts(int argc, char **argv,  char *f,char *o,structForAllConstants *allCons)
{
    int opt;
    int a;
    double b;
    while((opt = getopt(argc, argv, "f:o:r:m:i:t:w:W:s:v:RM:a:P:hep:g:b:S:C:c:Y:N:n:x:y:z:")) > 0)
    {
	switch(opt) 
        {
	case 'z':
	    allCons -> zprob = atof(optarg);
            break;
	case 'f':
	    strcpy(f,optarg);
            break;
	case 'o':
	    sprintf(o,"%s/",optarg);
	    break;
	case 'e':
	    allCons -> entropy = 1;
	    break;
	case 'r':
	    a = atoi(optarg);
	    if(a < 1)
		fprintf(stderr,"Number of regulatory programs %d not allowed, using default of %d\n", a,allCons -> r);
	    else
		allCons -> r = a;
	    break;
	case 'm':
	    a = atoi(optarg);
	    if(a < 1)
		fprintf(stderr,"Number of motifs %d not allowed, using default of %d\n", a, allCons -> m);
	    else
		allCons -> m = a;
	    break;
	case 'i':
	    a = atoi(optarg);
	    if(a < 0)
		fprintf(stderr,"Number of iterations %d not allowed, using default of %d\n", a, allCons -> iter);
	    else
		allCons -> iter = a;
	    break;
	case 't':
	    a = atoi(optarg);
	    if(a < 1)
		fprintf(stderr,"Number of runs %d not allowed, using default of %d\n", a,allCons -> trials);
	    else
		allCons -> trials = a;
	    break;
	case 'w':
	    a = atoi(optarg);
	    if(a < allCons -> minw | a > allCons -> maxw)
		fprintf(stderr,"Minimum width %d not allowed, using default of %d\n", a,allCons -> minw);
	    else
		allCons -> minw = a;
	    break;
	case 'W':
	    a = atoi(optarg);
	    if(a > allCons -> maxw |a < allCons -> minw )
		fprintf(stderr,"Maximum width %d not allowed, using default of %d\n", a, allCons -> maxw);
	    else
		allCons -> maxw = a;
	    break;
	case 'S':
	    a = atoi(optarg);
	    allCons -> minsites = a;
	    break;
	case 'C':
	    a = atoi(optarg);
	    allCons -> minseqs = a;
	    break;
	case 's':
	    a = atoi(optarg);
	    if(a < allCons -> minw  | a > allCons -> maxw)
		fprintf(stderr,"Starting width %d not allowed, using default of %d\n", a,allCons -> initw);
	    else
		allCons -> initw = a;
	    break;
	case 'v':
	    a = atoi(optarg);
	    if(a == 0 || a == 1 || a== 2)
		allCons -> verbose = a;
	    else
		allCons -> verbose = 1;
	    break;
	case 'R':
	    allCons -> rev = 0;
	    break;
	case 'c':
	    a = atoi(optarg);
	    allCons -> improve = a;
	    break;
	case 'n':
	    a = atoi(optarg);
	    if(a < 1 )
	    {
		fprintf(stderr,"Number of processors should be at least 1, setting it to 1\n");
		a = 1;
	    }
	    allCons -> threads = a;
	    break;
	case 'b':
	    allCons -> backorder = atoi(optarg);
	    break;

	case 'a':
	    allCons -> seed = atoi(optarg);
	    break;
	case 'p':
	    allCons -> pseudo = atof(optarg);
	    break;
	case 'Y':
	    allCons -> pseudoP = atof(optarg);
	    break;
	case 'N':
	    allCons -> pseudoA = atof(optarg);
	    break;
	case 'g':
	    allCons -> pseudoGamma = atof(optarg);
	    break;
	case 'h':
	    printFullUsage(argv[0]);
	    exit(0);
	    break;
	default:
            printUsage(argv[0]);
            exit(0);
	}
    }
}


void *forThreading(void *arguments)
{
    threadStruct *args = arguments;
    char trialprefix[1000] = {'\0'};
    char likefile[1000] = {'\0'};
    char infofile[1000] = {'\0'};
    char progressfile[1000] = {'\0'};
    char pssmfile[1000] = {'\0'};
    char command[10000] = {'\0'};
    modelStruct *M, *learned, *cleanedUp;
    FILE *info, *pssm, *sets, *prog;
    FILE **sites;
    double current;
    FILE *like;
    int k;
    char sitesfiles[100][1000];
    int i;
    unsigned int *seeds;
    info = NULL;
    pssm = NULL;
    sets = NULL;
    prog= NULL;
    seeds = ALLOC(sizeof(unsigned int) * (args -> tEnd - args -> tStart));
    for(i = args -> tStart; i< args -> tEnd; i++)
    {
	*seeds = i *args -> allCons -> seed ;
	sprintf(trialprefix,"%s/run%d",args -> prefix,i );
	printf("run %d, seed %d\n",i, *seeds);
	fflush(stdout);
	if(mkdir(trialprefix,0777) == -1)
	    printf("run %d directory already there\n",i);
 
	
	M = initializeRandomModel(args -> X,args ->allCons,seeds);

//	printf("Initial full posterior %lf at %d\n", fullPosterior(args -> X, M, args -> allCons),i);
	
	if(args -> allCons -> verbose == 2)
	{
	    sprintf(likefile,"%s/posteriors.txt",trialprefix);
	    if( NULL == (like = fopen(likefile, "w")))
		fileError(likefile, "is not writable");
	}
	else
	    like = NULL;
	if(args -> allCons -> verbose >= 1)
	{
	    sprintf(progressfile,"%s/progress.txt",trialprefix);
	    if( NULL == (prog = fopen(progressfile, "w")))
		fileError(progressfile, "is not writable");
	}
	else
	    like = NULL;
	sprintf(infofile,"%s/info.txt",trialprefix);
	if( NULL == (info = fopen(infofile, "w")))
	    fileError(infofile, "is not writable");
	sprintf(pssmfile,"%s/pssm.txt",trialprefix);
	if( NULL == (pssm = fopen(pssmfile, "w")))
	    fileError(pssmfile, "is not writable");

	learned = learnModel(args -> X, M, like,prog,args ->allCons,seeds);

	cleanedUp = cleanUpModel(learned,args ->allCons);
	
//	cleanedUp = makeACopyOfModel(learned);
	checkConsistencyAndPrint(args ->X, cleanedUp,"quitting at cleanedup",args ->allCons);

	freeModel(learned);

	printinfo(args -> X,cleanedUp,info,args -> allCons); 
	printv(cleanedUp, pssm,args ->allCons);
	printAllMotifs(cleanedUp, pssm,args ->allCons);
	printAllConsensus(cleanedUp, pssm,args ->allCons);
    
	current = fullPosterior(args -> X, cleanedUp,args ->allCons);
	fprintf(pssm,"\nFull posterior for run %d is %lf\n", i,current);
	if(i == args -> tStart)
	{
	    args -> best = current;
	    args -> bestRun = i;
	}
	else
	    if(current > args -> best)
	    {
		args -> best = current;
		args -> bestRun = i;
	    }
	if(globalBest < args -> best)
	{
	    globalBest = args -> best;
	    globalBestRun = i;
	} 
	if(args -> allCons -> verbose == 2)
	{
	    fclose(like);
	}

	fclose(pssm);
	fclose(info);
	sites = ALLOC(sizeof(FILE *) * cleanedUp -> numOfMotifs);
	memset(sites ,0,sizeof(FILE * ) * cleanedUp -> numOfMotifs);
    
	for(k =0 ;k< cleanedUp -> numOfMotifs; k++)
	{
	    if(fabs(cleanedUp -> motifs[k] -> total - args -> allCons -> pseudo * 4) < EPS) continue;
	    if(k < 9)
		sprintf(sitesfiles[k],"%s/sites_0%d.txt",trialprefix, k+1);
	    else
		sprintf(sitesfiles[k],"%s/sites_%d.txt",trialprefix, k+1);
	    if( NULL == (sites[k] = fopen(sitesfiles[k], "w")))
		fileError(infofile, "is not writable");
	}
	printSites(args -> X, cleanedUp,sites );
	for(k =0 ;k< cleanedUp -> numOfMotifs; k++)
	{
	    if(! sites[k] == 0)
		fclose(sites[k]);
	}
	memset(sites ,0,sizeof(FILE * ) * cleanedUp -> numOfMotifs);
	
	for(k =0 ;k< cleanedUp -> numOfMotifs; k++)
	{
	    if(fabs(cleanedUp -> motifs[k] -> total - args -> allCons -> pseudo * 4) < EPS) continue;
	    if(k < 9)
		sprintf(sitesfiles[k],"%s/revsites_0%d.txt",trialprefix, k+1);
	    else
		sprintf(sitesfiles[k],"%s/revsites_%d.txt",trialprefix, k+1);
	    
	    if( NULL == (sites[k] = fopen(sitesfiles[k], "w")))
		fileError(sitesfiles[k], "is not writable");
	}
	printRevSites(args -> X, cleanedUp,sites );
	for(k =0 ;k< cleanedUp -> numOfMotifs; k++)
	{
	    if(! sites[k] == 0)
		fclose(sites[k]);
	}
	free(sites);
	freeModel(cleanedUp);
	freeModel(M);
	sprintf(command,"python %s/createLogos.py %s %s",args -> allCons -> basepath,trialprefix, trialprefix);
	sprintf(command,"python %s/createLogosPNG.py %s %s",args -> allCons -> basepath,trialprefix, trialprefix);
	if(args -> allCons -> verbose >= 1) fprintf(prog,"%s\n",command);
	system(command);
	sprintf(command, "python %s/postProcessNoOptions.py %s %s %s %s ", args -> allCons -> basepath, fastafile, infofile, pssmfile, trialprefix);
	    
	if(args -> allCons -> verbose >= 1) fprintf(prog,"%s\n",command);
	system(command);
	sprintf(command, "Rscript  %s/createPlots.r %s",args -> allCons -> basepath,trialprefix);
	if(args -> allCons -> verbose >= 1) fprintf(prog,"%s\n",command);
	system(command);
	sprintf(command, "python  %s/createData.py %s %s %s",args -> allCons -> basepath,trialprefix,fastafile,args -> allCons -> basepath);
	if(args -> allCons -> verbose >= 1) fprintf(prog,"%s\n",command);
	system(command);
	if(args -> allCons -> verbose >= 1)
	{
	    fclose(prog);
	}
	seeds++;
    }

}

int main(int argc, char **argv)
{
    FASTA *X;
    
    double **backProbs;
    int i;
    char outputprefix[1000] = {'\0'};
    char infofile[1000] = {'\0'};
    char pssmfile[1000] = {'\0'};
    char setsfile[1000] = {'\0'};
    char command[1000];
    FILE  *source;
    FILE *info, *pssm, *sets;
    double current;
    structForAllConstants *allCons;
    pthread_t some_thread;
    threadStruct **argsForThread; //[NUM_THREADS];
    pthread_t *threads; //[ NUM_THREADS ];
    int result_code;
    int tasks, threadCount,tStart;
    double best = -DBL_MAX;
    int bestRun;
    modelStruct *M;
    
    char *lastslash=NULL;

    allCons = ALLOC(sizeof(structForAllConstants));

    setAllConst(allCons);

    allCons -> basepath = strdup(argv[0]);
    lastslash = strrchr(allCons -> basepath,'/');
    if (lastslash != NULL)
    {
	*lastslash = 0; // chop name off of string
    }
    else
    {
	allCons -> basepath = 0; // path is blank string
    }
    printf("path is %s\n", allCons -> basepath);
    parseOpts(argc, argv, fastafile, outputprefix, allCons);
    if(strlen(fastafile) == 0 || strlen(outputprefix) == 0 ) printUsage(argv[0]);
    if( NULL == (source = fopen(fastafile, "r")))
	fileError(fastafile, "is not readable");
    fclose(source);
    initialChecks(allCons);
    if(mkdir(outputprefix,0777) == -1) fileError(outputprefix, " directory cannot be created");
    sprintf(setsfile,"%ssettings.txt",outputprefix);
    if( NULL == (sets = fopen(setsfile, "w")))
	fileError(setsfile, "is not writable");
    printSettings(sets, fastafile, allCons);
    fclose(sets);
    X = readFASTA(fastafile,allCons);
    backProbs  = getBackground(X,allCons);
    setBackTable(X, backProbs,allCons -> maxw);

    M = initializeNullModel(X, allCons);
//    fprintf(stdout,"empty model posterior is %lf", fullPosterior(X, M,allCons));
//    fflush(stdout);
    freeModel(M);
    tasks =  allCons -> trials;
    threadCount = allCons -> threads;
    tStart = 0;
    argsForThread = ALLOC(sizeof(threadStruct *) *allCons -> threads);
    threads = ALLOC(sizeof(pthread_t) * allCons -> threads);
    for(i = 0; i< allCons -> threads; i++)
    {
	if(tasks <= 0)
	{
	    fprintf(stderr,"something wrong in bunching tasks\n");
	    exit(0);
	}
	argsForThread[i] = ALLOC(sizeof(threadStruct));
	
	argsForThread[i] -> prefix = outputprefix;
	argsForThread[i] -> X = X;
	argsForThread[i] -> allCons = allCons;
	argsForThread[i] -> tStart = tStart;
	argsForThread[i] -> tEnd = tStart + ceil(tasks/(0.0+threadCount));
	tStart = argsForThread[i] -> tEnd ;
	tasks = tasks -  ceil(tasks/(0.0+threadCount));
	
	threadCount--;
	result_code = pthread_create(&threads[i], NULL, &forThreading, argsForThread[i]);
	assert( !result_code );
	printf("Thread %d: %d to %d\n", i, argsForThread[i] -> tStart, argsForThread[i] -> tEnd);
    }
    for( i = 0; i< allCons -> threads; i++)
    {
	// block until thread 'index' completes
	result_code = pthread_join( threads[ i ], NULL );
	assert( !result_code );
	printf( "In main: thread %d has completed\n", i );

	if(best < argsForThread[i] -> best)
	{
	    best = argsForThread[i] -> best;
	    bestRun = argsForThread[i] -> bestRun;
	}
	printf( "globalbest is %lf at run %d\n", globalBest,globalBestRun );
    
    }

    sprintf(command," mv %s/run%d %s/bestSolution ", outputprefix,globalBestRun,outputprefix);
//    printf("%s",command);
    system(command);
    fprintf(stdout,"\nDone; exitting.\nOutput is in %s/bestSolution\n",outputprefix);
}
