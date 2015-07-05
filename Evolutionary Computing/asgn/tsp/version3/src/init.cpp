#include "stdafx.hpp"
#include <string.h>
#include <algorithm>
#include "type.hpp"

using namespace std;
void   statistics(POPULATION *p, IPTR pop);
void   rawStat(FILE *fp, POPULATION *p, IPTR pop);
double eval(POPULATION *p, IPTR pi);
void   initData(char *inputFile, POPULATION *p);
void   initPop(POPULATION *p);
void   initReport(POPULATION *p);

void initialize(char *argv[], POPULATION *p)
{
    /* initialize everything */
    char *Ifile = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
    strcpy(Ifile, argv[1]);
    initData(Ifile, p);
    printf("%s data loaded.\n", __PROMT_PASSED);
    initPop(p);
    printf("%s population initialized.\n", __PROMT_PASSED);
    statistics(p,p->op);
    initReport(p);
}

void initData(char *Ifile, POPULATION *p)
{
	//inittialize global params
    FILE *inpfl;
    inpfl = fopen(Ifile,"r");

    if( (inpfl = fopen(Ifile,"r")) == NULL){ printf("error in opening file %s \n", Ifile); exit(1);}
    fscanf(inpfl,"%d",&p->popSize);
    if(p->popSize % 2 != 0) p->popSize++;
    fscanf(inpfl,"%d",&p->lchrom);
    fscanf(inpfl,"%d",&p->maxGen);
    fscanf(inpfl,"%lf",&p->pCross);
    fscanf(inpfl,"%lf",&p->pMut);
    char tmp[1024];
    fscanf(inpfl,"%s", tmp);
    p->ofile = (char *) calloc (strlen(tmp)+1, sizeof(char));
	strcpy(p->ofile, tmp);
    std::remove(p->ofile);
    for(size_t row = 0; row < (size_t)p->lchrom; row++)
	{
        int elem;
        vector<int> r;
        for(size_t column = 0; column < (size_t)p->lchrom; column++)
        { fscanf(inpfl,"%i", &elem); r.push_back(elem); }
        p->city_costs.push_back(r);

    }
    fclose(inpfl);
    /* set progress indicators to zero */
    p->highestEverFitness = 0.0;
    p->highestEverGen = 0;
    p->highestEverIndex = 0;
}

void initPop(POPULATION *p)
{
	/* initialize a random population */
    IPTR pi, pj;
    //srand(time(0));
    p->op = (IPTR) calloc (p->popSize, sizeof(INDIVIDUAL));
    p->np = (IPTR) calloc (p->popSize, sizeof(INDIVIDUAL));
    for (int i = 0; i < p->popSize; i++){
        pi = &(p->op[i]);
        pi->chrom = (int *) calloc (p->lchrom, sizeof(int));
        pj = &(p->np[i]);
        pj->chrom = (int *) calloc (p->lchrom, sizeof(int));
        for (int j = 0; j<p->lchrom; j++) pi->chrom[j] = j;
        std::random_shuffle(pi->chrom, pi->chrom + p->lchrom);
        pi->fitness = eval(p, pi);
    }
}

void initReport(POPULATION* p)
{
    FILE *fp;
    printf("Population Size             ( popsize )    %d\n", p->popSize);
    printf("Chromosome length           ( lchrom  )    %d\n", p->lchrom);
    printf("Crossover Probability       ( pcross  )    %lf\n", p->pCross);
    printf("Mutation Probability        ( pmut    )    %lf\n", p->pMut);
    printf("Maximum num of generations  ( maxgen  )    %d\n\n\n", p->maxGen);
    if((fp = fopen(p->ofile, "w+")) == NULL)
    { fprintf(stderr, "error in opening file %s\n", p->ofile); exit(EXIT_FAILURE); }
#define XEF_HEADER " ", "HEG.", " ", "HEF.", " ", "HEI.", " "
#define header_format "%-10s%s%-7s%s%-8s%s%-8s%s%-10s%s%-11s%s%-13s%s%-12s\n", \
    " ", "gen#", " ", "max. fitness", " ", "avg. fitness", " ", "min. fitness", \
    XEF_HEADER
    // this header is pain in ass if want to plot log file's data
    /*fprintf(fp, header_format);*/
    rawStat(fp, p, p->op); fclose(fp);
#ifdef OS_UNIX
    printf("\033[7m");
#endif
    fprintf(stdout, header_format);
#ifdef OS_UNIX
    printf("\033[m");
#endif
    rawStat(stdout, p, p->op);
#undef header_format
#undef XEF_HEADER
}
