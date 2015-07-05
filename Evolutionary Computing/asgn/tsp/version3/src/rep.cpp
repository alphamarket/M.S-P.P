#include "stdafx.hpp"
#include "type.hpp"

void rawStat(FILE*, POPULATION*, IPTR);

void report(int __unused gen, POPULATION *p, IPTR pop)
{
	/* report generations stats */
    FILE *fp;
    if((fp = fopen(p->ofile, "a")) == NULL)
    { printf("error in opening file %s\n", p->ofile); exit(1); }
	else
    {
        rawStat(fp, p, pop);
        for(size_t i = 0; i < (size_t)p->popSize; i++) {
            for(size_t j = 0; j < (size_t)p->lchrom;j++)
                fprintf(fp, "%i ", pop[i].chrom[j]);
            fprintf(fp, " : %f\n", pop[i].fitness);
        }
        fclose(fp);
    }
    rawStat(stdout, p, pop);
}


void rawStat(FILE *fp, POPULATION *p, IPTR __unused pop)
{
    fprintf(fp,"%12d%-10s%10.2lf%-10s%10.2lf%-10s%10.2lf ",
        p->gen, " ", p->max, " ", p->avg, " ", p->min);
    fprintf(fp,"%12d%-10s%10.5lf%-1s%12d\n",
        p->highestEverGen, " ", p->highestEverFitness, " ", p->highestEverIndex);
}
