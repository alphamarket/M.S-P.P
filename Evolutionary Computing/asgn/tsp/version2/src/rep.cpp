#include <stdio.h>
#include <stdlib.h>
#include "deftypes.hpp"
void raw_stat(FILE *fp, const POPULATION* const, IPTR);
void report(
#ifdef XFILE_COMPATIBLE
    const char* ofile,
#endif
    __unused int gen, POPULATION *p, IPTR pop)
{
    FILE *fp;
    if( (fp = fopen(
#ifdef XFILE_COMPATIBLE
                ofile
#else
                p->ofile
#endif
                , "a")) == NULL){
        printf("error in opening file %s\n", p->ofile);
        exit(1);
    }else{
        raw_stat(fp, p, pop);
#ifdef XFILE_COMPATIBLE
        for(int chromIndex = 0; chromIndex < p->popSize; chromIndex++) {
            auto solution = pop[chromIndex];
            for(int i = 0; i < p->lchrom; i++)
                fprintf(fp, "%i ", solution.chrom[i]);
            fprintf(fp, " : %f\n", solution.fitness);
        }
#endif
        fclose(fp);
    }
    raw_stat(stdout, p, pop);
}
void raw_stat(FILE *fp, const POPULATION* const p, IPTR)
{
    fprintf(fp,"%12d%-10s%10.2lf%-10s%10.2lf%-10s%10.2lf ",
        p->gen, " ", p->max, " ", p->avg, " ", p->min);
#ifndef XFILE_COMPATIBLE
    fprintf(fp,"%12d%-10s%10.2lf%-1s%12d\n",
        p->largestEverGen, " ", p->largestEverFitness, " ", p->largestEverIndex);
#else
    fprintf(fp,"%12d%-10s%10.5lf%-1s%12d\n",
        p->largestEverGen, " ", p->largestEverFitness, " ", p->largestEverIndex);
#endif
}
