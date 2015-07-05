#include <stdio.h>
#include "deftypes.hpp"
void statistics(POPULATION *p, IPTR pop)
{
    p->sumFitness = pop[0].fitness;
    p->max = p->sumFitness;
    p->min = p->sumFitness;
    p->maxi = p->mini = 0;
    for(int i = 1; i < p->popSize; i++){
        IPTR pi = &(pop[i]);
        p->sumFitness += pi->fitness;
        if (p->max < pi->fitness) {
            p->max = pi->fitness;  p->maxi = i;
        }
        if (p->min > pi->fitness){
            p->min = pi->fitness;  p->mini = i;
        }
    }
    p->avg = p->sumFitness / (double) p->popSize;
    if(p->smallestEverFitness > p->min) {
        p->smallestEverFitness = p->min;
        p->smallestEverGen = p->gen;
        p->smallestEverIndex = p->mini;
    }
    if(p->largestEverFitness < p->max) {
        p->largestEverFitness = p->max;
        p->largestEverGen = p->gen;
        p->largestEverIndex = p->maxi;
    }
}
