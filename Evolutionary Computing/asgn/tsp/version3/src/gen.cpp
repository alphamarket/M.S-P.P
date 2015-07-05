#include "stdafx.hpp"
#include <vector>
#include <string.h>
#include <algorithm>
#include "type.hpp"

double  eval(POPULATION*, IPTR) __nonnull();
void    select(SELECT_TYPE, const POPULATION* const, size_t*, size_t*) __nonnull();
void    mutation(MUTATION_TYPE, const POPULATION* const, IPTR const, bool = false) __nonnull();
void    crossover(XOVER_TYPE, const POPULATION* const, IPTR, IPTR, IPTR, IPTR, bool = false) __nonnull();

void generation(POPULATION *p, int __unused t)
{
    struct indiv_cmp {
        bool operator() (const INDIVIDUAL& a,const INDIVIDUAL& b)
            { return a.fitness > b.fitness; }
    };
    size_t p1, p2, i = 2;
    IPTR pi, piPlus1, om1, om2;
    std::sort(p->op, p->op + p->popSize, indiv_cmp());
    for(size_t elit = 0; elit < i; elit++)
    {
        memcpy(p->np[elit].chrom, p->op[elit].chrom, sizeof(int) * p->lchrom);
        p->np[elit].fitness = p->op[elit].fitness;
    }
    for(; i < (size_t)p->popSize; i += 2)
    {
        select(SELECT_TYPE::TOURNAMENT, p, &p1, &p2);
        om1 = &(p->op[p1]);
        om2 = &(p->op[p2]);
		pi      = &(p->np[i]);
        piPlus1 = &(p->np[i+1]);
        crossover(XOVER_TYPE::TWO_POINT_PMX, p, om1, om2, pi, piPlus1);
        mutation(MUTATION_TYPE::SWAP_POINTS, p, piPlus1);
        mutation(MUTATION_TYPE::SWAP_POINTS, p, pi);
        pi->fitness = eval(p, pi);
        piPlus1->fitness = eval(p, piPlus1);
    }
}
