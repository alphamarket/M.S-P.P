#include "stdafx.hpp"
#include "type.hpp"

double eval(POPULATION *p, IPTR pj)
{
    /* Called from gen.c and init.c */
    double cost = 0;
    for(int i = 0; i < p->lchrom; i++) {
		int j = (i + 1) % p->lchrom;
        if(p->city_costs[pj->chrom[i]][pj->chrom[j]] < 0) return 0;
        cost += p->city_costs[pj->chrom[i]][pj->chrom[j]];
    }
    return 1 / cost;
}
