#include <math.h>
#include <stdio.h>
#include "deftypes.hpp"
#ifndef XFILE_COMPATIBLE
double distance(const cities* const, const size_t, const size_t);
double eval(
    const POPULATION* const p,
    const cities* const c,
    IPTR const i)
{
    double cost = 0.0;
    for(int gene = 0; gene <  p->lchrom - 1; gene++)
        cost += distance(c, i->chrom[gene], i->chrom[gene + 1]);
    return cost + distance(c, i->chrom[0], i->chrom[p->lchrom - 1]);
}
double distance(
    const cities* const c,
    const size_t index1,
    const size_t index2)
{
#ifdef CHECKSUM_IN_EFFECT
    assert(index1 < c->size() && index2 < c->size());
#endif
    const city  * const c1 = &c->at(index1),
                * const c2 = &c->at(index2);
    double xd = c1->x - c2->x, yd = c1->y - c2->y;
    return sqrt(xd * xd + yd * yd);
}
#else
#include <assert.h>
double distance(const POPULATION* const, const int** const, const size_t, const size_t) __nonnull();
double eval(const POPULATION* const p, const int** const dist, IPTR const i)  {
    double cost = 0.0;
    for(int gene = 0; gene <  p->lchrom - 1; gene++) {
        auto d = distance(p, dist, i->chrom[gene], i->chrom[gene + 1]);
        if(d < 0) return 0;
        cost += d;
    }
    auto d = distance(p, dist, i->chrom[p->lchrom - 1], i->chrom[0]);
    if(d < 0) return 0;
    return 1 / (cost + d);
}

double distance(
    const POPULATION* const p,
    const int** const dist,
    const size_t index1,
    const size_t index2)
{
    assert(index1 < (size_t)p->lchrom && index2 < (size_t)p->lchrom);
    return dist[index1][index2];
}
#endif
