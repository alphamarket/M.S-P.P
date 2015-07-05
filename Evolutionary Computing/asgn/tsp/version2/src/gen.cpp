#include <vector>
#include <stdio.h>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include "deftypes.hpp"
void select(SELECT_TYPE, const POPULATION* const, size_t*, size_t*, bool (*)(const double, const double)) __nonnull();
double eval(const POPULATION* const,
#ifndef XFILE_COMPATIBLE
    const cities* const
#else
    const int** const
#endif
    , IPTR const) __nonnull();
void mutation(MUTATION_TYPE, const POPULATION* const, IPTR const, bool = false) __nonnull();
void crossover(XOVER_TYPE, const POPULATION* const, IPTR, IPTR, IPTR, IPTR, bool = false) __nonnull();
std::vector<einfo> select_elit(
    const POPULATION* const,
    IPTR const,
    IPTR (*)(size_t),
    bool (*)(const IPTR, const IPTR)) __nonnull();

bool x_mut(const POPULATION* const p,
#ifndef XFILE_COMPATIBLE
    const cities* const c
#else
    const int** const c
#endif
    , IPTR pi, bool force = false) {
    if(!force && frand() < (1 - (p->pMut + p->pCross))) return false;
    size_t i1 = get_rand(0, p->lchrom), i2 = get_rand(0, p->lchrom);
    IPTR po = (IPTR)malloc(sizeof(INDIVIDUAL));
    po->chrom = (decltype(po->chrom))malloc(sizeof(int) * p->lchrom);
    memcpy(po->chrom, pi->chrom, sizeof(int) * p->lchrom);
    swap(po->chrom[i1], po->chrom[i2]);
    bool k = true;
__CHECK_FIT:
    po->fitness = eval(p, c, po);
    if(k && po->fitness < pi->fitness) {
        if(k)   swap(pi->chrom[i1], pi->chrom[i2]);
        else    memcpy(pi->chrom, po->chrom, sizeof(int) * p->lchrom);
        pi->fitness = po->fitness;
        k = true;
    } else if(k) {
        mutation(MUTATION_TYPE(get_rand(0, 4)), p, po, true);
        k = false;
        goto __CHECK_FIT;
    }
    free(po->chrom); free(po);
    return k;
}

void generation(POPULATION* const p,
#ifndef XFILE_COMPATIBLE
    const cities* const c
#else
    const int** const c
#endif
    , __unused const int genno)
{
    const auto select_op = SELECT_TYPE::TOURNAMENT;
    const auto xover_op  = XOVER_TYPE::TWO_POINT;
    const auto mut_op    = MUTATION_TYPE::SWAP_POINTS;
    const auto fitness_comparison =
        [](const IPTR a, const IPTR b)
#ifndef XFILE_COMPATIBLE
            { return a->fitness == numeric_limits<decltype(a->fitness)>::infinity() || a->fitness < b->fitness; };
#else
            { return  a->fitness > b->fitness; };
#endif

    size_t i = 0;
    vector<einfo> elit = select_elit(p, p->op,
        [](size_t) { return new INDIVIDUAL(); },
        fitness_comparison);

    for(auto& e : elit) {
        if(e.is_fake) continue;
        IPTR x = &p->np[i++];
        x->chrom = (int*)malloc(sizeof(int) * p->lchrom);
        memcpy(x->chrom, e.i->chrom, sizeof(int) * p->lchrom);
        x->fitness = e.i->fitness;
    }

    for( ;i < (size_t)p->popSize; i += 2){
        size_t ppi1 = -1, ppi2 = -1;
        select(select_op, p, &ppi1, &ppi2, [](const double a, const double b)
#ifndef XFILE_COMPATIBLE
            { return a < b; }
#else
            { return a > b; }
#endif
        );
        IPTR pc1 = &(p->np[i]), pc2 = &(p->np[i+1]);
        crossover(xover_op, p, &p->op[ppi1], &p->op[ppi2], pc1, pc2);
        mutation(mut_op, p, pc1);
        mutation(mut_op, p, pc2);
        pc1->fitness = eval(p, c, pc1);
        pc2->fitness = eval(p, c, pc2);
#if 0
        x_mut(p, c, pc1, (p->lm_best_fitness_repeat_count > 3));
        x_mut(p, c, pc2, (p->lm_best_fitness_repeat_count > 3));
#endif
    }
#if 0
    if(p->lm_best_fitness_repeat_count > 3) {
        auto pe = p->pElit;
        p->pElit *= 20;
        vector<einfo> x = select_elit(p, p->np,
                [](size_t) { return new INDIVIDUAL(); },
                fitness_comparison );
        p->pElit = pe;
        size_t t = sqrt(max(p->popSize, p->lchrom));
        for(size_t i = 0, j = 0; i < x.size(); i++, j = 0)
            if(!x[i].is_fake)
                while(j++ < t) x_mut(p, c, x[i].i, true);
    }
#endif
}
