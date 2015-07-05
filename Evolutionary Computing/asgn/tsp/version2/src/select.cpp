#include <cmath>
#include <stdexcept>
#include "stdafx.hpp"
#include "deftypes.hpp"
#include "select.helper.hpp"

void select(SELECT_TYPE st, const POPULATION* const p, size_t* selected1, size_t* selected2, bool (*cmp)(const double, const double)) {
    switch (st) {
    case SELECT_TYPE::ROULETTE:
        *selected1  = select_roulette(p->op, p->sumFitness, p->popSize),
        *selected2  = select_roulette(p->op, p->sumFitness, p->popSize);
    break;
    case SELECT_TYPE::TOP_TWO:
    {
#ifndef XFILE_COMPATIBLE
        double portion = 1 - double(4 * p->gen)/(3 * p->maxGen);
        select_two_top(p->op, p->popSize, (max(min(20.0, (double)p->popSize), portion * p->popSize)), selected1, selected2, cmp);
#else
        select_two_top(p->op, p->popSize, sqrt(p->popSize), selected1, selected2, cmp);
#endif
    }
    break;
    case SELECT_TYPE::TOURNAMENT:
    {
#ifndef XFILE_COMPATIBLE
        double portion = 1 - double(4 * p->gen)/(3 * p->maxGen);
        select_two_top(p->op, p->popSize, (max(min(20.0, (double)p->popSize), portion * p->popSize)), selected1, selected2, cmp);
#else
        size_t junk_pointer;
        select_two_top(p->op, p->popSize, sqrt(p->popSize), selected1, &junk_pointer, cmp);
        select_two_top(p->op, p->popSize, sqrt(p->popSize), &junk_pointer, selected2, cmp);
#endif
    }
    break;
    case SELECT_TYPE::RANDOM:
        select_random(p->op, p->popSize, selected1, selected2);
    break;
    default:
        throw std::invalid_argument("invalid input for selection op.");
    }
}
