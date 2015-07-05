#include "stdafx.hpp"
#include "type.hpp"
#include <stdexcept>
#include "select.helper.hpp"

void select(SELECT_TYPE st, const POPULATION* const p, size_t* selected1, size_t* selected2) {
    switch (st) {
    case SELECT_TYPE::ROULETTE:
        *selected1  = select_roulette(p->op, p->sumFitness, p->popSize),
        *selected2  = select_roulette(p->op, p->sumFitness, p->popSize);
    break;
    case SELECT_TYPE::TOP_TWO:
        select_two_top(p->op, p->popSize, sqrt(p->popSize), selected1, selected2);
    break;
    case SELECT_TYPE::TOURNAMENT:
    {
        size_t junk_index;
        select_two_top(p->op, p->popSize, sqrt(p->popSize), selected1, &junk_index);
        select_two_top(p->op, p->popSize, sqrt(p->popSize), selected2, &junk_index);
    }
    break;
    case SELECT_TYPE::RANDOM:
        select_random(p->op, p->popSize, selected1, selected2);
    break;
//    case SELECT_TYPE::RANDOM_TOURNAMENT:
//        select_random_tournament(p->op, p->popSize, sqrt(p->popSize), selected1, selected2);
//    break;
    default:
        throw std::invalid_argument("invalid input for selection op.");
    }
}
