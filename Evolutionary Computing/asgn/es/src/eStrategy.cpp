#include "stdafx.hpp"
#include <assert.h>
#include <algorithm>
#include "eStrategy.hpp"

eStrategy::eStrategy(fitnessFunc eval) {
    this->eval = eval;
}

void eStrategy::initialize() {
    for(size_t indiv = 0; indiv < CONF_SIZE_POP; indiv++) {
        chromosome::data_block genes, sigmas;
        for(size_t gene = 0; gene < CONF_DIM; gene++) {
            sigmas.push_back(frand());
            genes.push_back(frand() * (CONF_BOUND_UPPER - CONF_BOUND_LOWER) + CONF_BOUND_LOWER);
        }
        auto c = chromosome(genes, sigmas);
        c.set_fitness(this->eval(c));
        this->parent_pop.push_back(c);
    }
    this->sort_population();
}

void eStrategy::gen_children() {
    this->child_pop.clear();
    assert(this->parent_pop.size() == CONF_SIZE_POP);
    // make a room for x1 elit selection in `eStrategy::survivor_selection()`
    for(size_t child = 0; child < CONF_SIZE_CHILD - 1; child++) {
        size_t
            pindex1 = get_rand(0, this->parent_pop.size()),
            pindex2 = get_rand(0, this->parent_pop.size());
        chromosome c =
            this->xover_and_mut(
                &this->parent_pop.at(pindex1),
                &this->parent_pop.at(pindex2));
        c.set_fitness(this->eval(c));
        this->child_pop.push_back(c);
    }
}

chromosome eStrategy::xover_and_mut(const chromosome* const p1, const chromosome* const p2) const {
    p1->validate(), p2->validate();
    assert(p1->sigmas.size() == p2->sigmas.size() && p1->sigmas.size() == CONF_DIM);
    double prob = frand();
    auto sp = p1;
    if(prob > 0.5) sp = p2;
    chromosome::data_block new_sigmas, new_genes;
    if(prob < CONF_PROB_MUT) {
        double step = exp(PARAM_TAW_PRIME * normrnd() + CONST_STEP_GLOB);
        auto xd = normrnd();
        for(size_t index = 0; index < sp->sigmas.size(); index++) {
            auto s = sp->sigmas.at(index) * step;
            if(s < CONF_RANGE_SIGMA && s > -CONF_RANGE_SIGMA)
                s = s / abs(s) * CONF_RANGE_SIGMA;
            new_sigmas.push_back(s);
            auto g = sp->genes.at(index) + s * xd;
            if(g > CONF_BOUND_UPPER) g = CONF_BOUND_UPPER;
            else if(g < CONF_BOUND_LOWER) g = CONF_BOUND_LOWER;
            new_genes.push_back(g);
        }
        assert(new_genes.size() == new_sigmas.size());
    } else if((1 - prob) < CONF_PROB_XOVER) {
        double rp = frand(), rq = frand();
        for(size_t index = 0; index < sp->genes.size(); index++) {
            auto g = (p1->genes.at(index) * rp + p2->genes.at(index) * rq);
            if(g > CONF_BOUND_UPPER) g = CONF_BOUND_UPPER;
            else if(g < CONF_BOUND_LOWER) g = CONF_BOUND_LOWER;
            new_genes.push_back(g);
            if(frand() < 0.5) new_sigmas.push_back(p1->sigmas.at(index));
            else new_sigmas.push_back(p2->sigmas.at(index));
        }
    }
    auto c = chromosome(new_genes, new_sigmas);
    if(this->eval(c) < sp->get_fitness())
        return c;
    return chromosome(*sp);
}

void eStrategy::sort_population() {
    auto cmp =
        [](const chromosome& a, const chromosome& b)
            { return a.get_fitness() < b.get_fitness(); };
    std::sort(this->parent_pop.begin(), this->parent_pop.end(), cmp);
    std::sort(this->child_pop.begin(),  this->child_pop.end(),  cmp);
}

void eStrategy::survivor_selection() {
    this->sort_population();
    // do x1 elit selection
    this->child_pop.push_back(this->parent_pop.front());
    assert(this->parent_pop.size() == CONF_SIZE_POP);
    assert(this->child_pop.size()  == CONF_SIZE_CHILD);
    population new_pop;
    if(CONF_GEN_REP_TYPE == FLAG_MU_AND_LAMBDA) {
        new_pop = this->child_pop;
    } else if(CONF_GEN_REP_TYPE == FLAG_MU_PLUS_LAMBDA) {
        std::merge(
            this->parent_pop.begin(), this->parent_pop.end(),
            this->parent_pop.begin(), this->parent_pop.end(),
            std::back_inserter(new_pop),
            [](const chromosome& a, const chromosome& b)
                { return a.get_fitness() < b.get_fitness(); });

    } else throw std::runtime_error("Unknown ops.");
    this->parent_pop.clear();
    std::copy(
        new_pop.begin(),
        new_pop.begin() + CONF_SIZE_POP,
        std::back_inserter(this->parent_pop));
    assert(this->parent_pop.size() == CONF_SIZE_POP);
}
