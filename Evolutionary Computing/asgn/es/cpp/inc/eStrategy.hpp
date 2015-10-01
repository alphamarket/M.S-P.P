#ifndef ESTRATEGY_H
#define ESTRATEGY_H

#include "stdafx.hpp"
#include <vector>

#include "chromosome.hpp"

typedef std::vector<chromosome> population;
typedef double (*fitnessFunc)(const chromosome&);

class eStrategy
{
protected:
    fitnessFunc eval;
    population child_pop;
    population parent_pop;
public:
    eStrategy(fitnessFunc) __nonnull();
    void initialize();
    void gen_children();
    chromosome xover_and_mut(const chromosome* const, const chromosome* const) const __nonnull();
    void survivor_selection();
    inline const population& get_population() const { return this->parent_pop; }
protected:
    void sort_population();
};

#endif // ESTRATEGY_H
