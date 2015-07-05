#ifndef SELECT_HELPER_HPP
#define SELECT_HELPER_HPP
#include <limits.h>
#include "stdafx.hpp"
#include "deftypes.hpp"

__deprecated(void select_random_tournament(const IPTR,const size_t,const size_t,size_t*,size_t*));

size_t select_roulette(IPTR pop, double sumFitness, size_t popsize)
{
    /* select a single individual by roulette wheel selection */
    double r,partsum;
    size_t i;
    partsum = 0.0; i = 0;
    r = frand() * sumFitness;
    i = -1;
    do{
        i++;
        partsum += pop[i].fitness;
    } while (partsum < r && i < popsize - 1) ;
    return i;
}

void select_two_top(
    const IPTR pop,
    const size_t pop_size,
    const size_t election_size,
    size_t* output_index_1,
    size_t* output_index_2,
    bool (*cmp)(const double, const double))
{
    typedef double cost;
    size_t es = election_size;
    // create 2 main tournoment group
    // one with smallest cost value, the other is second small cost value
    size_t index = get_rand(0, pop_size);
    vector<pair<size_t, cost>>
        min1({{index,
#ifndef XFILE_COMPATIBLE
            std::numeric_limits<double>::max()
#else
            pop[index].fitness
#endif
         }});
    index = get_rand(0, pop_size);
    vector<pair<size_t, cost>>
        min2({{index,
#ifndef XFILE_COMPATIBLE
            std::numeric_limits<double>::max()
#else
            pop[index].fitness
#endif
        }});
    es -= 2;
    // while size of tournament
    while(es--) {
        // fetch an solution index
        size_t index = get_rand(0, pop_size);
        // get its cost value
        cost x = pop[index].fitness;
        // if it is the smallest one?
        if(cmp(x, min1.at(0).second))
            // empty the small container & add it
            { min1.clear(); min1.push_back({index, x}); }
        // if it is among other smallest costs?
        else if(x == min1.at(0).second)
            // add it
            min1.push_back({index, x});
        // if it does not belong to smallest group?
        // and is smaller than secondary small group?
        else if(cmp(x, min2.at(0).second))
            // purge them out, and place it into the second small group
            { min2.clear(); min2.push_back({index, x}); }
        // if it is among other secondary small costs?
        else if(x == min2.at(0).second)
            // add it
            min2.push_back({index, x});
    }
    /**
     * Purpose of this section is to combine two small groups
     * and properly fetch(randomly) 2 index as the tournament's winners
     */
    /**
     * Phase 1) if there is only an item in at least one of the above lists
     */
    if(min1.size() == 1 && min2.size() == 1)
    { *output_index_1 = min1.at(0).first, *output_index_2 = min2.at(0).first; return; }
    if(min1.size() == 1)
    { *output_index_1 = min1.at(0).first; *output_index_2 = min2.at(get_rand(0, min2.size())).first; return; }
    if(min2.size() == 1)
    { *output_index_2 = min2.at(0).first; *output_index_1 = min1.at(get_rand(0, min1.size())).first; return; }
    /**
     * Phase 2) if there are more than one item in both lists
     */
    const size_t t_size = min1.size() + min2.size() - 1;
    size_t p[2];
    // fetch two index
    do{ p[0] = get_rand(0, t_size), p[1] = get_rand(0, t_size); } while(p[0] == p[1]);
    // assign the two combined indexes as tournament's winners
    if(p[0] < min1.size()) *output_index_1 = min1.at(p[0]).first;
    else *output_index_1 = min2.at(p[0] - min1.size()).first;
    if(p[1] < min1.size()) *output_index_2 = min1.at(p[1]).first;
    else *output_index_2 = min2.at(p[1] - min1.size()).first;
    return;
}

void select_random(
    __unused const IPTR,
    size_t pop_size,
    size_t* output_index_1,
    size_t* output_index_2)
{
    *output_index_1 = get_rand(0, pop_size);
    *output_index_2 = get_rand(0, pop_size);
    if(frand() < 0.1) {
        *output_index_1 = get_rand(0, pop_size / 2);
        *output_index_2 = get_rand(pop_size / 2, pop_size / 2);
    }
}

void select_random_tournament(
    __unused const IPTR,
    const size_t pop_size,
    const size_t election_size,
    size_t* output_index_1,
    size_t* output_index_2)
{
    size_t es = election_size;
    auto e = new size_t[election_size];
    while(es--) e[es] = get_rand(0, pop_size);
    *output_index_1 = e[get_rand(0, election_size)];
    *output_index_2 = e[get_rand(0, election_size)];
}


#endif // SELECT_HELPER_HPP
