#ifndef DEFTYPES_HPP
#define DEFTYPES_HPP

#include <limits>
#include <vector>
#include "stdafx.hpp"

typedef struct INDIVIDUAL {
    int *chrom;  /* the chromosome */
    double fitness;
    INDIVIDUAL()
        : INDIVIDUAL(
            nullptr,
#ifdef XFILE_COMPATIBLE
            1 /
#endif
            std::numeric_limits<double>::max()) { }
    INDIVIDUAL(int* chrom, double fitness)
        : chrom(chrom), fitness(fitness) { }
} INDIVIDUAL, *IPTR;
typedef struct city {
    double x;
    double y;
    city() { }
    city(const city& c) : x(c.x), y(c.y)
    { }
} city;
typedef struct {
    IPTR op;  /* arrays of individuals form an evolving population*/
    IPTR np;
    int lchrom; /* chromosome length */
    int gen; /* current generation */
    double sumFitness; /* statistics parameters for selection and tracking*/
    double max; /* progress */
    double avg;
    double min;
    double pCross; /* probability of Xover */
    double pMut;  /* probability of Mutation */
    double pElit; /* proportion of elit selection */
    double randomseed;
    double smallestEverFitness;
    int smallestEverGen;
    int smallestEverIndex;
    double largestEverFitness;
    int largestEverGen;
    int largestEverIndex;
    double lm_lastgen_best_fitness;
    size_t lm_best_fitness_repeat_count;
    int maxi; /* index of best individual in current population*/
    int mini; /* index of worst individual in current population*/
    int maxGen; /* when to stop */
    int popSize;/* population size */
    char *ofile; /* output File name */
    char *cityfile; /* output File name */
} POPULATION, *POPULATION_PTR;

typedef struct einfo{
    IPTR i;
    int index;
    bool is_fake;
    einfo(IPTR i, int index, bool is_fake = false)
        : i(i), index(index), is_fake(is_fake) { }
} einfo;


typedef std::vector<city> cities;

enum class SELECT_TYPE  {
    ROULETTE  = 0,
    TOP_TWO,
    TOURNAMENT,
    RANDOM,
    RANDOM_TOURNAMENT
};
enum class XOVER_TYPE   {
    TWO_POINT = 0,
    ONE_POINT
};
enum class MUTATION_TYPE{
    SWAP_RANGE= 0,
    SWAP_POINTS,
    SWAP_REVERSE,
    REVERSE_RANGE,
    SWAP_POINTS_BOUNDED,
    REVERSE_RANGE_BOUNDED,
    SWAP_REVERSE_RANGE_BOUNDED,
};

#endif
