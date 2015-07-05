#ifndef TYPES_HPP
#define TYPES_HPP
#include <vector>
typedef struct 
{
    int    *chrom;         /* the chromosome */
    double fitness;
    int    parent1;
    int    parent2;
}INDIVIDUAL;
typedef struct 
{
    double x;
    double y;
} city;
typedef INDIVIDUAL *IPTR;
typedef struct
{
    IPTR   op;          /* arrays of individuals form an evolving population*/
    IPTR   np;
    int    lchrom;      /* chromosome length */
    int    gen;         /* current generation */
    double sumFitness;  /* statistics parameters for selection and tracking*/
    double max;         /* progress */
    double avg;
    double min;
    double pCross;      /* probability of Xover */
    double pMut;        /* probability of Mutation */
    double highestEverFitness;
    int    highestEverGen;
    int    highestEverIndex;
    int    maxi;        /* index of best individual in current population*/
    int    mini;        /* index of worst individual in current population*/
    int    maxGen;      /* when to stop */
    int    popSize;     /* population size */
    char   *ofile;      /* output File name */
    std::vector<std::vector<int>> city_costs;
} POPULATION;

enum class SELECT_TYPE  {
    ROULETTE  = 0,
    TOP_TWO,
    TOURNAMENT,
    RANDOM,
    RANDOM_TOURNAMENT
};
enum class XOVER_TYPE   {
    TWO_POINT_PMX = 0,
    ONE_POINT_PMX
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
