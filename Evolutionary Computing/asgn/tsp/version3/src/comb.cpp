//#define CHECKSUM_IN_EFFECT
#ifdef CHECKSUM_IN_EFFECT
#   include <set>
#   include <assert.h>
#endif
#include <vector>
#include <string.h>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include "type.hpp"
#include "comb.helper.hpp"
using namespace std;

void crossover(XOVER_TYPE xt, const POPULATION* const p, const IPTR p1, const IPTR p2, IPTR c1, IPTR c2, bool force = false)
{
    IPTR parents[2] = {p1, p2}, children[2] = {c1, c2};
    if(!force && frand() <= (1 - p->pCross)) {
        for(int index = 0; index < 2; index++)
        {
            children[index]->chrom = (decltype(c1->chrom))malloc(p->lchrom * sizeof(int));
            memcpy(children[index]->chrom, parents[index]->chrom, p->lchrom * sizeof(int));
            children[index]->fitness = parents[index]->fitness;
        }
        return;
    }
    switch(xt) {
    case XOVER_TYPE::TWO_POINT_PMX:
        xover_two_point(p, p1, p2, c1, c2);
    break;
    case XOVER_TYPE::ONE_POINT_PMX:
        xover_one_point(p, p1, p2, c1, c2);
    break;
    default:
        throw invalid_argument("invalid xover op. type requested!");
    }
#ifdef CHECKSUM_IN_EFFECT
    std::set<int> checksum1, checksum2;
    for(int i = 0; i < p->lchrom; i++) {
        auto x1 = c1->chrom[i], x2 = c2->chrom[i];
        checksum1.insert(x1), checksum2.insert(x2);
        assert(
                x1 < p->lchrom &&
                x2 < p->lchrom &&
                x1 * x2 >= 0
            );
    }
    assert(checksum1.size() == (size_t)p->lchrom && checksum2.size() == (size_t)p->lchrom);
#endif
}

void mutation(MUTATION_TYPE mt, const POPULATION* const p, IPTR const pi, bool force = false) {
    if(force && frand() < (1 - p->pMut)) return;
    switch(mt) {
    case MUTATION_TYPE::SWAP_RANGE:
        mut_swap_range(p, pi);
    break;
    case MUTATION_TYPE::REVERSE_RANGE:
        mut_reverse_range(p, pi);
    break;
    case MUTATION_TYPE::SWAP_POINTS:
        mut_swap_points(p, pi);
    break;
    case MUTATION_TYPE::SWAP_REVERSE:
        mut_swap_reverse(p, pi);
    break;
    default:
        throw invalid_argument("invalid mutation op. type requested!");
    }
#ifdef CHECKSUM_IN_EFFECT
    std::set<int> checksum;
    for(int i = 0; i < p->lchrom; i++) checksum.insert(pi->chrom[i]);
    assert(checksum.size() == (size_t)p->lchrom);
#endif
}
