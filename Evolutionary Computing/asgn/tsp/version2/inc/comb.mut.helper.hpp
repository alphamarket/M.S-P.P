#ifndef COMB_MUT_HELPER_HPP
#define COMB_MUT_HELPER_HPP
#include <vector>
#include <algorithm>
#include "stdafx.hpp"
#include "deftypes.hpp"
using namespace std;

__deprecated(void mut_swap_points_bounded(const POPULATION* const, IPTR const));
__deprecated(void mut_reverse_range_bounded(const POPULATION* const, IPTR const));
__deprecated(void mut_swap_reverse_range_bounded(const POPULATION* const, IPTR const));

vector<size_t> get_points(const POPULATION* const p, size_t* block, bool use_limited_range = true) {
    vector<size_t> range = {0, 0};
    if(use_limited_range)
        *block = floor(sqrt(p->lchrom));
    else
        *block = 0;
    while(!(range[0] + *block < range[1] && range[1] + *block < (size_t)p->lchrom))
        range = get_cross_points(2, 0, p->lchrom);
    return range;
}

void mut_swap_reverse(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block, false);
    for(size_t i = range[0], j = range[1]; i < j; i++, j--)
        std::swap(pi->chrom[i], pi->chrom[j]);
}

void mut_swap_range(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block);
    for(size_t i = 0; i < block; i++)
        std::swap(pi->chrom[range[0] + i], pi->chrom[range[1] + i]);
}

void mut_swap_points(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block, false);
    std::swap(pi->chrom[range[0]], pi->chrom[range[1]]);
}

void mut_reverse_range(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block, false);
    std::reverse(&pi->chrom[range[0]], &pi->chrom[range[1]]);
}

void mut_swap_points_bounded(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block);
    std::swap(pi->chrom[range[0]], pi->chrom[range[1]]);
}

void mut_reverse_range_bounded(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block);
    std::reverse(&pi->chrom[range[0]], &pi->chrom[range[1]]);
}

void mut_swap_reverse_range_bounded(const POPULATION* const p, IPTR const pi) {
    size_t block;
    auto range = get_points(p, &block);
    for(size_t i = 0; i < block; i++)
        std::swap(pi->chrom[range[0] + i], pi->chrom[range[1] + i]);
    for(size_t k = 0; k < 2; k++)
        for(size_t i = range[k], j = range[k] + block; i < j; i++, j--)
            std::swap(pi->chrom[i], pi->chrom[j]);
}

#endif // COMB_MUT_HELPER_HPP
