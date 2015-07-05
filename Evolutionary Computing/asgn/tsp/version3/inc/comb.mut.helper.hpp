#ifndef COMB_MUT_HELPER_HPP
#define COMB_MUT_HELPER_HPP
#include <vector>
#include <algorithm>
#include "stdafx.hpp"
#include "type.hpp"
using namespace std;

__deprecated(void mut_swap_points_bounded(const POPULATION* const, IPTR const));
__deprecated(void mut_reverse_range_bounded(const POPULATION* const, IPTR const));
__deprecated(void mut_swap_reverse_range_bounded(const POPULATION* const, IPTR const));

vector<size_t> get_points(const POPULATION* const p) {
    vector<size_t> range;
    do {
        size_t rand1 = get_rand(0, p->lchrom), rand2 = get_rand(0, p->lchrom);
        range = { rand1, rand2 };
        std::sort(range.begin(), range.end());
    } while(range[0] == range[1]);
    return range;
}

void mut_swap_reverse(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    for(size_t i = range[0], j = range[1]; i < j; i++, j--)
        std::swap(pi->chrom[i], pi->chrom[j]);
}

void mut_swap_range(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    size_t block = range[1] - range[0];
    for(size_t i = 0; i < block; i++)
        std::swap(pi->chrom[(range[0] + i) % p->lchrom], pi->chrom[(range[1] + i) % p->lchrom]);
}

void mut_swap_points(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    std::swap(pi->chrom[range[0]], pi->chrom[range[1]]);
}

void mut_reverse_range(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    std::reverse(&pi->chrom[range[0]], &pi->chrom[range[1]]);
}

void mut_swap_points_bounded(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    std::swap(pi->chrom[range[0]], pi->chrom[range[1]]);
}

void mut_reverse_range_bounded(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    std::reverse(&pi->chrom[range[0]], &pi->chrom[range[1]]);
}

void mut_swap_reverse_range_bounded(const POPULATION* const p, IPTR const pi) {
    auto range = get_points(p);
    size_t block = range[1] - range[0];
    for(size_t i = 0; i < block; i++)
        std::swap(pi->chrom[range[0] + i], pi->chrom[range[1] + i]);
    for(size_t k = 0; k < 2; k++)
        for(size_t i = range[k], j = range[k] + block; i < j; i++, j--)
            std::swap(pi->chrom[i], pi->chrom[j]);
}

#endif // COMB_MUT_HELPER_HPP
