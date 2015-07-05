#ifndef COMB_HELPER_HPP
#define COMB_HELPER_HPP
#include <vector>
#include "stdafx.hpp"
#include "deftypes.hpp"
using namespace std;

vector<size_t> get_cross_points(const size_t no, const size_t lower_bound, const size_t upper_bound) {
    vector<size_t> p;
    do {
        p.clear();
        size_t size = no;
        while(size--) p.push_back((size_t)get_rand(lower_bound, upper_bound));
        // sort them to get in ASC manner
        std::sort(p.begin(), p.end());
        // make sure there is going to a crossover actually happen
    } while(no != 1 && (size_t)(p.back() - p.front()) < no);
#ifdef CHECKSUM_IN_EFFECT
    assert(p.size() == no);
#endif
    return p;
}

#include "comb.mut.helper.hpp"
#include "comb.xover.helper.hpp"

#endif // COMB_HELPER_HPP
