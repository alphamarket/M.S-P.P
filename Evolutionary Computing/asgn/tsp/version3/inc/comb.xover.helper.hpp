#ifndef COMB_XOVER_HELPER_HPP
#define COMB_XOVER_HELPER_HPP
#include <vector>
#include "stdafx.hpp"
#include "type.hpp"
using namespace std;

void xover_inner_two_point(const POPULATION* const p, const IPTR p1, const IPTR p2, IPTR c1, IPTR c2, vector<size_t> points) {
    IPTR parents[2] = {p1, p2}, children[2] = {c1, c2};
    for(int index = 0; index < 2; index++) {
        // fetch both parents
        auto cpr = parents[index], opr = parents[(index + 1) % 2];
        auto c = children[index];
        c->chrom = (decltype(c->chrom))malloc(p->lchrom * sizeof(int));
        memcpy(
            &c   ->chrom[points[0]],
            &cpr ->chrom[points[0]],
            (points[1] - points[0]) * sizeof(int));
        // create an instance of map(for performance reasons unordered map has been used)
        std::unordered_map<int, size_t>       indexes;
        // store the copied genes into the child in reverse value-index mod into map
        for(size_t i = points[0]; i < points[1]; i++) indexes.insert({c->chrom[i], i});
        // foreach gens in other parent
        for(size_t
            /* for parent */ i = 0,
            /* for  child */ j = 0; i < (size_t)p->lchrom; /* only parent*/ i++) {
            // if already copied into child?
            if(indexes.count(opr->chrom[i])) { continue; }
            // if the copy location in child has reached to previously copied ones?
            if(j >= points[0] && j < points[1]) /* jump to other side */ j = points[1];
            // copy the other parent's genes into the child
            c->chrom[j++] = opr->chrom[i];
        }
    }
}

void xover_two_point(const POPULATION* const p, const IPTR p1, const IPTR p2, IPTR c1, IPTR c2) {
    xover_inner_two_point(p, p1, p2, c1, c2, get_cross_points(2, 0, p->lchrom));
}

void xover_one_point(const POPULATION* const p, const IPTR p1, const IPTR p2, IPTR c1, IPTR c2) {
    vector<size_t> points;
    do {
        points = get_cross_points(1, 0, p->lchrom);
        points.push_back(p->lchrom - 1);
    } while(points[0] == points[1]);
    xover_inner_two_point(p, p1, p2, c1, c2, points);
}

#endif // COMB_XOVER_HELPER_HPP
