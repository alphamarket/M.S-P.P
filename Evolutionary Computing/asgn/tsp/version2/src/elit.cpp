#ifndef ELIT_CPP
#define ELIT_CPP
#include <math.h>
#include <cstring>
#include <algorithm>
#include "stdafx.hpp"
#include "deftypes.hpp"

vector<einfo> select_elit(
    const POPULATION* const p,
    IPTR const op,
    IPTR (*allocator)(size_t),
    bool (*cmp)(const IPTR, const IPTR))
{
    vector<einfo> e;
    size_t es = floor(p->pElit * p->popSize);
    if(es == 0) es++;
    const size_t e_size =  es % 2 == 0 ? es : es + 1;
    if(e_size > (size_t)p->popSize) return vector<einfo>();
    for(size_t i = 0; i < e_size; i++) e.push_back(einfo(allocator(i), -1, true));
    IPTR best = e.front().i, worst = e.back().i;
    // pops from back of elit list
    auto pop_back = [&]() {
        // get the back
        einfo b = e.back();
        // dispose the end item
        e.pop_back();
        // if created by this function?
        if(b.is_fake) delete b.i;
        // update the best/worst
        best = e.front().i, worst = e.back().i;
    };
    for(size_t index = 0; index < (size_t)p->popSize; index++){
        IPTR s = &op[index];
        // don't count best/worst themselves, an also the solutions worst that the worst!!
        if(s == best || s == worst || !cmp(s, worst)) continue;
        // better than current best?
        if(cmp(s, best) || s->fitness == best->fitness) {
            // insert into front
            e.insert(e.begin(), einfo(s, index));
            // maintain size constraint
            pop_back();
            goto _continue;
        }
        else if(s->fitness == worst->fitness) {
            // [POLICY] tend to keep young solutions
            pop_back();
            // insert the young
            e.push_back(einfo(s, index));
            goto _continue;
        } else if(cmp(s, worst)) {
            // means `worst < s->fitness < best`
            // the worst has to be dismissed anyway, inorder to maintain the size constraint
            pop_back();
            // all we need to insert the `s` into a place that
            // our `best <= worst` constraint hold valid :)
            // performing a binary insertion
            e.insert(
                std::lower_bound(e.begin(), e.end(), einfo(s, index),
                    [&cmp](const einfo& i, const einfo& j)
                        { return cmp(i.i, j.i); }),
                einfo(s, index));
            goto _continue;
        }
        continue;
    _continue:
        continue;
    }
#ifdef CHECKSUM_IN_EFFECT
    assert(e.front().i->fitness <= e.back().i->fitness);
#endif
    return e;
}

#endif // ELIT_CPP
