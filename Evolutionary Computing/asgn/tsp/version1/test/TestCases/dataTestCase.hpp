#ifndef DATATESTCASE_HPP
#define DATATESTCASE_HPP
#include "../hpp/teststrap.hpp"
#include "../../inc/data.hpp"

using namespace CPP_TESTER;
using TSPGA::coordinate;
using TSPGA::coordinates;
typedef TSPGA::data data;

namespace CPP_TESTER { namespace TESTS {
    /**
     * The test class, which inherit from `CPP_TESTER::testCase`
     */
    class dataTestCase : public testCase {
    public:
        /**
         * Init your test class here.
         */
        bool __init() { return true; }
        /**
         * Free any resources used by your class here.
         */
        bool __dispose() { return true; }
        /**
         * Run your tests here
         */
        bool __run(int argc = 0, void** argv = NULL) {
            data* d = NULL;
            BESURE(this->test_allocation(d));
            BESURE(this->test_access(d));
            BESURE(this->test_measures(d));
            NOT_NULL(d);
            idelete(d);
            IS_NULL(d);
            return true;
        }
    private:
        bool test_allocation(data*& d) {
            IS_NULL(d);
            d = new data(new coordinates({ coordinate(1, 2), coordinate(3, 4) }));
            d->type("eUc_2d");
            NOT_NULL(d);
            SHOULD_BE(d->size(), 2);
            SHOULD_BE(d->type(), "EUC_2D");
            return true;
        }
        bool test_access(const data * const d) {
            coordinate x = (*d)[0];
            SHOULD_BE(x.longitude, 1);
            SHOULD_BE(x.latitude , 2);
            x = (*d)[1];
            SHOULD_BE(x.longitude, 3);
            SHOULD_BE(x.latitude , 4);
            x = d->get(0);
            SHOULD_BE(x.longitude, 1);
            SHOULD_BE(x.latitude , 2);
            x = d->get(1);
            SHOULD_BE(x.longitude, 3);
            SHOULD_BE(x.latitude , 4);
            SHOULD_THROW(d->get(100));
            /* should now change the original value  */
            (*d)[1] = coordinate(10, 20);
            x = (*d)[1];
            SHOULD_BE(x.longitude, 3);
            SHOULD_BE(x.latitude , 4);
            return true;
        }
        bool test_measures(const data* const d) {
            NOT_ZERO(d->getDistance(0, 1));
            IS_EQUAL(d->getDistance(0, 1), d->getDistance_EUC_2D(0, 1));
            SHOULD_THROW(d->getDistance(100, 200));
            return true;
        }
    };
} }

#endif // DATATESTCASE_HPP
