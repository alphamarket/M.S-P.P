#ifndef POPULATIONTESTCASE_HPP
#define POPULATIONTESTCASE_HPP
#include "../hpp/teststrap.hpp"
#include "../../inc/gaConfig.hpp"
#include "../../inc/chromosome.hpp"
#include "../../inc/population.hpp"

using namespace CPP_TESTER;
using TSPGA::chromosome;
using TSPGA::population;
using TSPGA::gaConfig;
typedef population population_t;

namespace CPP_TESTER { namespace TESTS {
    /**
     * The test class, which inherit from `CPP_TESTER::testCase`
     */
    class populationTestCase : public testCase {
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
            population_t* p = NULL;
            gaConfig* conf  = NULL;
            BESURE(this->population_gen_check(p, conf));
            idelete(p);
            IS_NULL(p);
        }
        static TSPGA::distance fitnessFunc(const TSPGA::solution * const c) {
            return 10.F;
        }
    private:

        bool population_gen_check(population_t*& p, gaConfig*& conf) {
            IS_NULL(p); IS_NULL(conf);
            conf = new gaConfig(10, 0.9F, 0.1F);
            conf->setGenes_CountLimit(0, 10);
            p = new population_t(conf);
            p->setFitnessFunction(populationTestCase::fitnessFunc);
            SHOULD_BE(p->getFitnessFunction(), populationTestCase::fitnessFunc);
            NOT_NULL(p);
            SHOULD_BE(p->getCurrentGeneration_NO(), 0);
            IS_ZERO(p->getCurrentPopulationSize());
            SHOULD_THROW(p->init());
            conf->option("combinations.crossover.method", "XOI");
            conf->option("combinations.mutation.method",  "M_SWAP_RANGE");
            SHOULD_NOT_THROW(p->init());
            SHOULD_BE(p->getCurrentGeneration_NO(), 1);
            SHOULD_BE(p->getCurrentPopulationSize(), conf->getPopulation_Size());
            for(auto it = 0; it < p->getCurrentPopulation()->size(); it++) {
                auto s = p->getCurrentPopulation()->at(it).get();
                SHOULD_BE(populationTestCase::fitnessFunc(s), s->getFitness());
            }
            return true;
        }
    };
} }

#endif // POPULATIONTESTCASE_HPP
