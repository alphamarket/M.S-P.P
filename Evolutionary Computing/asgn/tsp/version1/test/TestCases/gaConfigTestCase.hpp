#ifndef GACONFIGTESTCASE_HPP
#define GACONFIGTESTCASE_HPP
#include "../hpp/teststrap.hpp"
#include "../../inc/gaConfig.hpp"
#include <vector>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace CPP_TESTER;
typedef TSPGA::gaConfig gc;

namespace CPP_TESTER { namespace TESTS {
    /**
     * The test class, which inherit from `CPP_TESTER::testCase`
     */
    class gaConfigTestCase : public testCase {
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
            gc* c = NULL;
            BESURE(this->test_allocation(c));
            BESURE(this->test_access(c));
            BESURE(this->test_option(c));
            NOT_NULL(c);
            idelete(c);
            IS_NULL(c);
            return true;
        }
    private:
        /**
         * @brief generates random string
         * @param length The length of string
         * @return The random string
         */
        std::string getRandString(const size_t length = 5) {
            std::string chars(
                    "abcdefghijklmnopqrstuvwxyz"
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    "1234567890"
                    "!@#$%^&*()"
                    "`~-_=+[{]}\\|;:'\",<.>/? ");
            boost::random::random_device rng;
            boost::random::uniform_int_distribution<> index_dist(0, chars.size() - 1);
            std::stringstream ss;
            for(int i = 0; i < length; ++i) { ss<< chars[index_dist(rng)]; }
            return ss.str();
        }
        /**
         * @brief test_getRandString
         * @return The test results
         */
        bool test_getRandString() {
            for(int i = 0; i < 100; i++) {
                SHOULD_BE(this->getRandString(i).length(), i);
            }
        }
        /**
         * @brief test_allocation
         * @param The gaConfig instance
         * @return The test results
         */
        bool test_allocation(gc*& c) {
            IS_NULL(c);
            c = new gc(100, 0.8, 0.2);
            SHOULD_BE(c->getPopulation_Size(), 100);
            SHOULD_BE(c->getCrossover_prob(), 0.8);
            SHOULD_BE(c->getMutation_prob(), 0.2);
            SHOULD_THROW(new gc(100, 1 + 0.8, 0.2));
            SHOULD_THROW(new gc(100, 0.8, 1 + 0.2));
            return true;
        }
        /**
         * @brief test_access
         * @param The gaConfig instance
         * @return The test results
         */
        bool test_access(gc * const c) {
            NOT_NULL(c);
            /* invalid input testing */
            SHOULD_THROW(c->setCrossover_prob(1.1));
            SHOULD_THROW(c->setCrossover_prob(-.1));
            SHOULD_THROW(c->setMutation_prob(1.1));
            SHOULD_THROW(c->setMutation_prob(-.1));
            SHOULD_THROW(c->setElitism_ratio(1.1));
            SHOULD_THROW(c->setElitism_ratio(-.1));
            SHOULD_THROW(c->setGenes_CountLimit(120, 10));
            SHOULD_THROW(c->setGenes_ValueLimit(20, 2.1));
            /* normal set/get testing*/
            c->setCrossover_prob(.23234);
            c->setElitism_ratio(.123456);
            c->setGenerationMaxCount(12000);
            c->setGenes_CountLimit(10, 120);
            c->setGenes_ValueLimit(2.1, 20);
            c->setMutation_prob(0.2);
            c->setPopulation_size(10);
            SHOULD_BE(c->getCrossover_prob(), .23234);
            SHOULD_BE(c->getElitism_ratio(), .123456);
            SHOULD_BE(c->getGenerationMaxCount(), 12000);
            SHOULD_BE(c->getGenes_Count_MAX_Limit(), 120);
            SHOULD_BE(c->getGenes_Count_MIN_Limit(), 10);
            SHOULD_BE(c->getGenes_Value_MAX_Limit(), 20);
            SHOULD_BE(c->getGenes_Value_MIN_Limit(), 2.1);
            SHOULD_BE(c->getMutation_prob(), 0.2);
            SHOULD_BE(c->getPopulation_Size(), 10);
            return true;
        }
        /**
         * @brief test_option
         * @param The gaConfig instance
         * @return The test results
         */
        bool test_option(gc * const c) {
            // test random string generator
            this->test_getRandString();
            IS_ZERO(c->options_count());
            // test that nothing has been added yet
            for(int i = 0; i < 1000; i++) NOT_TRUE(c->option_exists(this->getRandString()));
            typedef std::pair<std::string, std::string> pair;
            // our test key/value sample container
            std::vector<pair> kva;
            // gen. 100 sample options, with signature that `key.length < value.length`
            for(int i = 0; i < 100; i++) { std::string key; do { key = this->getRandString(5); } while(c->option_exists(key)); kva.push_back({ key, this->getRandString(10) }); }
            // add sample options
            BOOST_FOREACH(pair kv, kva) c->option(kv.first, kv.second);
            // report back the size of imported options
            SHOULD_BE(c->options_count(), kva.size());
            // test every options
            BOOST_FOREACH(pair kv, kva) {
                IS_TRUE(c->option_exists(kv.first));
                SHOULD_BE(c->option(kv.first), kv.second);
                // since all keys have same size, so by adding a fix string, it should not exist
                IS_FALSE(c->option_exists(kv.first + "CORRUPTED"));
            }
            // initial functionality check of getOption()
            SHOULD_BE(c->getOptions().size(), c->options_count());
            // get origin options container
            auto m = c->getOptions();
            // alter the retrieved options
            m.insert({"DUMY", "INSERT"});
            // the alter should not effect the config's option container
            SHOULD_BE(c->options_count() + 1, m.size());
            return true;
        }
    };
} }

#endif // GACONFIGTESTCASE_HPP
