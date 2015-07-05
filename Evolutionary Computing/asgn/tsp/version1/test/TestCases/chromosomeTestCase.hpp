#ifndef CHROMOSOMETESTCASE_HPP
#define CHROMOSOMETESTCASE_HPP
#include <set>
#include <boost/foreach.hpp>
#include "../hpp/teststrap.hpp"
#include "../../inc/Timer.hpp"
#include "../../inc/chromosome.hpp"

using namespace CPP_TESTER;
using TSPGA::chromosome;
typedef chromosome<float> chromosome_t;

namespace CPP_TESTER { namespace TESTS {
    /**
     * The test class, which inherit from `CPP_TESTER::testCase`
     */
    class chromosomeTestCase : public testCase {
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
            chromosome_t* c = NULL;
            BESURE(this->alloc_check(c));
            BESURE(this->gene_set_check(c));
            BESURE(this->clone_check(c));
            BESURE(this->fitness_check(c));
            BESURE(this->crossover_check(c));
            BESURE(this->mutation_check(c));
            NOT_NULL(c);
            idelete(c);
            IS_NULL(c);
        }
    private:
        /**
         * @brief Generates random genes
         * @param size The size of genes
         */
        vector<float>* getGenes(size_t size) const {
            vector<float>* g = new vector<float>();
            updateseed();
            for(int i = 0; i < size; i++) g->push_back(i);
            std::random_shuffle(g->begin(), g->end());
            return g;
        }
        /**
         * @brief Checks if a chromosome is valid
         * @param c The chromosome
         * @return The test result
         */
        bool is_valid_chromosome(const chromosome_t* const c) const {
            std::set<float> s(c->_genes->begin(), c->_genes->end());
            return s.size() == c->_genes->size();
        }
        /**
         * @brief Checks if a list of chromosomes are valid
         * @param cc The chromosome collection
         * @return The test result
         */
        bool is_valid_chromosome(vector<chromosome_t*>* cc, bool auto_dispose = false) const {
            if(!cc) return true;
            bool flag = true;
            BOOST_FOREACH(chromosome_t* c, *cc) {
                if(!this->is_valid_chromosome(c))
                    flag = false;
                if(auto_dispose)
                    idelete(c);
            }
            if(auto_dispose)
                idelete(cc);
            return flag;
        }
        /**
         * @brief Tests is_valid_cross_over()
         * @return test result
         */
        bool is_valid_chromosome_check() {
            auto v = vector<float>({1, 2, 3, 4, 5});
            auto g1 = new chromosome_t(10), g2 = new chromosome_t();
            g2->genes(&v);
            auto vc = vector<vector<chromosome_t*>*>({ new vector<chromosome_t*>({ g1 }),
                       new vector<chromosome_t*>({ g2 }),
                       new vector<chromosome_t*>({ g1, g2 }),
                       new vector<chromosome_t*>({ g2, g2 }) });
            IS_FALSE(this->is_valid_chromosome(vc[0]));
            IS_TRUE(this->is_valid_chromosome(vc[1]));
            IS_FALSE(this->is_valid_chromosome(vc[2]));
            IS_TRUE(this->is_valid_chromosome(vc[3]));
            delete g1, g2;
            return true;
        }
        /**
         * @brief tests on allocation
         * @param c the chromosome to init, should be initially null
         * @return test result
         */
        bool alloc_check(chromosome_t*& c) {
            IS_NULL(c);
            c = new chromosome<float>();
            NOT_NULL(c);
            NOT_NULL(c->_genes);
            IS_ZERO(c->getSize());
            IS_EQUAL(c->_genes->size(), c->getSize());
            IS_EQUAL(c->getFitness(), float());
            NOT_ZERO(c->to_str().length());
            IS_NULL(c->fitFunc);
            return true;
        }
        /**
         * @brief gene_set_check
         * @param c the chromosome assign, should be initialized
         * @return
         */
        bool gene_set_check(chromosome_t* const c) {
            // chromosome should be initialized
            NOT_NULL(c);
            // create a sample gene vect
            auto g = vector<float>({1, 2, 3, 4, 5});
            // assign a initial genes
            c->genes(&g);
            // check allocation address
            NOT_SAME(&g, c->_genes);
            // validate the gene's size
            SHOULD_BE(c->getSize(), 5);
            /**
             * check gene's access ops.
             * indexes [0...+4] should be accessible
             * indexes [-5..-1] are equivalent to [0...+4] indexes
             * method `chromosome<>::gene()` should work as getter/setter for genes,
             *    i.e by assigning `chromosome<>::gene(index)` the `chromosome<>::genes[index]` value should change
             */
            for(int gene = 0; gene < c->getSize(); gene++) {
                // validate the genes' value
                SHOULD_BE(c->gene(gene), gene + 1);
                // validate the returned typeof `chromosome<>::gene()` for the index
                IS(c->gene(gene), float);
                // cout<<"IS_EQUAL(c->gene("<<gene<<"), c->gene("<<gene - int(c->getSize())<<"))"<<endl;
                // test negative indexing ops
                IS_EQUAL(c->gene(gene), c->gene(gene - int(c->getSize())));
                // test assigning feature of `chromosome<>::gene()`
                // just for concrete test using negative indexing for settings
                c->gene(gene - int(c->getSize())) = 10 + gene;
                // and positive indexing for getings
                SHOULD_BE(c->gene(gene), 10 + gene);
            }
            // final check for genes size
            SHOULD_BE(c->getSize(), 5);
            return true;
        }
        /**
         * @brief clone_check
         * @param c the chromosome assign, should be initialized
         * @return test result
         */
        bool clone_check(const chromosome_t* const) {
            auto gen_genes = [](size_t size) {
                vector<float>* g = new vector<float>();
                for(int i=0;i<size;i++)
                    g->push_back(getRand(INT_MIN, INT_MAX));
                return g;
            };
            vector<chromosome_t*> k1, k2;
            for(size_t i = 0; i < 1000; i++) {
                chromosome_t* c = new chromosome_t();
                c->genes(gen_genes(1000));
                k1.push_back(c);
                k2.push_back(c->clone());
                SHOULD_BE(k1.size(), i + 1);
                SHOULD_BE(k2.size(), i + 1);
                PSAME_POINT(k1.back(), c);
            }
            BOOST_FOREACH(chromosome_t* c, k1) {idelete(c); IS_NULL(c);}
            k1.clear();
            BOOST_FOREACH(chromosome_t* c, k2) {
                BESURE(c);
                BESURE(c->_genes);
                SHOULD_BE(c->_genes->size(), 1000);
                BOOST_FOREACH(float x, *c->_genes);
                idelete(c);
                IS_NULL(c);
            }
            k2.clear();
            return true;
        }
        /**
         * @brief fitness_check
         * @param c the chromosome assign, should be initialized
         * @return test result
         */
        bool fitness_check(chromosome_t* const c) {
            // chromosome should be initialized
            NOT_NULL(c);
            // fitness calc. func. should be initially null
            IS_NULL(c->fitFunc);
            // the init. fitness value should be float's default value
            SHOULD_BE(c->getFitness(), float());
            // define/set a `sum(genes)` linear fit. func.
            c->setFitnessFunc([](const chromosome<float>  * const c) {
                float f = 0;
                for(int index = 0; index < c->getSize(); index++)
                    f += (*c->_genes)[index];
                return f;
            });
            // fitness calc. func. should NOT be null NOW
            NOT_NULL(c->fitFunc);
            // considering our knowledge on genes setted in previous checks
            // the updated fitness value should not be the float's default ANYMORE
            SHOULD_NOT_BE(c->calcFitness(), float());
            IS_EQUAL(c->getFitness(), c->calcFitness());
            return true;
        }
        /**
         * @brief crossover_check
         * @param c the chromosome assign, should be initialized
         * @return test result
         */
        bool crossover_check(const chromosome_t* const c) {
            NOT_NULL(c);
            BESURE(this->is_valid_chromosome_check());
            auto xcheck = [](
                    const CPP_TESTER::TESTS::chromosomeTestCase* const $this,
                    const TSPGA::CROSSOVER_OPERATOR x)  {
                for(int testing = 0; testing < 100; testing++) {
                    chromosome_t *c1, *c2;
                    vector<chromosome_t*>* cc;
                    c1 = new chromosome_t();
                    c2 = new chromosome_t();
                    c1->setDefaultCrossoverOP(x);
                    c2->setDefaultCrossoverOP(x);
                    delete c1->_genes, c2->_genes;
                    c1->_genes = $this->getGenes(1000);
                    c2->_genes = $this->getGenes(c1->getSize());
                    IS_EQUAL(c1->getSize(), c2->getSize());
                    cc = *c1 + *c2;
                    NOT_ZERO(cc->size());
                    BESURE($this->is_valid_chromosome(cc));
                    idelete(c1); idelete(c2);
                    IS_NULL(c1); IS_NULL(c2);
                    // test after deleting parents, does the children remain intacted?
                    NOT_ZERO(cc->size());
                    BOOST_FOREACH(chromosome_t* c, *cc) {
                        NOT_ZERO(c->getSize());
                        NOT_NULL(c->_genes);
                        idelete(c);
                    }
                }
            };
            TSPGA::Timer_mil timer;
            // xcheck(this, TSPGA::CROSSOVER_OPERATOR::XO);
            cout<<"[SKIPPED] Crossover of 100x1000 chromosomes with XO method took: "<<timer.elapsed()<<"ms"<<endl;
            timer.reset();
            xcheck(this, TSPGA::CROSSOVER_OPERATOR::XOI);
            cout<<"Crossover of 100x1000 chromosomes with XOI method took: "<<timer.elapsed()<<"ms"<<endl;
            return true;
        }
        /**
         * @brief mutation_check
         * @param c the chromosome assign, should be initialized
         * @return test result
         */
        bool mutation_check(const chromosome_t* const c) {
            NOT_NULL(c);
            BESURE(this->is_valid_chromosome_check());
            auto xcheck = [](
                    const CPP_TESTER::TESTS::chromosomeTestCase* const $this,
                    const TSPGA::MUTATION_OPERATOR m)  {
                for(int testing = 0; testing < 100; testing++) {
                    chromosome_t *c1 = new chromosome_t();
                    delete c1->_genes;
                    c1->_genes = $this->getGenes(1000);
                    chromosome_t* c2 = new chromosome_t();
                    delete c2->_genes;
                    c2->_genes = new vector<float>(*c1->_genes);
                    c1->setDefaultMutationOP(m);
                    c2->setDefaultMutationOP(m);
                    BESURE($this->is_valid_chromosome(c1));
                    BESURE($this->is_valid_chromosome(c2));
                    PNOT_SAME_POINT(c1->_genes, c2->_genes);
                    IS_EQUAL(*c1->_genes, *c2->_genes);
                    int i = 100;
                    while(i--) {
                        NOT_EQUAL(*((*c1)++)._genes, *c2->_genes);
                        delete c2->_genes;
                        c2->_genes = new vector<float>(*c1->_genes);
                    }
                    BESURE($this->is_valid_chromosome(c1));
                    BESURE($this->is_valid_chromosome(c2));
                    delete c1;
                    delete c2;
                }
            };
            TSPGA::Timer_mil timer;
            xcheck(this, TSPGA::MUTATION_OPERATOR::M_SWAP);
            cout<<"Mutation of 100x1000 chromosomes with M_SWAP method took: "<<timer.elapsed()<<"ms"<<endl;
            timer.reset();
            xcheck(this, TSPGA::MUTATION_OPERATOR::M_SWAP_RANGE);
            cout<<"Mutation of 100x1000 chromosomes with M_SWAP_RANGE method took: "<<timer.elapsed()<<"ms"<<endl;
            return true;
        }
    };
} }
#endif
