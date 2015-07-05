#ifndef CHROMOSOME_HPP
#define CHROMOSOME_HPP
#include "stdafx.hpp"
#include <vector>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <boost/lexical_cast.hpp>
using namespace std;

namespace TSPGA {
    typedef enum {
        XO, XOI
    } CROSSOVER_OPERATOR;
    typedef enum {
        M_SWAP, M_SWAP_RANGE
    } MUTATION_OPERATOR;
    /**
     * The chromosome template class
     * template desc:
     *      G = Genes' Type
     *      F = Fitness Value's Type (default: float)
     */
    template<class G, class F = float>
    class chromosome
    {
        CROSSOVER_OPERATOR  default_crossover_op;
        MUTATION_OPERATOR   default_mutation_op;
    protected:
        /**
         * @brief The related fitness value
         */
        F fitness;
    public:
        /**
         * Typedef of a fitness function
         */
        typedef F (*fitnessFunc_t)(const chromosome<G, F> * const c);
        /**
         * @brief The fitness function's pointer
         */
        fitnessFunc_t fitFunc;
        /**
         * @brief The genes
         */
        vector<G>* _genes = NULL;
    public:
        /**
         * Definition/Implementation of relational operators of chromosome class
         */
        inline bool operator< (const chromosome<G, F>& c) const { return this->fitness < c.fitness; }
        inline bool operator> (const chromosome<G, F>& c) const { return this->fitness > c.fitness; }
        inline bool operator==(const chromosome<G, F>& c) const { return this->fitness == c.fitness; }
        /**
         * @brief operator << To stream operator
         */
        inline friend std::ostream& operator<< (std::ostream& stream, const chromosome<G, F>& c) { stream<<c.to_str(); return stream; }
        inline friend std::ostream& operator<< (std::ostream& stream, const chromosome<G, F>* c) { stream<<c->to_str(); return stream; }
        /**
         * @brief Default constructor
         */
        chromosome(CROSSOVER_OPERATOR default_crossover_op = XOI, MUTATION_OPERATOR default_mutation_op = M_SWAP)
            : default_crossover_op(default_crossover_op),
              default_mutation_op(default_mutation_op)
        { this->_genes = new vector<G>(); this->fitFunc = NULL; this->fitness = F(); }
        /**
         * @brief Vector constructor
         * @param The genes vector
         */
        chromosome(size_t genes_size, CROSSOVER_OPERATOR default_crossover_op = XOI, MUTATION_OPERATOR default_mutation_op = M_SWAP)
            : chromosome(default_crossover_op, default_mutation_op)
        { idelete(this->_genes); this->_genes = new vector<G>(genes_size); }
        /**
         * @brief Copy constructor
         */
        chromosome(const chromosome<G, F>& c)
            : chromosome(c.default_crossover_op, c.default_mutation_op)
        { this->genes(c._genes); this->fitness = c.fitness; this->fitFunc = c.fitFunc; }
        /**
         * @brief explicit version of copy ctor
         * @return The cloned chromosome of current instance
         */
        inline chromosome<G, F>* clone() const { return new chromosome<G, F>(*this); }
        /**
         * @brief Destructore
         */
        ~chromosome() { idelete(this->_genes); }
        /**
         * @brief gset Get/Set genes' value
         * @param index The gene's index
         * @return The gene
         */
        inline G    gene(const int index) const { if(index >=(int)this->getSize())throw overflow_error("Index#" + std::to_string(index) + " of " + std::to_string(this->getSize())+" overflow"); if(index >= 0) return this->_genes->at(index); return this->_genes->at(this->getSize() + index); }
        inline G&   gene(const int index)       { if(index >=(int)this->getSize())throw overflow_error("Index#" + std::to_string(index) + " of " + std::to_string(this->getSize())+" overflow"); if(index >= 0) return this->_genes->at(index); return this->_genes->at(this->getSize() + index); }
        /**
         * @brief Assigns genes to this chromosome's genes(clears previous genes)
         * @param genes The target geness
         */
        inline void genes(const vector<G>* const _genes) __nonnull() {
            if(this->_genes) idelete(this->_genes);
            this->_genes = new vector<G>(_genes->size());
            std::copy(_genes->begin(), _genes->end(), this->_genes->begin());
        }
        /**
         * @brief Sets fitness calculator function related to current chromosome instance
         */
        inline void setFitnessFunc(fitnessFunc_t ff) __nonnull() { this->fitFunc = ff; }
        /**
         * @brief Calculate the fitness value based on assigned fitness function
         * @return Returns the calculated fitness
         */
        inline F calcFitness() { if(this->fitFunc == NULL) throw std::invalid_argument("The fitness function is null"); this->setFitness(this->fitFunc(this)); return this->getFitness(); }
        /**
         * @brief Get size of chromosome
         * @return The size of chromosome
         */
        inline size_t getSize() const { return this->_genes->size(); }
        /**
         * @brief Set fitness value
         * @param The fitness value
         */
        inline void setFitness(F fitness) { this->fitness = fitness; }
        /**
         * @brief Get fitness value
         * @return Returns the fitness value
         */
        inline F getFitness() const { return this->fitness; }
        /**
         * @brief Set default crossover operators
         */
        inline void setDefaultCrossoverOP(CROSSOVER_OPERATOR x) { this->default_crossover_op = x; }
        /**
         * @brief Get default crossover operators
         */
        inline CROSSOVER_OPERATOR getDefaultCrossoverOP() const { return this->default_crossover_op; }
        /**
         * @brief Set default mutation operators
         */
        inline void setDefaultMutationOP(MUTATION_OPERATOR x) { this->default_mutation_op = x; }
        /**
         * @brief Get default mutation operators
         */
        inline MUTATION_OPERATOR getDefaultMutationOP() const { return this->default_mutation_op; }
        /**
         * @brief   Creates a string version of current choromosome
         * @param   gene_to_str(OPTIONAL) A function handler to convert genes to proper string values
         * @return  Returns a string with format of `[genes value]#{fitness}`
         */
        inline const string to_str(char* (*gene_to_str)(G) = NULL) const { stringstream s; s<<"[ "; for(int i = 0; i < this->getSize(); i++) { if(gene_to_str == NULL) s<<(*this->_genes)[i]; else s<<gene_to_str((*this->_genes)[i]); s<<(i + 1 != this->getSize() ? ", ": ""); } s<<" ]#{"<<this->fitness<<"}"; return s.str(); }
        /**
         * @brief The default crossover operator
         * @param The second chromosome
         * @return The vector of children chromosomes
         */
        inline vector<chromosome<G, F>*>* operator+ (const chromosome<G, F>& c) const { return this->crossover(&c, this->getDefaultCrossoverOP()); }
        /**
         * @brief The default mutation operator
         */
        inline chromosome<G, F>& operator++ (int/* postfix */)  { this->mutation(this->getDefaultMutationOP()); return *this; }
        inline chromosome<G, F>& operator++ (/* prefix */)      { return (*this)++; }
        /**
         * @brief The mutation operator
         * @param m The desired mutation operator
         */
        void mutation(const MUTATION_OPERATOR m) {
            // for genes' size less than 2, mutation is not allowed
            if(this->getSize() < 2) return;
            switch(m) {
                case MUTATION_OPERATOR::M_SWAP:
                {
                    size_t p1 = 0, p2 = 0;
                    while(p1 == p2) { p1 = getRand(0, this->getSize()); p2 = getRand(0, this->getSize()); }
                    // swap two genes in current chromosome
                    iter_swap(this->_genes->begin() + p1, this->_genes->begin() + p2);
                    break;
                }
                case MUTATION_OPERATOR::M_SWAP_RANGE:
                {
                    // cache data for large scale mutation process
                    static double size_sqrt = -1; static size_t size = -1;
                    // validate the cache data an update them if necessary
                    if(size_sqrt < 0 || size != this->getSize()) { size = this->getSize(); size_sqrt = sqrt(this->getSize()); }
                    // alloc and init
                    size_t p1 = 0, p2 = 0, range = size_t(size_sqrt);
                    // if the size of current chromosome is less or eqaul than its sqrt value?
                    // M_SWAP_RANGE will perform bad, execute M_SWAP mutation procedure instead
                    if(ceil(size_sqrt) >= this->getSize()) { this->mutation(MUTATION_OPERATOR::M_SWAP); return; }
                    // find some proper points
                    while(abs(long(p1) - long(p2)) < range) { p1 = getRand(0, this->getSize() - range); p2 = getRand(0, this->getSize() - range); }
                    // swap two ranges in current chromosome
                    swap_ranges(this->_genes->begin() + p1, this->_genes->begin() + p1 + range, this->_genes->begin() + p2);
                    break;
                }
                default:
                    throw std::out_of_range("Undefined mutation operation.");
            }
        }
        /**
         * @brief   The crossover operator
         * @param   c The second chromosome
         * @param   x The desired crossover operator
         * @return  The vector of children chromosomes
         */
        vector<chromosome<G, F>*>* crossover(const chromosome<G, F>* const c, CROSSOVER_OPERATOR x) const __nonnull((1)) {
            if(this->getSize() != c->getSize()) throw std::logic_error("Crossover of chromosomes with different gene's size are now allowed!");
            typedef chromosome<G, F>    child;
            typedef child               parent;                 /* every parent was a child once:) */
            vector<child*>* children =  new vector<child*>();   /* the children container */
            auto min_size = min(this->getSize(), c->getSize()); /* get min size of parrents */
            auto dif_size = abs((long)this->getSize() - (long)c->getSize());
            switch(x) {
                case CROSSOVER_OPERATOR::XO:
                case CROSSOVER_OPERATOR::XOI:
                {
                    vector<int> p;
                    do {
                        // get 2 random break point
                        p = vector<int>({ getRand(0, (int)min_size), getRand(0, (int)min_size) });
                        // sort them to get in ASC manner
                        std::sort(p.begin(), p.end());
                        // make sure there is going to a crossover actually happen
                    } while(abs(p[0] - p[1]) < 2);
                    // make an iterate between parents
                    const parent* const parents[] = {this, c};
                    for(int pindex = 0; pindex < 2; pindex++) {
                        // fetch both parents
                        auto cpr = parents[pindex], opr = parents[(pindex + 1) % 2];
                        // create a child with size in range of [min_size, max_size]
                        child* _c = new child(min_size + (int)(dif_size * frand()), this->default_crossover_op, this->default_mutation_op);
                        // copy a sub-gene from a parrent to child accordingly
                        std::copy(cpr->_genes->begin() + p[0], cpr->_genes->begin() + p[1], _c->_genes->begin() + p[0]);
                        // branch into XO verions
                        switch(x) {
                            case CROSSOVER_OPERATOR::XO:
                            {
                                // define a smart copy flag
                                auto x = p[0] - 1;
                                // iterate on other parent's genes
                                for(size_t i = 0; i < opr->getSize(); i++) {
                                    // find a gene in other parent which is not already in child
                                    for(size_t j = p[0]; j < p[1]; j++) if(_c->gene(j) == opr->_genes->at(i)) goto __next;
                                    // check for upper-bound copy?
                                    if(x < 0) x = p[1];
                                    // for lower-bound copy
                                    if(x >= 0 && x < p[0]) _c->gene(p[0] -1 - x--) = opr->_genes->at(i);
                                    // for upper-bound copy
                                    else _c->gene(x++) = opr->_genes->at(i);
                                __next:
                                    /* next gene in other parent */;
                                }
                                break;
                            }
                            case CROSSOVER_OPERATOR::XOI:
                            {
                                // create an instance of map(for performance reasons unordered map has been used)
                                std::unordered_map<G, size_t>       indexes;
                                // store the copied genes into the child in reverse value-index mod into map
                                for(size_t i = p[0]; i < p[1]; i++) indexes.insert({_c->gene(i), i});
                                // foreach gens in other parent
                                for(size_t
                                    /* for parent */ i = 0,
                                    /* for  child */ j = 0; i < opr->getSize(); /* only parent*/ i++) {
                                    // if already copied into child?
                                    if(indexes.count(opr->gene(i))) { continue; }
                                    // if the copy location in child has reached to previously copied ones?
                                    if(j >= p[0] && j < p[1]) /* jump to other side */ j = p[1];
                                    // copy the other parent's genes into the child
                                    _c->gene(j++) = opr->gene(i);
                                }
                                break;
                            }
                            default:
                                throw std::out_of_range("Unexpected XO version");
                        }
                        // pass along the fitness function
                        _c->setFitnessFunc(this->fitFunc);
                        // push the generated child into array
                        children->push_back(_c);
                    }
                    break;
                }
                default:
                    throw std::out_of_range("Undefined crossover operation.");
            }
            return children;
        }
    };
}
#endif // CHROMOSOME_H
