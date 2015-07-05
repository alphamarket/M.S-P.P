#include "../inc/Timer.hpp"
#include "../inc/population.hpp"
#include <thread>
#include <assert.h>
#include <limits.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace TSPGA;

#define replace_population(old_, new_) while(!old_->empty()) { old_->back().reset(); old_->pop_back(); } delete old_; old_ = new_
#define fetch_conf(x, c) if(!this->config->option_exists(c)) throw std::runtime_error("Expecting option <"+std::string(c)+">, but does not exists!"); static const std::string x = boost::to_lower_copy(this->config->option(c))
#define not_anomaly(c) (c->_genes)
/**
 * @brief The ctor
 * @param gac The configuration for genetic ops
 */
population::population(TSPGA::gaConfig* gac) {
    // init the config
    this->config = gac;
    // no fitness function assigned initially
    this->fitness_func = NULL;
    // init/clear the population
    this->clear_populations();
    // create a thread to dynamically update the seed
    this->thread_pool.push_back(new std::thread([](){ while(1) { updateseed(); usleep(100); } }));
}
/**
 * @brief dtor, release every resouces that current instance has used
 */
population::~population() {
    // clear out the popultions
    this->clear_populations();
    // delete the popultions
    // because in `clear_populations` after de-allocation
    // these will be init. too(for convention)
    idelete(this->__generation);
    idelete(this->__population);
}
/**
 * @brief clears(de-alloc) the popultions and re-alloc them
 */
void population::clear_populations(void) {
    // delete and replace the popultions with NULL
    if(this->__generation) replace_population(this->__generation, NULL);
    if(this->__population) replace_population(this->__population, NULL);
    // re-alloc the popultions
    this->__population = new tsp_population_ptr_t();
    this->__generation = new tsp_generations_ptr_t();
}
/**
 * @brief init. the popultion
 * @param tpr The progress report handler
 * @return The current popultion instance
 */
population* population::init(TSPGA::tsp_progress_report_t tpr){
    // clear stuffs
    this->clear_populations();
    // calc. `gene_count_max - gene_count_min` constant for further ops.
    const size_t gene_count_diff = this->config->getGenes_Count_MAX_Limit() - this->config->getGenes_Count_MIN_Limit();
    // generate constant template solution with is a solution with <1..gene_count_max> values
    vector<int> __genes; for(size_t gene = 0; gene < this->config->getGenes_Count_MAX_Limit(); gene++) __genes.push_back(gene); const vector<int> genes = std::move(__genes);
    // while population is not fulfilled?
    while(this->getCurrentPopulationSize() < this->config->getPopulation_Size() && (tpr == NULL || tpr(this))) {
        // create a default solution
        solution* c = this->op_generate_solution();
        // assign the genes
        c->genes(&genes);
        // shuffle the default solution, to generate new solution
        std::random_shuffle(c->_genes->begin(), c->_genes->end());
        std::shuffle(c->_genes->begin(), c->_genes->end(), mt19937(dynamic_seed));
        // if gene_count_min != gene_count_max
        if(gene_count_diff)
            // randomly erase some part of current solution
            c->_genes->erase(c->_genes->end() - rand() % gene_count_diff, c->_genes->end());
        // if any fitness function defined?
        if(c->fitFunc != NULL)
            // calculate the fitness of the solution
            c->calcFitness();
        // push into current population
        this->__population->push_back(solution_ptr(c));
    }
    // we have our first generation
    this->generation_cout++;
    return this;
}
/**
 * @brief generates a solution instance, based on passed configurations
 * @return The solution instance
 */
solution* population::op_generate_solution() const {
    solution* c = new solution();
    // if any fitness function defined?
    if(this->fitness_func != NULL)
        // bind the solution to the fitness function
        c->setFitnessFunc(this->fitness_func);
    fetch_conf(method_X, "combinations.crossover.method");
    fetch_conf(method_M, "combinations.mutation.method");
    bool mached[2] = {false};
#define FOR_METHOD(m, M) if(boost::to_upper_copy(std::string(m)) == boost::to_upper_copy(std::string(M)))
#define FOR_CROSSOVER_METHOD(M, o) FOR_METHOD(method_X, M) { o; mached[0] = true; }
#define FOR_MUTATION_METHOD(M, o) FOR_METHOD(method_M, M) { o; mached[1] = true; }
    FOR_CROSSOVER_METHOD("XO", c->setDefaultCrossoverOP(CROSSOVER_OPERATOR::XO));
    FOR_CROSSOVER_METHOD("XOI", c->setDefaultCrossoverOP(CROSSOVER_OPERATOR::XOI));
    FOR_MUTATION_METHOD ("M_SWAP", c->setDefaultMutationOP(MUTATION_OPERATOR::M_SWAP));
    FOR_MUTATION_METHOD ("M_SWAP_RANGE", c->setDefaultMutationOP(MUTATION_OPERATOR::M_SWAP_RANGE));
    if(!mached[0]) throw runtime_error("Option <combinations.crossover.method> not implemented!");
    if(!mached[1]) throw runtime_error("Option <combinations.mutation.method> not implemented!");
#undef FOR_MUTATION_METHOD
#undef FOR_CROSSOVER_METHOD
#undef FOR_METHOD
    return c;
}
/**
 * @brief performs evolutionary operations on init. population
 * @param The handler to determine if the evolution should be terminated or not? (checked at every new generation)
 * @return The best found solution
 */
const solution* population::evolve(TSPGA::tsp_termination_validator_t should_terminate) {
    if(this->getCurrentPopulationSize() == 0) throw std::logic_error("The population is empty.");
    typedef solution                child;
    typedef child                   parent;     /* every parent once was a child too :) */
    typedef vector<solution*>       children;
    // for each generations
    while(!should_terminate(this) && (++this->generation_cout)) {
        /*
         [THE PSEUDO CODE]
            > elit selection into new generation
            > create new pop
                > select two solution
                > [prob] crossover the two solution, and get new solutions
                > foreach new solutions
                   > [prob] mutate them
            > replace the generatations
        */
        // init working generation
        tsp_population_ptr_t* const cur_gen = this->__population, *new_gen = new tsp_population_ptr_t;
        // extract elit popultion from current generation and insert into new generation
        this->op_select_elit(cur_gen, new_gen);
        // a lambda for checking if the new generation has reached its limit?
        auto proceed = [&]() { return new_gen->size() /*+ elit->size()*/ < this->config->getPopulation_Size(); };
        // fail-safe for population's size
        assert(cur_gen->size() == this->config->getPopulation_Size());
        // while new generation needs to be proceed
        while(proceed()) {
            // init parent indexes
            size_t pindex1 = 0, pindex2 = 0;
            // select two parent's indexes from current population
            this->op_select_parent(cur_gen, &pindex1, &pindex2);
            // obtain the selected parents
            parent *p1 = cur_gen->at(pindex1).get(), *p2 = cur_gen->at(pindex2).get();
            // fail-safe check for not-anomaly phenomena
            assert(not_anomaly(p1) && not_anomaly(p2));
            // the children container
            children* cc = NULL;
            // with a prob. crossover the parents and generate new children
            if(frand() <= this->config->getCrossover_prob()) cc = *p1 + *p2;
            // if crossover prob. failed? consider the parent as member of
            else cc = new children({ p1->clone(), p2->clone() });
            // fail-safe check point for children size
            assert(cc->size());
            // for each child
            BOOST_FOREACH(child* c, *cc) {
                // a through checkpoint for colonization process logic
                assert(p1 != c && p2 != c && c->_genes != p1->_genes && c->_genes != p2->_genes);
                // with a prob. apply for mutation
                if(frand() <= this->config->getMutation_prob()) (*c)++;
                // calculate the fitness of current child
                c->calcFitness();
                // add it to new popultion if allowed
                if(proceed()) new_gen->push_back(solution_ptr(c));
            }
            // de-alloc the children container
            idelete(cc);
        }
        // replace the gens.
        replace_population(this->__population, new_gen);
    }
    // fail-safe check point for children's size
    assert(this->__population->size());
    // return the best found solution
    return this->__population->at(0).get();
}
/**
 * @brief selects elit from current popultion, and store it into pass elit popultion
 * @note  The new population should be initialized and empty.
 * Fail to satisfy any of these constraints, the program's execution will halt by assertions!!
 * And also the elit solutions are CLONED into the new solution so manipulation in either current or new popultion
 * as far as `population::op_select_elit` concerns, won't effect each other entirely at all.
 * The function has been written in such manner that to perform very low-cost `O(n)` which `n` is number of
 * solutions in the population which elits will select from. -- That is totally AWESOME :)
 * @param current_population The population to select elits from
 * @param new_population The new population to insert the elit solutions(should always be initialized and empty)
 */
void population::op_select_elit(const tsp_population_ptr_t* const current_population, tsp_population_ptr_t* const new_population) {
    fetch_conf(elit_method, "population.elitism.method");
    static const _ratio elit_ratio = this->config->getElitism_ratio();
    /**
     * Get elit size based on popultion size and configuration
     */
    auto get_elit_size = [&](size_t pop_size) {
        if(elit_method == "%")
            return double(pop_size * elit_ratio) / 100;
        throw runtime_error("Invalid elit selction method: "+ elit_method);
    };
    // the elits size from configuration
    static const size_t e_size = get_elit_size(this->config->getPopulation_Size());
    // check point
    assert(new_population->size() == 0 && current_population->size() >= e_size);
    // typedef elit instances
    typedef const solution* elit;
    // init target population [ located on stack ]
    // this popultion will only contain the pointers of current generation's solutions
    // no need to delete its elements at the end
    vector<elit> elit_population = vector<elit>();
    // load the target elit population the worst population ever!!
    for(size_t i = 0; i < e_size; i++) { solution* s = new solution(); s->setFitness(ULONG_LONG_MAX); elit_population.push_back(s); }
    // fetch the best and worst solution in the elits' list
    elit best = elit_population.front(), worst = elit_population.back();
    // fail-safe check for worst is atleast the best
    assert(best->getFitness() <= worst->getFitness());
    // pops from back of elit list
    auto pop_back = [&]() {
        // dispose the end item
        elit_population.pop_back();
        // update the best/worst
        best = elit_population.front(), worst = elit_population.back();
    };
    // search through
    for(auto it = current_population->begin(); it != current_population->end(); it++) {
        // go no-change mode :)
        elit s = (*it).get();
        // assert on anomaly solution
        assert(not_anomaly(s));
        // don't count best/worst themselves, an also the solutions worst that the worst!!
        if(s == best || s == worst || s->getFitness() > worst->getFitness()) continue;
        // better than current best?
        if(s->getFitness() <= best->getFitness()) {
            // insert into front
            elit_population.insert(elit_population.begin(), s);
            // maintain size constraint
            pop_back();
            goto _continue;
        }
        else if(s->getFitness() == worst->getFitness()) {
            // [POLICY] tend to keep young solutions
            pop_back();
            // insert the young
            elit_population.push_back(s);
            goto _continue;
        } else if(s->getFitness() < worst->getFitness()) {
            // means `worst < s->getFitness() < best`
            // the worst has to be dismissed anyway, inorder to maintain the size constraint
            pop_back();
            // all we need to insert the `s` into a place that
            // our `best <= worst` constraint hold valid :)
            // performing a binary insertion
            elit_population.insert(std::lower_bound(elit_population.begin(), elit_population.end(), s, [](const solution* i, const solution* j) { return i->getFitness() <= j->getFitness(); }), s);
            goto _continue;
        }
        continue;
    _continue:
        continue;
    }
    // fail-safe check on size constraint
    assert(new_population->size() == 0 && elit_population.size() == e_size && elit_population.front()->getFitness() <= elit_population.back()->getFitness());
    // release memory used by elit popultion
    BOOST_FOREACH(elit s, elit_population) new_population->push_back(solution_ptr(s->clone()));
    // check point to make sure elits are emerged into the new generation
    assert(new_population->size() == elit_population.size());
}
/**
 * @brief selects parent indexes based on assigned configuration
 * @param p The population to select parent from
 * @param output_index_1 [OUTPUT] The first parent index
 * @param output_index_2 [OUTPUT] The second parent index
 */
void population::op_select_parent(const tsp_population_ptr_t* const p, size_t* const output_index_1, size_t* const output_index_2) const {
#define _case(x, y) if(x == y)
#define FOR_SELECTION_METHOD(x) _case(current_select_method, x)
#define _panic(x) throw runtime_error("Invalid setting `"+std::string(x)+"`");
    fetch_conf(current_select_method, "population.selection.method_type");
    // for tournoment selection
    FOR_SELECTION_METHOD("tournament") {
        fetch_conf(size_method, "population.selection.method_details.tournament.size_method");
        fetch_conf(ssize,       "population.selection.method_details.tournament.size");
        _case(size_method, "%") {
            // fetch the size of tournoment
            size_t size = double(boost::lexical_cast<size_t>(ssize)) / 100 * this->config->getPopulation_Size();
            // create 2 main tournoment group
            // one with smallest cost value, the other is second small cost value
            vector<pair<size_t, cost>> min1({{0, ULONG_MAX}}), min2({{0, ULONG_MAX - 1}});
            // while size of tournament
            while(size--) {
                // fetch an solution index
                size_t index = getRand(0, p->size());
                // get its cost value
                cost x = p->at(index).get()->getFitness();
                // if it is the smallest one?
                if(x < min1.at(0).second)
                    // empty the small container & add it
                    { min1.clear(); min1.push_back({index, x}); }
                // if it is among other smallest costs?
                else if(x == min1.at(0).second)
                    // add it
                    min1.push_back({index, x});
                // if it does not belong to smallest group?
                // and is smaller than secondary small group?
                else if(x < min2.at(0).second)
                    // purge them out, and place it into the second small group
                    { min2.clear(); min2.push_back({index, x}); }
                // if it is among other secondary small costs?
                else if(x == min2.at(0).second)
                    // add it
                    min2.push_back({index, x});
            }
            /**
             * Purpose of this section is to combine two small groups
             * and properly fetch(randomly) 2 index as the tournament's winners
             */
            /**
             * Phase 1) if there is only an item in at least one of the above lists
             */
            if(min1.size() == 1 && min2.size() == 1)
            { *output_index_1 = min1.at(0).first, *output_index_2 = min2.at(0).first; return; }
            if(min1.size() == 1)
            { *output_index_1 = min1.at(0).first; *output_index_2 = min2.at(getRand(0, min2.size())).first; return; }
            if(min2.size() == 1)
            { *output_index_2 = min2.at(0).first; *output_index_1 = min1.at(getRand(0, min1.size())).first; return; }
            /**
             * Phase 2) if there are more than one item in both lists
             */
            const size_t t_size = min1.size() + min2.size() - 1;
            size_t p[2];
            // fetch two index
            do{ p[0] = getRand(0, t_size), p[1] = getRand(0, t_size); } while(p[0] == p[1]);
            // assign the two combined indexes as tournament's winners
            if(p[0] < min1.size()) *output_index_1 = min1.at(p[0]).first;
            else *output_index_1 = min2.at(p[0] - min1.size()).first;
            if(p[1] < min1.size()) *output_index_2 = min1.at(p[1]).first;
            else *output_index_2 = min2.at(p[1] - min1.size()).first;
            return;
        };
        _panic("population.selection.method_details.tournament.size_method");
    };
    /* PANIC */
    throw std::runtime_error("The method `"+current_select_method+"` has no implementation!");
#undef  _panic
#undef  FOR_SELECTION_METHOD
#undef  _case
}
#undef replace_population
#undef not_anomaly
#undef fetch_conf
