#ifndef POPULATION_HPP
#define POPULATION_HPP
#include "stdafx.hpp"
#include <thread>
#include <time.h>
#include <vector>
#include <cstdlib>
#include <stdexcept>
using namespace std;

#include "data.hpp"
#include "utils.hpp"
#include "gaConfig.hpp"
#include "chromosome.hpp"

namespace TSPGA {
    // the population class signature[used in below typedefs]
    class population;
    typedef distance                                    cost;
    typedef chromosome<int, cost>                       solution;
    typedef unique_ptr<solution>                        solution_ptr;
    typedef vector<solution*>                           tsp_population_t;
    typedef vector<solution_ptr>                        tsp_population_ptr_t;
    typedef vector<unique_ptr<tsp_population_ptr_t>>    tsp_generations_ptr_t;
    typedef bool (*tsp_progress_report_t)               (const population* const);
    typedef tsp_progress_report_t                       tsp_termination_validator_t;
    class population
    {
        /**
         * @brief The generation#
         */
        size_t                  generation_cout = 0;
        /**
         * @brief The configuration
         */
        gaConfig                *config;
        /**
         * @brief The current population
         */
        tsp_population_ptr_t*   __population = NULL;
        /**
         * @brief The generations
         */
        tsp_generations_ptr_t*  __generation = NULL;
        /**
         * @brief The thread pool launched by current instance
         */
        vector<std::thread*>    thread_pool;
        /**
         * @brief The fitness calc. function
         */
        solution::fitnessFunc_t fitness_func;
    protected:
        /**
        * @brief clears(de-alloc) the popultions and re-alloc them
        */
        void clear_populations(void);        
        /**
         * @brief selects elit from current popultion, and store it into pass elit popultion
         */
        void op_select_elit(
            const tsp_population_ptr_t* const,
            tsp_population_ptr_t* const)                __nonnull();
        /**
        * @brief selects parent indexes based on assigned configuration
        */
        void op_select_parent(
            const tsp_population_ptr_t* const,
            size_t* const, size_t* const)               const __nonnull();
        /**
        * @brief generates a solution instance, based on passed configurations
        */
        solution* op_generate_solution()                const;
    public:
        /**
        * @brief The ctor
        */
        population(gaConfig* const);
        /**
        * @brief dtor, release every resouces that current instance has used
        */
        virtual                         ~population();
        /**
        * @brief init. the popultion
        */
        population*                     init(tsp_progress_report_t = NULL);
        /**
         * @brief performs evolutionary operations on init. population
         */
        const TSPGA::solution*          evolve(tsp_termination_validator_t);
        /**
         * @brief Get binded configuration to current population instance
         */
        inline const gaConfig*          getConfig()                 const   { return this->config; }
        /**
         * @brief Get generations evolved so for
         */
        inline tsp_generations_ptr_t*   getGenerations()            const   { return this->__generation; }
        /**
         * @brief Get current popultion in evolution
         */
        inline tsp_population_ptr_t*    getCurrentPopulation()      const   { return this->__population; }
        /**
         * @brief Get current generation number
         */
        inline size_t                   getCurrentGeneration_NO()   const   { return this->generation_cout; }
        /**
         * @brief Get current population size
         */
        inline size_t                   getCurrentPopulationSize()  const   { return this->__population->size(); }
        /**
         * @brief Set fitness fucntion to use on calc. solutions' fitness value
         */
        inline void                     setFitnessFunction(solution::fitnessFunc_t ff) __nonnull() { this->fitness_func = ff; }
        /**
         * @brief Get setted fitness calc. function
         */
        inline solution::fitnessFunc_t  getFitnessFunction()        const   { return this->fitness_func; }
    };
}

#endif // POPULATION_HPP
