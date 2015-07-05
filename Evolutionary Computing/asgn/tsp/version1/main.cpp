#include "inc/Timer.hpp"
#include "inc/stdafx.hpp"
#include "main.helper.hpp"
#include <cstdio>
#include <limits>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
using namespace TSPGA;
/**
 * @brief global data instance
 */
data* _data;
/**
 * The signatures
 */
void            __intro             (void);
void            __extro             (void);
cost            costFunc            (const TSPGA::solution*   const) __nonnull();
bool            progress_init       (const TSPGA::population* const) __nonnull();
bool            progress_evolve     (const TSPGA::population* const) __nonnull();
void            analyze_evolution   (const TSPGA::population* const) __nonnull();
/**
 * @brief main entry point
 * @return exit flag
 */
int main (int argc, char** argv)
{
#ifdef __DEBUG__
    /**
     * Do debug init stuff
     */
#   include <string>
#   include <unistd.h>
    /* Change pwd to the main.cpp's directory */
    chdir(string(const_cast<char*>(__FILE__)).replace(string(const_cast<char*>(__FILE__)).find(basename(__FILE__)), string(basename(__FILE__)).length(), "").c_str());
#endif
    try{
        // init
        __intro();
        // validate and extract args
        auto vm         = validate_arg(argc, argv);
        // load json configuration file
        auto jconfig    = load_configuration(vm["config-file"].as<string>().c_str());
        // load data considering cache file
        /*^*/_data      = load_data_cached(vm);
        // build ga config instance based on json configs and loaded data
        auto gaconf     = build_gaConfig(jconfig, _data);
        // tick the timer
        Timer_sec       pTimer;
        // init. a population instance
        auto pop        = new TSPGA::population(&gaconf);
        // init. cost function
        pop ->setFitnessFunction(costFunc);
        cout<<"[ 0% ] Starting to initialize population"<<flush;
        // init. the population
        pop ->init(progress_init);
        cout<<"\r"<<__PROMT_PASSED<<" Population initialized                            "<<endl<<flush;
        // evolve the population
        const solution* const best = pop ->evolve(progress_evolve);
        cout<<"\r"<<__PROMT_PASSED<<" Evolution has done!                               "<<endl<<flush;
        // analyze the population
        analyze_evolution(pop);
        const char* outputfile = "out.res";
        cout<<"[.] rendering best solution to `"<<outputfile<<"`";
        std::ofstream os(outputfile);
        BOOST_FOREACH(size_t city, *best->_genes) {
            auto loc = _data->get(city);
            os<<loc.longitude<<" "<<loc.latitude<<endl;
        }
        os.close();
        cout<<"\r"<<__PROMT_PASSED<<" best solution rendered to `"<<outputfile<<"`"<<endl;
        // release resources
        delete pop;
        delete _data;
        // output the execution time
        cout<<endl<<"------------------------------"<<endl
            <<"Algorithm's execution took: "<<pTimer.elapsed()<<"sec."<<endl;
    } catch(const std::exception& e) {
        cerr<<endl<<"Error happen:"<<endl
            <<"------------------------"<<endl
            <<e.what()<<endl
            <<"------------------------"<<endl;
        __extro();
        return EXIT_FAILURE;
    }
    // deinit
    __extro();
    return EXIT_SUCCESS;
}
/**
 * @brief Inits env.
 */
void __intro(void) {
    system("setterm -cursor off");
}
/**
 * @brief Undoes init env.
 */
void __extro(void) {
    system("setterm -cursor on");
}
/**
 * @brief The progress function for init. the population
 * @param The population beeing init.
 * @return Should the init. continue?
 */
bool progress_init  (const TSPGA::population * const p) {
    printf("\r[%2.2f%%] Of population initialized        ", (double(p->getCurrentPopulationSize()) / p->getConfig()->getPopulation_Size()) * 100);
    fflush(stdout);
    return true;
}
/**
 * @brief The progress function for population evolution
 * @param The population beeing evolved
 * @return Should the evolution stop?
 */
bool progress_evolve(const TSPGA::population* const p) {
    cout<<"\r["<<fixed<<setprecision(2)<<double(p->getCurrentGeneration_NO())/p->getConfig()->getGenerationMaxCount() * 100<<"%] Of evolutions has done..."<<flush;
    return p->getCurrentGeneration_NO() > p->getConfig()->getGenerationMaxCount();
}
/**
 * @brief Analyzes evolution happened to a population
 * @param The population beeing evolved
 */
void analyze_evolution(const population * const p) {
    cost b = ULONG_LONG_MAX, w = -ULONG_LONG_MAX;
    for(auto it = p->getCurrentPopulation()->begin(); it != p->getCurrentPopulation()->end(); it++) {
        cost f = (*it).get()->getFitness();
        if(f < b) b = f;
        else if(f > w) w = f;
    }
    cout<<"The best found solution's fitness: "<<b<<endl;
    cout<<"The worst found solution's fitness: "<<w<<endl;
}

/**
 * @brief The cost calculator
 * @param The solution to calculate the cost for
 * @return The cost value
 */
cost costFunc(const solution* const s) {
    cost d = 0; /* aka the distance */
    /*
     * By setting [ i = -1, j = 0 ] we automatically count the distance between
     * The first and the last gene, since the path is circular tour.
     */
    for(int i = -1, j = 0; i < (int)s->getSize() && j < (int)s->getSize(); i++, j++) d += _data->getDistance(s->gene(i), s->gene(j));
    // return the cost value of current solution
    return d;
}
