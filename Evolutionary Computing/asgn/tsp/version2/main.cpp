#include "stdafx.hpp"
#include <cstdio>
#include <thread>
#include <string>
#include <math.h>
#include <cstdint>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <signal.h>
#include "deftypes.hpp"
#include "main.helper.hpp"
using namespace std;
void statistics(POPULATION*, IPTR) __nonnull();
void report(
#ifdef XFILE_COMPATIBLE
    const char*,
#endif
    int, POPULATION*, IPTR) __nonnull();
void swap_populations(POPULATION* const) __nonnull();
void generation(POPULATION* const,
#ifndef XFILE_COMPATIBLE
    const cities* const
#else
    const int** const
#endif
    , int) __nonnull();
void initialize(jsoncons::json&, POPULATION* const,
#ifndef XFILE_COMPATIBLE
    cities* const
#else
    int**&
#endif
);
volatile sig_atomic_t tsp_evol_continue = 1;
void sigint_func(int){ tsp_evol_continue = 0; }
void swap_populations(POPULATION* const p) {
    auto tmp = p->op;
    p->op = p->np;
    p->np = tmp;
    if(p->min == p->lm_lastgen_best_fitness)
        p->lm_best_fitness_repeat_count++;
    else {
        p->lm_lastgen_best_fitness = p->min;
        p->lm_best_fitness_repeat_count = 0;
    }
}
int get_rand(ulong min, ulong range) {
    static size_t z = 0;
    if(z++ % 100 == 0)
        updateseed();
    return min + rand() % range;
}
int main(int argc, char *argv[])
{
#ifdef __DEBUG__
    /**
     * Do debug init stuff
     */
    /* Change pwd to the main.cpp's directory */
    chdir(string(const_cast<char*>(__FILE__)).replace(string(const_cast<char*>(__FILE__)).find(basename(__FILE__)), string(basename(__FILE__)).length(), "").c_str());
#endif
    bool stop = false;
    setbuf(stdout, NULL);
    signal(SIGINT, sigint_func);
    new std::thread([&stop](){ while(!stop) { updateseed(); usleep(25); } });
    // validate and extract args
    auto vm         = validate_arg(argc, argv);
#ifndef XFILE_COMPATIBLE
    auto conf_file  = vm["config-file"].as<string>();
#else
    auto conf_file  = vm["input-file"].as<string>();
#endif
    // load json configuration file
    auto jconfig    = load_configuration(conf_file.c_str());
#ifndef XFILE_COMPATIBLE
    cities          data;
    cities*         data_ptr = &data;
#else
    int**           data_ptr = NULL;
#endif
    jconfig["input-file"]   = vm["input-file"].as<string>();
#ifndef XFILE_COMPATIBLE
    jconfig["output-file"]  = vm["output-file"].as<string>();
#endif
    POPULATION population;
    POPULATION *p = &population;
    p->gen = 0;
    if(jconfig["population_size"].as_int() == 0) { printf("Empty population!\n[ ABORT ]"); exit(EXIT_SUCCESS); }
    initialize(jconfig, p, data_ptr);
#ifdef XFILE_COMPATIBLE
    std::remove(jconfig["output-file"].as_string().c_str());
#endif
#ifndef STAT_ENABLED
    printf("[%2.0f%%] of generations done...", (1/(double)p->maxGen)*100);
#endif
    while(tsp_evol_continue && ++p->gen < p->maxGen){
        generation(p,
#ifndef XFILE_COMPATIBLE
            data_ptr
#else
            const_cast<const int** const>(data_ptr)
#endif
            , p->gen);
#ifdef STAT_ENABLED
        statistics(p, p->np);
        report(
#ifdef XFILE_COMPATIBLE
            jconfig["output-file"].as_string().c_str(),
#endif
            p->gen, p, p->np);
#else
        printf("\r[%2.0f%%]", (p->gen/(double)p->maxGen)*100);
#endif
        swap_populations(p);
    }
#ifndef STAT_ENABLED
    printf("\r%s all generations processed.\n", __PROMT_PASSED);
    printf("best found solution cost: %f", p->smallestEverFitness);
#endif
    std::sort(p->op, p->op + p->popSize, [](INDIVIDUAL a, INDIVIDUAL b) { return a.fitness < b.fitness; });
#ifndef XFILE_COMPATIBLE
    IPTR best_found = &p->op[0];
#endif
    ofstream out;
    out.open(jconfig["output-file"].as_string()
#ifdef XFILE_COMPATIBLE
            , std::ofstream::out | std::ofstream::app
#endif
            );
#ifdef XFILE_COMPATIBLE
    out<<"------------------------------------------"<<endl;
#endif
    for(int soluIndex = 0; soluIndex < p->popSize; soluIndex++) {
        auto solution = p->op[soluIndex];
        for(int i = 0; i < p->lchrom; i++) {
#ifndef XFILE_COMPATIBLE
            city c = data.at(best_found->chrom[i]);
            out<<c.x<<" "<<c.y<<endl;
#else
            out<<solution.chrom[i]<<" ";
#endif
        }
        out<<" : "<<solution.fitness<<endl;
    }
    out.close();
    printf("\n\n%s the TSP solution poured into `%s`.\n\n", __PROMT_PASSED, jconfig["output-file"].as_string().c_str());
    stop = true;
    exit(EXIT_SUCCESS);
}
