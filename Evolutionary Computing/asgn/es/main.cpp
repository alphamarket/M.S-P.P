#include "stdafx.hpp"
#include <cmath>
#include <limits>
#include <thread>
#include <assert.h>
#include <signal.h>
#include "eStrategy.hpp"

void    glob_print_configurations();
double  eval_fitness_f1(const chromosome&);
double  eval_fitness_f2(const chromosome&);
double  eval_fitness_f3(const chromosome&);
double  eval_fitness_f4(const chromosome&);
double  eval_fitness_f5(const chromosome&);

volatile size_t iter = 0;

int main(__unused int argc, __unused char** argv) {
#ifndef QTCREATOR
    bool stop = false;
    new std::thread([&stop](){ while(!stop) { updateseed(); usleep(25); } });
#else
    updateseed();
#endif
    signal(SIGINT, [](int) { iter = CONF_GEN_MAX; });
    clear_screen();
    cout<<"Running with below configurations:"<<endl;
    glob_print_configurations();
    size_t best_iter_ever = 0;
    chromosome* best_solution = nullptr;
    eStrategy es(eval_fitness_f1);
    es.initialize();
    printf("\033[7m    Iter#    Fitness    \n\033[m");
    for(iter = 0; iter < CONF_GEN_MAX; iter++) {
        es.gen_children();
        es.survivor_selection();
        if(best_solution == nullptr ||
            best_solution->get_fitness() > es.get_population().front().get_fitness())
        {
            best_iter_ever = iter;
            if(best_solution) delete best_solution;
            best_solution = new chromosome(es.get_population().front());
        }
        printf("    %zu    %e    \n", iter, es.get_population().front().get_fitness());
    }
#ifndef QTCREATOR
    stop = true;
#endif
    best_solution->validate();
    assert(best_solution->genes.size() == CONF_DIM);
    printf("\nBest solution with fitness of `%e` found at `%zu` iteration\n", best_solution->get_fitness(), best_iter_ever);
    printf("The best found solution is described as following schema:");
    printf("\n\nGenes :");
    for(size_t index = 0; index < CONF_DIM; index++)
    { if(index % 10 == 0) printf("\n\t"); printf("%e ", best_solution->genes.at(index)); }
    printf("\n\nSigmas:");
    for(size_t index = 0; index < CONF_DIM; index++)
    { if(index % 10 == 0) printf("\n\t"); printf("%e ", best_solution->sigmas.at(index)); }
    printf("\n\nFitness:\n\t%e\n\n", best_solution->get_fitness());
}

double eval_fitness_f1(const chromosome& c) {
    // ackley
    assert(c.genes.size() == CONF_DIM);
    double sum_1 = 0, sum_2 = 0;
    for(size_t i = 0; i < c.genes.size(); i++) {
        auto g = c.genes.at(i);
        sum_1 += pow(g, 2);
        sum_2 += cos(2 * M_PI * g);
    }
	return -20 * exp(-0.2 * sqrt(sum_1 / CONF_DIM)) - exp(sum_2 / CONF_DIM) + 20 + exp(1);
}

double eval_fitness_f2(const chromosome& c) {
    assert(c.genes.size() == CONF_DIM);
    double sum_1 = 0;
    for(size_t i = 0; i < CONF_DIM; i++) {
        double sum_2 = 0;
        for(size_t j = 0; j < i; j++)
            sum_2 += c.genes.at(i);
        sum_1 += pow(sum_2, 2);
    }
    return sum_1;
}

double eval_fitness_f3(const chromosome& c) {
    assert(c.genes.size() == CONF_DIM);
    double fitv = 418.9829 * CONF_DIM;
    for(size_t i = 0; i < CONF_DIM; i++)
        fitv -= c.genes.at(i) * sin(sqrt(abs(c.genes.at(i))));
    return fitv;
}

double eval_fitness_f4(const chromosome& c) {
    assert(CONF_DIM == 2);
    assert(c.genes.size() == CONF_DIM);
    double x1 = c.genes.at(0), x2 = c.genes.at(1);
    return -abs(sin(x1) * cos(x2) * exp(abs(1 - (sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI))));
}

double eval_fitness_f5(const chromosome& c) {
    assert(c.genes.size() == CONF_DIM);
    double sum = 0;
    const int m = 10;
    for(size_t i = 0; i < CONF_DIM; i++)
        sum += sin(c.genes.at(i)) * pow(sin(i * pow(c.genes.at(i), 2) / M_PI), 2 * m);
    return -sum;
}
