#include "stdafx.hpp"
#include <thread>
#include <string>
#include <unistd.h>
#include <libgen.h>
#include <signal.h>
#include "type.hpp"
using namespace std;

void initialize(char *argv[], POPULATION *p);
void generation(POPULATION *p, int gen);
void report(int gen, POPULATION *p, IPTR pop);
void statistics(POPULATION *p, IPTR pop);

volatile sig_atomic_t tsp_evol_continue = 1;
void sigint_func(int){ tsp_evol_continue = 0; }

int main(int __unused argc, char *argv[])
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
#ifndef QTCREATOR
    new std::thread([&stop](){ while(!stop) { updateseed(); usleep(25); } });
#else
    updateseed();
#endif
    argv[1]="infile";
    POPULATION population;
    POPULATION *p = &population;
    p->gen = 0;
    initialize(argv, p);
    while(tsp_evol_continue && p->gen < p->maxGen)
    {
         p->gen++;
         generation(p, p->gen);
         statistics(p, p->np);
         report(p->gen, p, p->np);
         IPTR tmp = p->op;
         p->op = p->np;
         p->np = tmp;
    }
    stop = true;
}
