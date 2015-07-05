#ifndef STDAFX_HPP
#define STDAFX_HPP

#include <time.h>
#include <unistd.h>

#define idelete(i) delete i; i = NULL;
#define adelete(a) delete[] a; a = NULL;

#define getRand(min, max)   (min) + rand() % (max)
#define frand()             (double(rand()) / RAND_MAX)
#define dynamic_seed        time(NULL) + rand() + getpid() * 100
#define updateseed()        srand(dynamic_seed)

#endif // STDAFX_HPP
