#ifndef STDAFX_HPP
#define STDAFX_HPP

#ifdef __unix__
#   define OS_UNIX
#else
#   define OS_WIN
#   warning "The compilation of this source in any other OS than UNIX may cause many errors!"
#endif

#ifdef __GNUC__
#   define __unused __attribute__((unused))
#else
#   define __unused
#endif

#include <time.h>
#include <random>
#include <cstdlib>
using namespace std;

#define get_rand(min, range)((min) + rand() % (range))
#define frand()             (double(rand()) / RAND_MAX)
#ifdef OS_UNIX
#   include <unistd.h>
#   define dynamic_seed     (time(NULL) + rand() + getpid() * 100)
#else
#   define dynamic_seed     (time(NULL) + rand())
#endif
#define updateseed()        srand(dynamic_seed)

#ifdef __GNUC__
#	define __deprecated(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#	define __deprecated(func) __declspec(deprecated) func
#else
#	pragma message("WARNING: You need to implement `__deprecated` for this compiler!")
#	define __deprecated(func) func
#endif

#ifdef OS_UNIX
#   define clear_screen() system("clear")
#else
#   ifdef OS_WIN
#       define clear_screen() system("cls")
#   else
#       define clear_screen()
#   endif
#endif

#ifdef OS_UNIX
#   define PROMT_WAIT   "[.] "
#   define PROMT_FAILED "[\u2613] "
#   define PROMT_PASSED "[\u221A] "
#else
#   define PROMT_WAIT   "[.] "
#   define PROMT_FAILED "[!] "
#   define PROMT_PASSED "[X] "
#endif

typedef float scalar;

#endif // STDAFX_HPP
