#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

inline double normrnd(double mean = 0.0, double stddev = 1.0) {
    static size_t counter = 0;
    static std::mt19937 rand_gen;
    std::normal_distribution<double> norm(mean, stddev);
    if(counter++ % 1000 == 0)
        rand_gen.seed(time(0) + rand());
    return norm(rand_gen);
}

#endif // RANDOM_HPP
