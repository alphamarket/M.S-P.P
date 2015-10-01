#ifndef CHROMOSOME_HPP
#define CHROMOSOME_HPP
#include "stdafx.hpp"
#include <assert.h>
#include <stdexcept>

class chromosome {
protected:
    double fitness;
    bool valid_fitness = false;
public:
    typedef std::vector<double> data_block;
    data_block genes;
    data_block sigmas;

    chromosome(data_block genes, data_block sigmas) {
        this->genes = genes;
        this->sigmas = sigmas;
        this->valid_fitness = false;
    }

    chromosome(data_block genes, data_block sigmas, double fitness)
        : chromosome(genes, sigmas)
    { this->set_fitness(fitness); }

    chromosome(const chromosome& c)
        : chromosome(c.genes, c.sigmas, c.fitness)
    { this->valid_fitness = c.valid_fitness; }

    void set_fitness(double fitness)
    { this->fitness = fitness; this->valid_fitness = true; }

    double get_fitness() const
    { return this->fitness; }

    bool is_fitness_valid() const
    { return this->valid_fitness; }

    bool validate(bool kill_on_error = true) const {
        bool size_match = (this->genes.size() == this->sigmas.size());
        if(kill_on_error) assert(size_match);
        else if(!size_match) throw runtime_error("genes and sigmas miss-matched!");
        return size_match;
    }
};

#endif // CHROMOSOME_HPP

