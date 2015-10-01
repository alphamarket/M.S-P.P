#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cmath>
#include <iostream>
#include "random.hpp"

#define FLAG_MU_AND_LAMBDA      0x00
#define FLAG_MU_PLUS_LAMBDA     0x01

#define CONF_SIZE_POP           20
#define CONF_SIZE_CHILD         200
#define CONF_GEN_MAX            500
#define CONF_PROB_XOVER         0.5f
#define CONF_PROB_MUT           0.5f
#define CONF_RANGE_SIGMA        1e-5
#define CONF_GEN_REP_TYPE       FLAG_MU_AND_LAMBDA
#define CONF_DIM                30
#define CONF_BOUND_LOWER        -5
#define CONF_BOUND_UPPER        abs(CONF_BOUND_LOWER)

#define PARAM_TAW               (1 / sqrt(2 * sqrt(CONF_DIM)))
#define PARAM_TAW_PRIME         (1 / sqrt(2 * CONF_DIM))

const double CONST_STEP_GLOB =  (PARAM_TAW * normrnd());

inline void glob_print_configurations() {
    using namespace std;
#define conf(l, v) cout<<"\t"<<l<<" : "<<v<<endl
    cout<<"config: {"<<endl;
    conf("Population size", CONF_SIZE_POP);
    conf("Children   size", CONF_SIZE_CHILD);
    conf("Maximum   gen.", CONF_GEN_MAX);
    conf("Crossover prob.", CONF_PROB_XOVER);
    conf("Mutation  prob.", CONF_PROB_MUT);
    conf("Sigma     range", "[ -" << CONF_RANGE_SIGMA << " , " << CONF_RANGE_SIGMA <<" ]");
    conf("Selection type", (CONF_GEN_REP_TYPE == FLAG_MU_AND_LAMBDA ? "mu, lambda" : "mu + lambda"));
    conf("Boundries", "[ " << CONF_BOUND_LOWER << " , " << CONF_BOUND_UPPER <<" ]");
    conf("Stepsize  [taw]", PARAM_TAW);
    conf("Stepsize  [taw']", PARAM_TAW_PRIME);
    conf("Glob.     stepsize", CONST_STEP_GLOB);
    cout<<"}"<<endl;
#undef conf
}

#endif // CONFIG_HPP

