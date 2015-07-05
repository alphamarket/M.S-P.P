/*
 * File:   manifest.hpp
 * Author: dariush
 *
 * Created on March 13, 2015, 21:00
 */
#ifndef MANIFEST_HPP
#define	MANIFEST_HPP
#include "hpp/teststrap.hpp"
#include "hpp/registery.hpp"
/*
 * Include test case files
 */
#include "TestCases/chromosomeTestCase.hpp"
#include "TestCases/populationTestCase.hpp"
#include "TestCases/gaConfigTestCase.hpp"
#include "TestCases/dataTestCase.hpp"
using namespace CPP_TESTER::TESTS;
namespace CPP_TESTER {
    /**
     * bootstrap the test suite for testing
     */
    void __bootstrap() {
        registery::__register("Data Tester",        new dataTestCase());
        registery::__register("GAConfig Tester",    new gaConfigTestCase());
        registery::__register("Chromosome Tester",  new chromosomeTestCase());
        registery::__register("Population Tester",  new populationTestCase());
    }
}
#endif	/* MANIFEST_HPP */
