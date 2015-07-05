/*
 * File:   tester.hpp
 * Author: dariush
 *
 * Created on March 31, 2014, 10:51 PM
 */
#ifndef TESTER_HPP
#define	TESTER_HPP
#include "teststrap.hpp"
namespace CPP_TESTER {
    class testCase {
    public:
        /**
         * Initiates the testcase
         * @return true if it execute successfully; otherwise false
         */
        virtual bool __init()  =0;
        /**
         * Disposes the testcase
         * @return true if it execute successfully; otherwise false
         */
        virtual bool __dispose()  =0;
        /**
         * Runs t the testcase
         * @return true if it execute successfully; otherwise false
         */
        virtual bool __run(int argc = 0, void** argv = NULL) = 0;
    };
}
#endif	/* TESTER_HPP */