/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __TEST_GOTOPOINTDIFF__
#define __TEST_GOTOPOINTDIFF__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define __GPD(p, pg) STOP = false ; $ = gotoPointDiff(p, pg, &STOP)

// Tests differential function
// @return	int		The exit status
int __test_gotoPointDiff() {
	bool STOP = false;
	int i = 0;
	vw_t $;
	// same place same angel
	__GPD(pos_make(0,0,0),pos_make(0,0,0));
	// it should stop
	assert(STOP);
	// same place different angel
	__GPD(pos_make(0,0,0),pos_make(0,0,PI));
	// it should stop
	assert(STOP);
	// if in same direction
	__GPD(pos_make(0,0,PI/2), pos_make(0,1,PI/2));
	// the angular speed should be zero
	assert($.w == 0);
	// linear speed should be greater than zero
	assert($.v > 0);
	// if same direction but different angel
	__GPD(pos_make(0,0,0), pos_make(0,1,-PI/2));
	// it should be and confirmed with matlab simulink
	assert($.w > 0);
	assert($.v > 0);
	// if same direction but different angel
	__GPD(pos_make(0,0,0), pos_make(0,-1,-PI/2));
	// it should be and confirmed with matlab simulink
	assert($.w < 0);
	assert($.v > 0);
	// same place but different angel
	return EXIT_SUCCESS;
}

#undef __GPD
#endif