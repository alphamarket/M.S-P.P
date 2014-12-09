/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __DIFFERENTIAL_TEST__
#define __DIFFERENTIAL_TEST__

#include <stdlib.h>
#include <assert.h>

// define a schema for launching the test suites
#define __MAKE_DIFF_RUN(__ASSERTS) \
		p = pos_make(init_x, init_y, init_theta); \
		for(i = 0; i < TIME_STEP; i+=0.02f) { \
			position_t np = differential(left, right, (length_t)WHEEL_RADIUS_CM/100, (length_t)ROBOT_RADIUS_CM/200, &p, (dtime_t)0.02); \
			__ASSERTS \
			p = np; \
		}

// Tests differential function
// @return	int		The exit status
int __test_differential(){
	float i = 0;
	phidot_t left = 4, right = 4;
	radian_t init_theta = PI/4;
	point_t init_x = 0, init_y = 0;
	position_t p;
	// assert equal speed
	for(; left < 10 && right < 10; left++, right++) {
		__MAKE_DIFF_RUN({
			assert(p.x == p.y);
			assert(p.theta == init_theta);
			assert(np.x != init_x);
			assert(np.y != init_y);
		});
	}
	// assert hold directions
	init_theta = 0;
	__MAKE_DIFF_RUN({
		assert(init_x != np.x);
		assert(p.x < np.x);
		assert(init_y == np.y);
		assert(np.theta == init_theta);
	});
	// go arround
	left = 0;
	right = 4;
	init_theta = 0;
	init_x = init_y = 0;
	__MAKE_DIFF_RUN({
		assert(np.theta != init_theta);
		assert(np.x != np.y);
	});
	return EXIT_SUCCESS;
}

#endif