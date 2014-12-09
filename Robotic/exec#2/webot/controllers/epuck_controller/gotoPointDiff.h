/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __GOTOPOINTDIFF__
#define __GOTOPOINTDIFF__

#include <math.h>
#include "Differential.h"

#define ERROR_THRESHOLD		0.05f
#define K_RHO				+3.0f
#define K_ALPHA				+8.0f
#define K_BETA				-1.5f

// linear and angular velocity configuration container
typedef struct vw{
	float v;
	float w;
} vw_t;
// returns a new linear and angular velocities based on
// @param  p	current position of robot
// @param  pg	target position of robot
// @output STOP has the robot reached the target destination?
// @return vw_t a new velocity configuration
vw_t gotoPointDiff(position_t p, position_t pg, bool* STOP) {
	vw_t $;
	length_t rho	= pos_calc_destance(p, pg);
	float __atan2	= atan2f(pg.y - p.y, pg.x - p.x);
	length_t beta	= pg.theta -__atan2;
	length_t alpha	= __atan2 - p.theta;
	*STOP			= (rho <= ERROR_THRESHOLD);
	$.v				= K_RHO * rho;
	$.w				= K_ALPHA * alpha + K_BETA * beta;
	return $;
}

#undef ERROR_THRESHOLD
#undef K_RHO
#undef K_ALPHA
#undef K_BETA

#endif