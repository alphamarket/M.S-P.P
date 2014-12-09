/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __DIFFERENTIAL__
#define __DIFFERENTIAL__

#include <math.h>

// gets the speed of left and right wheel with their details and calculates the new location of the robot
// @param	left		The left wheel's angular velocity
// @param	right		The right wheel's angular velocity
// @param	r			The wheels' radius
// @param	l			The destance between 2 wheels
// @param	current		The current position of robot
// @param	dt			The time granity used for integral purposes
// @return	position_t	The new position of robot based on inputs
position_t differential(phidot_t left, phidot_t right, length_t r, length_t l, position_t* current, dtime_t dt){
	// define locals
	dot_t tdot, p, xdot, ydot;
	radian_t theta;
	// calc. the fix part of eq.
	phidot_t
		phi_sum = right + left,
		phi_sub = right - left;
	// if the pev. is null
	if(current == NULL){
		// assume all zero
		current = pos_new(0, 0, 0);
	}
	// calc. the \dot{\theta}
	tdot = (thetadot_t)(r * phi_sub / ( 2 * l ));
	// make the integral
	theta = tdot * dt  + current->theta;
	// make the fixed part of {x|y}dot eq.
	p = r * phi_sum / 2;
	// calc. the x dot
	xdot = p * (radian_t)cos(theta);
	// calc. the y dot
	ydot = p * (radian_t)sin(theta);
	// return an integrated new pos.
	return pos_make(xdot * dt + current->x, ydot * dt + current->y, __ANGDIFF(theta));
}
#endif