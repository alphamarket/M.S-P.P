/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __EXECUTE__
#define __EXECUTE__

// The webot's includes
#include <webots/robot.h>
#include <webots/pen.h>
#include <webots/differential_wheels.h>

// C-Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

// include internal vendors
#include "bootstrap.h"
#include "Differential.h"
#include "gotoPointDiff.h"
#include "__test_differential.h"
#include "__test_gotoPointDiff.h"

// get phi-dots according to the inputs
// @param	$			The velocity configuration
// @param	r			The wheels' radius
// @param	l			The robot's radius
// @output	phildot		The left wheel's phi-dot value
// @output	phildot		The right wheel's phi-dot value
void get_phidots(vw_t $, length_t r, length_t l, phidot_t* const phildot, phidot_t* const phirdot) {
	*phildot = ($.v - l * $.w) / r;
	*phirdot = ($.v + l * $.w) / r;
}
// normalizes the velocities
// @input/output	phildot		The left wheel's phi-dot value
// @input/output	phildot		The right wheel's phi-dot value
void normalize_velocities(phidot_t* phildot, phidot_t* phirdot) {
#define L *phildot
#define R *phirdot
	// normalize top limit
	while(L > 1000) L /= 1.01f;
	while(R > 1000) R /= 1.01f;
	printf("\t\tDiff: %f\t\t", R - L);
#undef L
#undef R
}
// load taget position from file
// @output	x0	The initial position
// @output	xg	The target  position 
void init_pos(position_t* x0, position_t* xg){
	// open-up the file
	FILE* tf = fopen(POS_TARGET_FILE, "r");
	// create a position instance
	position_t tpos, ipos;
	point_t x, y;
	radian_t t;
	char l;
	int i;
	// define the default initial position
	ipos = pos_make(0,0,PI/2);
	// define the default target position
	tpos = pos_make(1, 3, 7);
	// if file not opened?
	if(!tf) goto __ASSIGNMENTS;
	// read two line of file
	for(i=0; tf && i<2 ;i++){
		if(fscanf(tf,"%c : { %f %f %f }\n",&l, &x, &y, &t) == 4) {
			switch(tolower(l)){
				case 'i':
					ipos = pos_make(x, y, deg2rad(t));
					break;
				case 't':
					tpos = pos_make(x, y, deg2rad(t));
					break;
				default:
					fprintf(stderr, "Invalid file content....");
					i = 2;
					break;
			}
		}
	}
	fclose(tf);
__ASSIGNMENTS:
	// assign the outputs
	*x0 = ipos;
	*xg = tpos;
	// promt the result
	printf("From ");
	__PRINT_POSITION_PTR(stdout, x0);
	printf(" going to ");
	__PRINT_POSITION_PTR(stdout, xg);
	printf("\n");
}
// An executive interface for robot controller
// @param	argc	The arguments count#
// @param	argv	The arguments values
// @return	int		The exit status
int execute(int argc, char *argv[]) {
	bool STOP = false;
	position_t x0, xg, x;
	phidot_t phildot, phirdot;
	WbDeviceTag pen;
	printf("TESTING....");
	if(__test_differential() == EXIT_FAILURE){
		fprintf(stderr, "Differential test failed cannot continue with program!");
		return EXIT_FAILURE;
	}
	if(__test_gotoPointDiff() == EXIT_FAILURE){
		fprintf(stderr, "gotoPointDiff test failed cannot continue with program!");
		return EXIT_FAILURE;
	}
	printf("\r[ OK ] TESTED\n");
	wb_robot_init();
	{
		// get the pen device
		pen = wb_robot_get_device("pen");
		// init x0 and xg
		init_pos(&x0, &xg);
		// define the initial position as current position
		x = x0;
		// drive the robot
		while(!STOP && wb_robot_step(TIME_STEP) != -1) {
			// trace the path
			wb_pen_write(pen, true);
			// print current position
			printf("Position: ");
			__PRINT_POSITION(stdout, x);
			// calc the new phi-dots value
			get_phidots(gotoPointDiff(x, xg, &STOP), (length_t)WHEEL_RADIUS_CM / 100, (length_t)ROBOT_RADIUS_CM / 100, &phildot, &phirdot);
			// if should stop?
			if(STOP) { phildot = 0; phirdot = 0; }
			// normalize the velocities
			normalize_velocities(&phildot, &phirdot);
			// print status
			printf("\t\tWheels velocity: [ %.4f, %.4f ]\n", phildot, phirdot);
			// set the wheels' velocity
			wb_differential_wheels_set_speed(phildot, phirdot);
			// get the new location of the robot
			x = differential(phildot, phirdot, (length_t)WHEEL_RADIUS_CM/100, (length_t)ROBOT_RADIUS_CM/100, &x, DT);
		}
		// stop tracing the path
		wb_pen_write(pen, false);
		printf("Execution terminated ...");
	}
	wb_robot_cleanup();
	return EXIT_SUCCESS;
}

#endif