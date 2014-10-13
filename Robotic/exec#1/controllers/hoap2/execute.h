#ifndef __EXECUTE__
#define __EXECUTE__

#include <webots/robot.h>
#include <webots/motor.h>
#include <webots/touch_sensor.h>
#include <webots/emitter.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#include "utils.h"
#include "movement.h"


int execute(int argc, char *argv[]) {
	int i, sampling ;
	int com_interval;
	int pos_from_cvs[22];
	char l[500];
	int file_ended = 0;
	//int tempMotor[] = { 0,      /* body */
	//	0, 0, 0, 0, 0, 0,       /* lleg */
	//	0, 0, 0, 0, 0, 0,       /* rleg */
	//	0, 0, 0, 0, 0,          /* larm */
	//	0, 0, 0, 0, 0,          /* rarm */
	//	0, 0					/* head */
	//};                          

	/* initialize Webots */
	wb_robot_init();

	int control_step = 50;
	
	for (i = 0; i < 25; i++) {
		joint[i] = wb_robot_get_device(joint_number_to_name(i));
		wb_motor_enable_position(joint[i], control_step);
	}
    left_touch = wb_robot_get_device("left touch");
    right_touch = wb_robot_get_device("right touch");
    wb_touch_sensor_enable(left_touch, control_step);
    wb_touch_sensor_enable(right_touch, control_step);

	//tempMotor[larm_joint_5] = 0;    /* currently not needed, */
	//tempMotor[rarm_joint_5] = 0;    /* because we don't use the fingers */

	//next_position[larm_joint_5] = 0.0;
	//next_position[rarm_joint_5] = 0.0;

	//for (i = 0; i < 23; i++) {
	//	motor_position[i] = tempMotor[i] * (M_PI / 180.0) / pulse[i];
	//}

	//motor_position[larm_joint_5] = 0.0;
	//motor_position[rarm_joint_5] = 0.0;

	/* We wait a little bit before starting. */
	wb_robot_step(control_step);

	for (i = 0; i < 25; i++) {
		//next_position[i] = wb_motor_get_position(joint[i]);
		//printf("%i : %f\n", i, next_position[i]);
	}
	//next_position[24] = (double){100, 2, -2, -1, 4178, 8357, -4598, -1, 18772, -2091, -2090, 5222, -2, -1, -4180, -8360, 4595, 0, -18847, 2087, 2087, -5228, 3, 60};

	for (i=0;i<100;i++) {

		move_forward(control_step);

		double left_force = wb_touch_sensor_get_value(left_touch) / 10.0;
		double right_force = wb_touch_sensor_get_value(right_touch) / 10.0;
		double sum_force = left_force + right_force;

		printf("Touch sensors: left force: %4.1f N right force: %4.1f N -> sum: %4.1f N\n", left_force, right_force, sum_force);
	}

	return EXIT_SUCCESS;
}

#pragma GCC diagnostic pop

#endif

