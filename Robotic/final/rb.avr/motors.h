/*
 * motors.h
 *
 *  Created on: Jan 14, 2015
 *      Author: dariush
 */

#ifndef MOTORS_H_
#define MOTORS_H_

#define MOTOR_DDR	DDRA
#define MOTOR_PORT 	PORTA

#define MOTOR_1 1
#define MOTOR_2 2
#define MOTOR_3 3
#define MOTOR_4 4

typedef unsigned short int flag;
typedef flag MOTOR;
typedef flag ROTATE_TYPE;
typedef flag bool;
typedef flag bit;

#define true 	1
#define false 	0

#define MOTOR_HALT			0b00
#define MOTOR_CLOCKWISE 			0b01
#define MOTOR_COUNTER_CLOCKWISE 	0b10

#define motor_active_1(wise) __motor_active(MOTOR_1, wise)
#define motor_active_2(wise) __motor_active(MOTOR_2, wise)
#define motor_active_3(wise) __motor_active(MOTOR_3, wise)
#define motor_active_4(wise) __motor_active(MOTOR_4, wise)

#define motor_halt(m)			__motor_active(m, MOTOR_HALT);
#define motor_active(m) 		__motor_active(m, MOTOR_CLOCKWISE);
#define motor_active_rev(m) 	__motor_active(m, MOTOR_COUNTER_CLOCKWISE);

#define __is_valid(c, m) __is_valid_##c(m)

bool __is_valid_motor(MOTOR m) { if(m > 0 && m < 5) return true; return false; }

bool __is_valid_rotate_type(ROTATE_TYPE r) { switch(r) { case 0b00: case 0b01: case 0b10: return true; default: return false; } }

void __get_motor_port_bits(MOTOR m, bit* b1, bit* b2) {
	switch(m) {
		case MOTOR_1:
			*b1 = PIN(3);
			*b2 = PIN(2);
			break;
		case MOTOR_2:
			*b1 = PIN(0);
			*b2 = PIN(1);
			break;
		case MOTOR_3:
			*b1 = PIN(6);
			*b2 = PIN(7);
			break;
		case MOTOR_4:
			*b1 = PIN(4);
			*b2 = PIN(5);
			break;
	}
}

void __motor_active(MOTOR m, ROTATE_TYPE r) {
	bit b1, b2;
	if(!__is_valid(motor, m)) return;
	if(!__is_valid(rotate_type, r)) return;
	__get_motor_port_bits(m, &b1, &b2);
	switch(r) {
		case MOTOR_HALT:
		    // halt motor
		    UNSET	(MOTOR_PORT, b1);
		    UNSET	(MOTOR_PORT, b2);
		    break;
		case MOTOR_CLOCKWISE:
		    // rotate clockwise motor
		    UNSET	(MOTOR_PORT, b1);
		    SET		(MOTOR_PORT, b2);
		    break;
		case MOTOR_COUNTER_CLOCKWISE:
		    // rotate counter clockwise motor
		    SET		(MOTOR_PORT, b1);
		    UNSET	(MOTOR_PORT, b2);
			break;
	}
}
void motors_rotate(ROTATE_TYPE r1, ROTATE_TYPE r2, ROTATE_TYPE r3, ROTATE_TYPE r4){
	motor_active_1(r1);
	motor_active_2(r2);
	motor_active_3(r3);
	motor_active_4(r4);
}
												/* MOTOR 1 */				/* MOTOR 2 */				/* MOTOR 3 */				/* MOTOR 4 */
#define motors_move_north() 	motors_rotate(	MOTOR_CLOCKWISE,			MOTOR_HALT, 				MOTOR_CLOCKWISE,			MOTOR_HALT				);
#define motors_move_east()  	motors_rotate(	MOTOR_HALT, 				MOTOR_COUNTER_CLOCKWISE, 	MOTOR_HALT, 				MOTOR_COUNTER_CLOCKWISE	);
#define motors_move_south() 	motors_rotate(	MOTOR_COUNTER_CLOCKWISE,	MOTOR_HALT, 				MOTOR_COUNTER_CLOCKWISE, 	MOTOR_HALT				);
#define motors_move_west()  	motors_rotate(	MOTOR_HALT, 				MOTOR_CLOCKWISE, 			MOTOR_HALT, 				MOTOR_CLOCKWISE			);
#define motors_turn_left()		motors_rotate(	MOTOR_CLOCKWISE,			MOTOR_HALT,					MOTOR_HALT,					MOTOR_CLOCKWISE			);
#define motors_turn_left_rev()	motors_rotate(	MOTOR_COUNTER_CLOCKWISE,	MOTOR_HALT,					MOTOR_HALT,					MOTOR_COUNTER_CLOCKWISE	);
#define motors_turn_right()		motors_rotate(	MOTOR_HALT, 				MOTOR_CLOCKWISE, 			MOTOR_CLOCKWISE,			MOTOR_HALT 				);
#define motors_turn_right_rev()	motors_rotate(	MOTOR_HALT, 				MOTOR_COUNTER_CLOCKWISE, 	MOTOR_COUNTER_CLOCKWISE,	MOTOR_HALT 				);


#endif /* MOTORS_H_ */
