/*
 * init.h
 *
 *  Created on: Jan 14, 2015
 *      Author: dariush
 */

#ifndef INIT_H_
#define INIT_H_

#define F_CPU 8000000 // 8MHz
#include "lcd.h"
#include <stdlib.h>
#include <ctype.h>
#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>
#include "stdafx.h"
#include "interrupts.h"
#include "motors.h"
/**
 * general init handler
 */
void __invoke_init(void);
/**
 * ISR initializer
 */
void __init_isr(void);
/**
 * motors initializer
 */
void __init_motors(void);
/**
 * general init handler
 */
void __invoke_init(void){
	__init_isr();
	__init_motors();
	lcd_init_4d();
}
/**
 * ISR initializer
 */
void __init_isr(void) {
	GICR = (1<<INT0 | 1<<INT1);								// Enable INT0/1
	MCUCR = (1<<ISC01 | 1<<ISC00 | 1<<ISC10 | 1<<ISC11);	// Trigger INT0/1 on rising edge
	sei();													// enable interrupts
}
/**
 * motors initializer
 */
void __init_motors(void) {
	MOTOR_DDR = 0xFF;			// motors' pins are output
}

#endif /* INIT_H_ */
