/*
 * interrupts.h
 *
 *  Created on: Jan 14, 2015
 *      Author: dariush
 */

#ifndef INTERRUPTS_H_
#define INTERRUPTS_H_
#include "lcd.h"
#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>
#include "stdafx.h"


// typedef a interrupt ID
typedef unsigned short int INT_ID;
/**
 * the global ISR handler
 */
void __glob_isr(INT_ID);

// a macro for defining a ISR for an INT0/1/2
#define DEFINE_INT_ISR(N) ISR(INT##N##_vect){ __glob_isr(N); }

// Enable ISR for INT0
DEFINE_INT_ISR(0);
// Enable ISR for INT1
DEFINE_INT_ISR(1);


// the global ISR handler
void __glob_isr(INT_ID __IID) {
	// convert __IID to string
	LCD_ATOI("INT# %d", __IID);
    // clear LCD
	lcd_clear_line_4d(lcd_LineOne);
	// display the first line of information
	lcd_write_string_line_4d(lcd_atoi, lcd_LineOne);
	_delay_ms(500);
    // clear LCD
	lcd_clear_line_4d(lcd_LineOne);
}

#endif /* INTERRUPTS_H_ */
