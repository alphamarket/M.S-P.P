#define F_CPU 8000000 // 8MHz
#include "lcd.h"
#include <avr/io.h>
#include <util/delay.h>
#include "stdafx.h"

/******************************* Main Program Code *************************/
int main(void)
{
// enable motor#1
    DDRA = 0b00001100;
    
// initialize the LCD controller as determined by the defines (LCD instructions)
    lcd_init_4d();                                  // initialize the LCD display for a 4-bit interface

// endless loop
    while(1) {
		lcd_clear_line_4d(lcd_LineOne);
		_delay_ms(500);
    	// display the first line of information
		lcd_write_string_line_4d("Motor Clockwise", lcd_LineOne);
		// rotate motor
	    UNSET(PORTA, PIN(2));
	    SET(PORTA, PIN(3));
	    // clear LCD
		lcd_clear_line_4d(lcd_LineTwo);
		_delay_ms(500);
    	// display the first line of information
		lcd_write_string_line_4d("Counter Clockwiz", lcd_LineTwo);
		// rotate motor
	    SET(PORTA, PIN(2));
	    UNSET(PORTA, PIN(3));
    }
    return 0;
}
/******************************* End of Main Program Code ******************/
