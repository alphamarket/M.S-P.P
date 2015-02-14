#define F_CPU 8000000 // 8MHz
#include "lcd.h"
#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>
#include "stdafx.h"

//Interrupt Service Routine for INT0
//ISR(INT0_vect)
//{
//	_delay_ms(500);
//    // clear LCD
//	lcd_clear_line_4d(lcd_LineOne);
//	// display the first line of information
//	lcd_write_string_line_4d("Interupt", lcd_LineOne);
//}

/******************************* Main Program Code *************************/
int main(void)
{
// enable motor#1
    DDRA = 0b00001100;
    DDRB = 0b11111011;

	char c[lcd_MaxChars];

//    GICR = /*1<<INT0 | */1<<INT2;					// Enable INT0
//	//MCUCR = 1<<ISC01 | 1<<ISC00;
//	MCUCSR &= 0<<ISC2;
//    sei();
    
// initialize the LCD controller as determined by the defines (LCD instructions)
    lcd_init_4d();                                  // initialize the LCD display for a 4-bit interface
// endless loop
    while(1) {
		lcd_clear_line_4d(lcd_LineOne);
		_delay_ms(500);
		//
		memset(c, ' ', lcd_MaxChars);
		c[lcd_MaxChars] = '\0';
		//sprintf(c, "%x", PORTB);
    	// display the first line of information
		lcd_write_string_line_4d(IS_SET(PORTB, PIN(2)) ? "1" : "0", lcd_LineOne);
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
