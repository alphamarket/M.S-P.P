#define F_CPU 8000000 // 8MHz
#include "lcd.h"
#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>
#include "stdafx.h"

// typedef a interrupt ID
typedef unsigned short int INT_ID;
// the global ISR handler
void GLOB_ISR(INT_ID);
// a macro for defining a ISR for an INT0/1/2
#define DEFINE_INT_ISR(N) ISR(INT##N##_vect){ GLOB_ISR(N); }
// a macro for LCD atoi function
#define LCD_ATOI(S, I) char lcd_atoi[lcd_MaxChars]; sprintf(lcd_atoi, S, I)

// Enable ISR for INT1
DEFINE_INT_ISR(0);
// Enable ISR for INT1
DEFINE_INT_ISR(1);

// the global ISR handler
void GLOB_ISR(INT_ID __IID) {
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

/******************************* Main Program Code *************************/
int main(void)
{
// enable motor#1
    DDRA = 0b00001100;
    DDRB = 0x00;
    GICR = (1<<INT0 | 1<<INT1);								// Enable INT0/1
	MCUCR = (1<<ISC01 | 1<<ISC00 | 1<<ISC10 | 1<<ISC11);	// Trigger INT0/1 on rising edge
    sei();
    
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
