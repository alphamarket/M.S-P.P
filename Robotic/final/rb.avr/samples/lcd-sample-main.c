#define F_CPU 8000000 // 8MHz
#include "lcd.h"

/******************************* Main Program Code *************************/
int main(void)
{

// initialize the LCD controller as determined by the defines (LCD instructions)
    lcd_init_4d();                                  // initialize the LCD display for a 4-bit interface

// endless loop
    while(1) {
		lcd_clear_line_4d(lcd_LineOne);
		_delay_ms(500);
    	// display the first line of information
		lcd_write_string_line_4d("LINE# 1", lcd_LineOne);
		lcd_clear_line_4d(lcd_LineTwo);
		_delay_ms(500);
    	// display the first line of information
		lcd_write_string_line_4d("LINE# 2", lcd_LineTwo);
    }
    return 0;
}
/******************************* End of Main Program Code ******************/
