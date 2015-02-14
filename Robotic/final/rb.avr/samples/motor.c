#define F_CPU 8000000 // 8MHz
#include <avr/io.h>
#include <util/delay.h>
#include "stdafx.h"

int main(void)
{
    DDRA = 0x0B;
    DDRB = 0xff*0;
    DDRC = 0xff*0;
    DDRD = 0xff*0;
	while(1)
	{
	    UNSET(PORTA, PIN(2));
	    SET(PORTA, PIN(3));
	    _delay_ms(1000);
	    SET(PORTA, PIN(2));
	    UNSET(PORTA, PIN(3));
	    _delay_ms(1000);
	}
	return 0;
}
