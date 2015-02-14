#ifndef __STDAFX__
#define __STDAFX__

#define __PORT volatile uint8_t
#define __MASK uint8_t
#define SET(p, v) p |= v
#define UNSET(p, v) p &= ~v
#define PIN(p) _BV(p)
#define IS_SET(p, pin) bit_is_set(p, pin)
#define IS_UNSET(p, pin) bit_is_clear(p, pin)

// a macro for LCD atoi function
#define LCD_ATOI(S, I) char lcd_atoi[lcd_MaxChars]; sprintf(lcd_atoi, S, I)

#endif
