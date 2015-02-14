#include "init.h"

#define __run_forever(c) while(1) c
#define sleep(m) _delay_ms(m)
#define wr_line(s, l) lcd_clear_line_4d(l); lcd_write_string_line_4d(s, l)
#define wr_line_1(s) wr_line(s, lcd_LineOne)
#define wr_line_2(s) wr_line(s, lcd_LineTwo)

/************************* Main Program Code *************************/
int main(void)
{
	__invoke_init();
    __run_forever({
    	wr_line_1("NORTH");
    	motors_move_north();
    	sleep(2500);
    	wr_line_1("SOUTH");
    	motors_move_south();
    	sleep(2500);
    	wr_line_1("EAST");
    	motors_move_east();
    	sleep(2500);
    	wr_line_1("WEST");
    	motors_move_west();
    	sleep(2500);
    	wr_line_1("TURN LEFT");
    	motors_turn_left();
    	sleep(2500);
    	wr_line_1("TURN LEFT - REV");
    	motors_turn_left_rev();
    	sleep(2500);
    	wr_line_1("TURN WRITE");
    	motors_turn_right();
    	sleep(2500);
    	wr_line_1("TURN WRITE - REV");
    	motors_turn_right_rev();
    	sleep(2500);
    });
    return 0;
}
/********************* End of Main Program Code **********************/