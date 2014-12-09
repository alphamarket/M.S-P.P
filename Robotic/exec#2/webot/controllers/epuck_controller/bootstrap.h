/**
* Date: Nov. 25 2014
* Author: dariush <b.g.dariush@gmail.com>
* Development/Compilation Env.: Visual Studio 2008
* Execution Platform: Webot V7.3.0
*/

#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__

// define constants
#define TIME_STEP 100
#define DT 0.001f
#define POS_TARGET_FILE "init.pos"

// define wheel radius in both cm and mm range
#define WHEEL_RADIUS_CM 2.05
#define WHEEL_RADIUS_MM 20.5

// define robots's radius used for wheel's distance from each other
#define ROBOT_RADIUS_CM 3.7
#define ROBOT_RADIUS_MM	37

// define PI#
#define PI (radian_t)3.14159265358979323846
// deg2rad * degrees = radians
#define deg2rad(d) d * (float)(PI/180.0)
// rad2deg * radians = degrees
#define rad2deg(r) r * (float)(180.0/PI)
// define a matlab angdiff
#define __ANGDIFF(_ang) (radian_t)deg2rad(fmod(rad2deg(_ang), 360.0f))

// define types
typedef float dot_t;
typedef dot_t phidot_t;
typedef dot_t dtime_t;
typedef dot_t thetadot_t;
typedef dot_t xdot_t;
typedef dot_t ydot_t;
typedef float radian_t;
typedef float length_t;
typedef float point_t;

// define a position struct
typedef struct position{
	point_t x;
	point_t y;
	radian_t theta;
} position_t;

// define a position printer
#define __SGN(x) x >= 0 ? "+" : ""
#define __PRINT_POSITION(__stdio, p) fprintf(__stdio, "[ %s%.5f %s%.5f %s%.5f ]", __SGN(p.x), p.x , __SGN(p.y), p.y, __SGN(p.theta), p.theta * (180/PI))
#define __PRINT_POSITION_PTR(__stdio, p) fprintf(__stdio, "[ %s%.5f %s%.5f %s%.5f ]", __SGN(p->x), p->x , __SGN(p->y), p->y, __SGN(p->theta), p->theta * (180/PI))

// creates new position on stack
// @param	x			The x point value
// @param	y			The y point value
// @param	t			The theta value
// @return	position_t	The created position instance on stack
position_t pos_make(point_t x, point_t y, radian_t t) {
	position_t p;
	p.x = x;
	p.y = y;
	p.theta = __ANGDIFF(t);
	return p;
}
// creates new position on heaps
// @param	x			The x point value
// @param	y			The y point value
// @param	t			The theta value
// @return	position_t	The created position instance on head
position_t* pos_new(point_t x, point_t y, radian_t t) {
	position_t* p = (position_t*)malloc(sizeof(position_t));
	p->x = x;
	p->y = y;
	p->theta = __ANGDIFF(t);
	return p;
}
// calculates the destance between 2 positions
// @param	p1			The first position
// @param	p2			The second position
// @return	length_t	The length of destance
length_t pos_calc_destance(const position_t p1, const position_t p2) {
	return (length_t)sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

#endif