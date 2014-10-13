/*
 * File:         avoidance_with_lidar.c
 * Date:         24 Aug 2011
 * Description:  Example of Sick LMS 291.
 *               The velocity of each wheel is set
 *               according to a Braitenberg-like algorithm which takes the values returned by the Sick as input.
 * Author:       luc.guyot@cyberbotics.com, adapted from the original code of fabien.rohrer@cyberbotics.com
 */

#include <webots/robot.h>
#include <webots/motor.h>
#include <webots/camera.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TIME_STEP 32
#define MAX_SPEED 20
#define INIT_SPEED 10
#define DRIVING_MODE MODE_DEFAULT
#ifdef _WIN32
  #define CLEAR "cls"
#else
  #define CLEAR "clear"
#endif

WbDeviceTag front_left_wheel,
  front_right_wheel,
  back_left_wheel,
  back_right_wheel;
void set_speed(int flw, int frw, int blw, int brw) {
  if(flw != 0 && abs(flw) > MAX_SPEED)
    flw = flw * MAX_SPEED / abs(flw);
  if(frw != 0 && abs(frw) > MAX_SPEED)
    frw = frw * MAX_SPEED / abs(frw);
  if(brw != 0 && abs(brw) > MAX_SPEED)
    brw = brw * MAX_SPEED / abs(brw);
  if(blw != 0 && abs(blw) > MAX_SPEED)
    blw = blw * MAX_SPEED / abs(blw);
  wb_motor_set_velocity(front_left_wheel,  flw);
  wb_motor_set_velocity(front_right_wheel, frw);
  wb_motor_set_velocity(back_left_wheel,   blw);
  wb_motor_set_velocity(back_right_wheel,  brw);
}
void move_forward(int s) { set_speed(s, s, s, s); }
void move_backward(int s) { move_forward(-1 * abs(s)); }
void rotate_right(int s) { set_speed(s,0,s,0); }
void rotate_left(int s){ set_speed(0,s,0,s); }
void suspend() { move_forward(0); }

typedef enum move{
  FORWARD = 315,
  BACKWARD = 317,
  ROTATE_RIGHT = 316,
  ROTATE_LEFT = 314,
  NONE = 0
} move;
    
typedef enum mode {
  MODE_DEFAULT = 1,
  MODE_AUTO = 2
} mode;
    
mode invoke_mode = DRIVING_MODE;
    
int main(int argc, char **argv)
{
  // init webots stuff
  wb_robot_init();

  // get devices
  front_left_wheel  = wb_robot_get_device("front_left_wheel");
  front_right_wheel = wb_robot_get_device("front_right_wheel");
  back_left_wheel   = wb_robot_get_device("back_left_wheel");
  back_right_wheel  = wb_robot_get_device("back_right_wheel");
  
  // init motors
  wb_motor_set_position(front_left_wheel, INFINITY);
  wb_motor_set_position(front_right_wheel,INFINITY);
  wb_motor_set_position(back_left_wheel,  INFINITY);
  wb_motor_set_position(back_right_wheel, INFINITY);
  wb_robot_keyboard_enable(TIME_STEP);
  int k = FORWARD, s = INIT_SPEED;
  move move_com = NONE, prev_com = NONE;
  // control loop
  while (wb_robot_step(TIME_STEP) != -1) {
      k = wb_robot_keyboard_get_key();
      if(k == 0){
          if(invoke_mode != MODE_AUTO) {
            suspend();
            continue;
          }
          if(move_com == NONE)
            k = FORWARD;
          else continue;
        }
        switch(k) {
          case 'T': suspend(); exit(EXIT_SUCCESS);
          case 'W':
            // speed-up
            if(s < MAX_SPEED)
              s += 1;
            else continue;
            goto __SPEED_FINILIZER;
          case 'S':
            // speed-down
            if(s > 0)
              s -= 1;
            else if(s != 0) s = 0;
            else continue;
            goto __SPEED_FINILIZER;
          case ' ':
            // STOP
            s = 0;
__SPEED_FINILIZER:
            prev_com = NONE;
            break;
          // UP-ARROW key
          case (int)FORWARD:
          // RIGHT-ARROW key 
          case (int)ROTATE_RIGHT:
          // DOWN-ARROW key
          case (int)BACKWARD:
          // LEFT-ARROW key
          case (int)ROTATE_LEFT:
            prev_com = move_com;
            move_com = (move)k;
            break;
          default: continue;
        }
        // in auto-mode skip action if the command has already taken
        if(invoke_mode == MODE_AUTO && prev_com == move_com) continue;
        // apply an action
        switch(move_com){
          case FORWARD:
            move_forward(s);
            break;
          case ROTATE_RIGHT:
            if(prev_com == BACKWARD)
              rotate_right(-s);
            else
              rotate_right(s);
            break;
          case BACKWARD:
            move_backward(s);
            break;
          case ROTATE_LEFT:
            if(prev_com == BACKWARD)
              rotate_left(-s);
            else
              rotate_left(s);
            break;
          default:
            exit(EXIT_FAILURE);
        }
  }

  wb_robot_cleanup();
  
  return 0;
}
