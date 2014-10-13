#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std;


#include <webots/Robot.hpp>
#include <webots/DifferentialWheels.hpp>
#define TIME_STEP 10
using namespace webots;

#define echo(o) std::cout<<o<<std::endl
// All the webots classes are defined in the "webots" namespace
#define die(o) echo(o); exit(EXIT_SUCCESS)
#ifdef _WIN32
  #define CLEAR "cls"
#else
  #define CLEAR "clear"
#endif

#define DRIVING_MODE MODE_AUTO
#define INIT_SPEED 400

class epack : public DifferentialWheels {
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
    
    mode invoke_mode;
    
  protected:

    inline void move_forward(int s = 200) { this->setSpeed(s, s); }
    inline void move_backward(int s = 200) { this->move_forward(-abs(s)); }
    
  public:

    // epack constructor
    epack() : DifferentialWheels() {
      this->keyboardEnable(TIME_STEP);
      this->invoke_mode = DRIVING_MODE;
    }
    
    inline void suspend() { this->move_forward(0); }

    // epack destructor
    virtual ~epack() { this->suspend(); }

    void run() {
      // set up the speeds
      int k = FORWARD, s = INIT_SPEED;
      move move_com = NONE, prev_com = NONE;
      do {
        // read key
        k = this->keyboardGetKey();
        if(k == 0){
          if(this->invoke_mode != MODE_AUTO) {
            this->suspend();
            continue;
          }
          if(move_com == NONE)
            k = FORWARD;
          else continue;
        }
        switch(k) {
          case 'T': this->~epack(); exit(EXIT_SUCCESS);
          case 'W':
            // speed-up
            // to prevent [ Warning: wb_differential_wheels_set_speed(?, ?) overflows maxSpeed/speedUnit: 6.28/0.00628=1000 ]
            if(s < 1000)
              s += 20;
            else continue;
            goto __SPEED_FINILIZER;
          case 'S':
            // speed-down
            if(s > 20)
              s -= 20;
            else if(s != 0) s = 0;
            else continue;
            goto __SPEED_FINILIZER;
          case ' ':
            // STOP
            s = 0;
__SPEED_FINILIZER:
            prev_com = NONE;
            system(CLEAR);
            echo("Speed updated: "<<s);
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
        if(this->invoke_mode == MODE_AUTO && prev_com == move_com) continue;
        // apply an action
        switch(move_com){
          case FORWARD:
            this->move_forward(s);
            break;
          case ROTATE_RIGHT:
            if(prev_com == BACKWARD)
              this->setSpeed(-s, 0);
            else
              this->setSpeed(s, 0);
            break;
          case BACKWARD:
            this->move_backward(s);
            break;
          case ROTATE_LEFT:
            if(prev_com == BACKWARD)
              this->setSpeed(0, -s);
            else
              this->setSpeed(0, s);
            break;
          default:
            std::cerr<<"Undefined move command ["<<move_com<<" ]"<<endl;
            exit(EXIT_FAILURE);
        }
      } while (this->step(TIME_STEP) != -1);
    }	
};

// This is the main program of your controller.
// It creates an instance of your Robot subclass, launches its
// function(s) and destroys it at the end of the execution.
// Note that only one instance of Robot should be created in
// a controller program.
// The arguments of the main function can be specified by the
// "controllerArgs" field of the Robot node
int main(int argc, char **argv)
{
  epack* controller = new epack();
  controller->run();
  delete controller;
  return 0;
}
