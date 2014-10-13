#ifndef __UTILS__
#define __UTILS__

typedef enum {
    body_joint_1,
    lleg_joint_1,
    lleg_joint_3,
    lleg_joint_2,
    lleg_joint_4,
    lleg_joint_5,
    lleg_joint_6,
    rleg_joint_1,
    rleg_joint_3,
    rleg_joint_2,
    rleg_joint_4,
    rleg_joint_5,
    rleg_joint_6,
    larm_joint_1,
    larm_joint_2,
    larm_joint_3,
    larm_joint_4,
    larm_joint_5,
    rarm_joint_1,
    rarm_joint_2,
    rarm_joint_3,
    rarm_joint_4,
    rarm_joint_5,
    head_joint_2,
    head_joint_1
} joints;

static WbDeviceTag left_touch, right_touch;

//static double motor_position[] = { 0.0, /* body */
//    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /* lleg */
//    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /* rleg */
//    0.0, 0.0, 0.0, 0.0, 0.0,   /* larm */
//    0.0, 0.0, 0.0, 0.0, 0.0,   /* rarm */
//    0.0, 0.0                  /* head */
//};
//static double next_position[] = { 0.0,  /* body */
//    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /* lleg */
//    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /* rleg */
//    0.0, 0.0, 0.0, 0.0, 0.0,   /* larm */
//    0.0, 0.0, 0.0, 0.0, 0.0,   /* rarm */
//    0.0, 0.0                  /* head */
//};
static WbDeviceTag joint[25];     /* all the motors */

#define __dopen() freopen("debug", "w", stdout)
          
#define __debug(o) printf("%s\n", o); fflush(stdout)

#define __dclose() fclose(stdout)

static const char *joint_number_to_name(int num)
{
    switch (num) {
    case body_joint_1:
        return "body_joint_1";
    case lleg_joint_1:
        return "lleg_joint_1";
    case lleg_joint_3:
        return "lleg_joint_3";
    case lleg_joint_2:
        return "lleg_joint_2";
    case lleg_joint_4:
        return "lleg_joint_4";
    case lleg_joint_5:
        return "lleg_joint_5";
    case lleg_joint_6:
        return "lleg_joint_6";
    case rleg_joint_1:
        return "rleg_joint_1";
    case rleg_joint_3:
        return "rleg_joint_3";
    case rleg_joint_2:
        return "rleg_joint_2";
    case rleg_joint_4:
        return "rleg_joint_4";
    case rleg_joint_5:
        return "rleg_joint_5";
    case rleg_joint_6:
        return "rleg_joint_6";
    case larm_joint_1:
        return "larm_joint_1";
    case larm_joint_2:
        return "larm_joint_2";
    case larm_joint_3:
        return "larm_joint_3";
    case larm_joint_4:
        return "larm_joint_4";
    case larm_joint_5:
        return "larm_joint_5";
    case rarm_joint_1:
        return "rarm_joint_1";
    case rarm_joint_2:
        return "rarm_joint_2";
    case rarm_joint_3:
        return "rarm_joint_3";
    case rarm_joint_4:
        return "rarm_joint_4";
    case rarm_joint_5:
        return "rarm_joint_5";
    case head_joint_2:
        return "head_joint_2";
    case head_joint_1:
        return "head_joint_1";
    default:
        return "none";
    }
}

#endif