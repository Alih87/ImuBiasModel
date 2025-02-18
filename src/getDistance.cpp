#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>
#include <Eigen/Dense>
#include <cmath>
#include <tf/tf.h>
#include <iostream>

using namespace Eigen;

class calculateDistance {
    public:
        calculateDistance() {

        }
};

int main(int argc, char** argv) {
    ros::init(argc, argv, "getDistance");
    ros::NodeHandle nh;

    while(ros::ok()) {
        ros::spinOnce();

    }

    return 0;
}