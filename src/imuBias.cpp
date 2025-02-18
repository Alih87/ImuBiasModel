#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>
#include <Eigen/Dense>
#include <cmath>
#include <tf/tf.h>
#include <iostream>

using namespace Eigen;

// Global variables
double a1 = 0.001;
double a2 = 0.001;
double a3 = 0.001;
double a4 = 0.001;

double sigma_r = 0.1;
double sigma_h = M_PI / 180.0; // 1 degree in radians
double sigma_s = 0.00001;

// Normalize an angle in the state vector x, at index `index`
void normalize_angle(Vector3d& x, int index) {
    if (x(index) > M_PI) {
        x(index) -= 2 * M_PI;
    }
    if (x(index) < -M_PI) {
        x(index) += 2 * M_PI;
    }
}

Vector3d control_update(const Vector3d& x, const Vector2d& u, double dt) {
    double v = u(0);
    double w = u(1);
    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    Vector3d x_new;
    x_new << -r * sin(x(2)) + r * sin(x(2) + w * dt),
             r * cos(x(2)) - r * cos(x(2) + w * dt),
             w * dt;

    return x + x_new;
}

void ekfloc_predict(Vector3d& x, Matrix3d& P, const Vector2d& u, double dt) {
    double h = x(2);
    double v = u(0);
    double w = u(1);

    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    double sinh = sin(h);
    double sinhwdt = sin(h + w * dt);
    double cosh = cos(h);
    double coshwdt = cos(h + w * dt);

    Matrix3d G;
    G << 1, 0.5 * (dt * dt), -r * cosh + r * coshwdt,
         0, 1, -r * sinh + r * sinhwdt,
         0, 0, 1;

    Matrix<double, 3, 2> V;
    V << (-sinh + sinhwdt) / w, v * (sin(h) - sinhwdt) / (w * w) + v * coshwdt * dt / w,
         (cosh - coshwdt) / w, -v * (cosh - coshwdt) / (w * w) + v * sinhwdt * dt / w,
         0, dt;

    Matrix2d M;
    M << a1 * v * v + a2 * w * w, 0,
         0, a3 * v * v + a4 * w * w;

    x += Vector3d(
        -r * sinh + r * sinhwdt,
         r * cosh - r * coshwdt,
         w * dt
    );
    P = G * P * G.transpose() + V * M * V.transpose();
}

void ekfloc(Vector3d& x, Matrix3d& P, const Vector2d& u, const Vector3d& z, double dt) {
    double h = x(2);
    double v = u(0);
    double w = u(1);

    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    double sinh = sin(h);
    double sinhwdt = sin(h + w * dt);
    double cosh = cos(h);
    double coshwdt = cos(h + w * dt);

    Matrix3d F;
    F << 1, 0.5 * (dt * dt), -r * cosh + r * coshwdt,
         0, 1, -r * sinh + r * sinhwdt,
         0, 0, 1;

    Matrix<double, 3, 2> V;
    V << (-sinh + sinhwdt) / w, v * (sin(h) - sinhwdt) / (w * w) + v * coshwdt * dt / w,
         (cosh - coshwdt) / w, -v * (cosh - coshwdt) / (w * w) + v * sinhwdt * dt / w,
         0, dt;

    Matrix2d M;
    M << a1 * v * v + a2 * w * w, 0,
         0, a3 * v * v + a4 * w * w;

    x += Vector3d(
        -r * sinh + r * sinhwdt,
         r * cosh - r * coshwdt,
         w * dt
    );
    P = F * P * F.transpose() + V * M * V.transpose();

    Matrix3d R = Matrix3d::Zero();
    R.diagonal() << sigma_r * sigma_r, sigma_h * sigma_h, sigma_s * sigma_s;

    Vector3d z_est = x;
    Matrix3d H = Matrix3d::Identity();
    Matrix3d S = H * P * H.transpose() + R;
    Matrix3d K = P * H.transpose() * S.inverse();

    Vector3d y = z - z_est;
    normalize_angle(y, 2);
    x += K * y;

    Matrix3d I = Matrix3d::Identity();
    Matrix3d I_KH = I - K * H;
    P = I_KH * P * I_KH.transpose() + K * R * K.transpose();
}

class ImuBias {
public:
    ImuBias(double t_start, bool use_ekf_odom = false)
        : t_start_(t_start), use_ekf_odom_(use_ekf_odom), prev_acq_time_(0.0)
    {
        // Subscribe once
        odom_sub_ = nh_.subscribe("odom", 1, &ImuBias::callback_odom, this);
        imu_sub_  = nh_.subscribe("imu/data", 1, &ImuBias::callback_imu, this);

        // (Optional) Wait for first Odom message
        try {
            ros::topic::waitForMessage<nav_msgs::Odometry>("odom", nh_, ros::Duration(1));
        } catch (const ros::Exception& e) {
            ROS_WARN("No odom message received within 1 second.");
            v_ = 0.0;
            w_ = 0.0;
            if (use_ekf_odom_) {
                odom_sub_ = nh_.subscribe("odometry/filtered", 1, &ImuBias::callback_odom, this);
            } else {
                x_ = 0.0;
            }
        }

        // Publish corrected odometry
        pub_ = nh_.advertise<nav_msgs::Odometry>("bias_corrected_odom", 1);
        // Publish predicted IMU
        pub_imu_ = nh_.advertise<sensor_msgs::Imu>("bias_corrected_imu", 1);
    }

    // IMU callback
    void callback_imu(const sensor_msgs::Imu::ConstPtr& data)
    {
        theta_ = tf::getYaw(data->orientation);
        a_ = data->linear_acceleration.y;

        double current_time = data->header.stamp.toSec();
        if (prev_acq_time_ == 0.0) {
            dt_ = 0.0; // First measurement
        } else {
            dt_ = current_time - prev_acq_time_;
        }
        acq_time_      = current_time;
        prev_acq_time_ = current_time;
    }

    // Odom callback
    void callback_odom(const nav_msgs::Odometry::ConstPtr& data)
    {
        x_ = data->pose.pose.position.x;
        y_ = data->pose.pose.position.y;
        v_ = data->twist.twist.linear.x;
        w_ = data->twist.twist.angular.z;
    }

    // Main update
    void pub_posteriori()
    {
        // Prepare state & control vectors
        X_ << x_, a_, theta_;
        U_ << v_, w_;
        P_ = Matrix3d::Identity();

        // Predict step
        Vector3d xp = X_;
        ekfloc_predict(xp, P_, U_, dt_);

        // Update step
        z_ << x_, a_, theta_;
        ekfloc(X_, P_, U_, z_, dt_);

        t_start_ = acq_time_;  // store last IMU time

        // Publish corrected odom
        publishOdom();

        // Publish predicted IMU
        publishImu();
    }

private:
    void publishOdom()
    {
        odom_.header.stamp = ros::Time::now();
        odom_.pose.pose.position.x = x_;
        odom_.pose.pose.position.y = y_;

        tf::Quaternion q = tf::createQuaternionFromYaw(X_(2));
        odom_.pose.pose.orientation.x = q.x();
        odom_.pose.pose.orientation.y = q.y();
        odom_.pose.pose.orientation.z = q.z();
        odom_.pose.pose.orientation.w = q.w();

        odom_.twist.twist.linear.x = v_;
        odom_.twist.twist.angular.z = w_;

        pub_.publish(odom_);
    }

    void publishImu()
    {
        // Construct an Imu message with EKF-predicted orientation & acceleration
        sensor_msgs::Imu imu_msg;
        imu_msg.header.stamp = ros::Time::now();

        // Orientation from X_(2) = yaw
        tf::Quaternion q = tf::createQuaternionFromYaw(X_(2));
        imu_msg.orientation.x = q.x();
        imu_msg.orientation.y = q.y();
        imu_msg.orientation.z = q.z();
        imu_msg.orientation.w = q.w();

        // We only have linear accel in X_(1) = a_
        // For a simple example, assume all accel in X direction
        imu_msg.linear_acceleration.x = 0.0;
        imu_msg.linear_acceleration.y = a_;
        imu_msg.linear_acceleration.z = 0.0;

        // (Optional) Covariances
        // imu_msg.orientation_covariance[0] = 1e-3;
        // imu_msg.angular_velocity_covariance[0] = -1.0; // no estimate
        // imu_msg.linear_acceleration_covariance[0] = 1e-3;

        pub_imu_.publish(imu_msg);
    }

    ros::NodeHandle nh_;
    // Publishers
    ros::Publisher pub_;      // bias_corrected_odom
    ros::Publisher pub_imu_;  // bias_corrected_imu

    // Subscribers
    ros::Subscriber odom_sub_;
    ros::Subscriber imu_sub_;

    double t_start_;
    bool   use_ekf_odom_;

    // Timing
    double acq_time_ = 0.0;
    double prev_acq_time_ = 0.0;
    double dt_ = 0.0;

    // EKF state
    double x_ = 0.0;
    double y_ = 0.0;
    double theta_ = 0.0;
    double a_ = 0.0;
    double v_ = 0.0;
    double w_ = 0.0;

    nav_msgs::Odometry odom_;

    Vector3d X_;
    Vector3d z_;
    Vector2d U_;
    Matrix3d P_;
};

int main(int argc, char** argv)
{
    ros::init(argc, argv, "imu_bias");
    ros::NodeHandle nh;

    double t_now = ros::Time::now().toSec() * 1e9 + ros::Time::now().toNSec();
    ImuBias imu_bias_obj(t_now);

    ros::Rate rate(30);
    while (ros::ok()) {
        // Process callbacks (IMU + Odom)
        ros::spinOnce();

        // Perform the EKF update & publish
        imu_bias_obj.pub_posteriori();

        rate.sleep();
    }
    return 0;
}
