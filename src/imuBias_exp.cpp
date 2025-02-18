#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>
#include <Eigen/Dense>
#include <cmath>
#include <tf/tf.h>
#include <iostream>

using namespace Eigen;

Matrix<double, 7, 1> x_temp_prev; 

// Global variables
double a1 = 0.001;
double a2 = 0.001;
double a3 = 0.001;
double a4 = 0.001;

double sigma_x = 0.09, sigma_y = 0.09;
double sigma_vx = 0.005, sigma_vy = 0.005;
double sigma_ax = 0.05, sigma_ay = 0.05;
double sigma_h = M_PI*0.75 / 180.0; // 0.5 degree in radians

class LowPassFilter {
private:
    double alpha;  // Filter coefficient (tuning parameter)
    double prev_output;

public:
    LowPassFilter(double cutoff_freq, double dt) {
        // Compute alpha using the RC filter formula
        double RC = 1.0 / (2 * M_PI * cutoff_freq);
        alpha = dt / (RC + dt);
        
        // Initialize previous output to 0
        prev_output = 0.0;
    }

    double apply(double input) {
        double output = alpha * input + (1 - alpha) * prev_output;
        
        // Update previous output
        prev_output = output;

        return output;
    }
};

// Normalize an angle in the state vector x, at index `index`
void normalize_angle(Matrix<double, 7, 1> & x, int index) {
    if (x(index) > M_PI) {
        x(index) -= 2 * M_PI;
    }
    if (x(index) < -M_PI) {
        x(index) += 2 * M_PI;
    }
}

double get_vy(double v, double w, double dt, double e=0.000001) {
    if(abs(w*dt) > e) {
        return v*(1-cos(w*dt))/sin(w*dt);
    }
    else {
        return 0.5*v*w*dt;
    }
}

double get_vx(double v, double e=2) {
    if(abs(v) > e || abs(v) > -e) {
        return 0;
    }
    else {
        return v;
    }
}

Matrix<double, 7, 1> control_update(const Matrix<double, 7, 1>& x, const Vector2d& u, double dt) {
    double v = u(0);
    double w = u(1);
    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    double vy = get_vy(v, w, dt);

    Matrix<double, 7, 1> x_new;
    x_new << -r * sin(x(6)) + r * sin(x(6) + w * dt),
             r * cos(x(6)) - r * cos(x(6) + w * dt),
             v,
             vy,
             v / dt,
             vy / dt,
             w * dt;

    return x + x_new;
}

void ekfloc_predict(Matrix<double, 7, 1>& x, Matrix<double, 7, 7>& P, const Vector2d& u, const Vector2d& u_prev, double dt) {
    double h = x(6);
    double v = u(0);
    double w = u(1);
    double v_x_prev = u_prev(0);
    double v_y_prev = u_prev(1);

    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    double sinh = sin(h);
    double sinhwdt = sin(h + w * dt);
    double cosh = cos(h);
    double coshwdt = cos(h + w * dt);

    Matrix<double, 7, 7> G;
    // G << 1, 0.5 * (dt * dt), -r * cosh + r * coshwdt,
    //      0, 1, -r * sinh + r * sinhwdt,
    //      0, 0, 1;
    G << 1, 0, dt, 0, 0.5*(dt*dt), 0, -r*cosh+r*coshwdt,
         0, 1, 0, dt, 0, 0.5*(dt*dt), -r*sinh+r*sinhwdt,
         0, 0, 1, 0, dt, 0, 0,
         0, 0, 0, get_vy(v, w, dt), 0, dt, 0,
         0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 1;

    Matrix<double, 7, 4> V;
    V << (-sinh + sinhwdt) / w, v * (sin(h) - sinhwdt) / (w * w) + v * coshwdt * dt / w, 0, 0,
         (cosh - coshwdt) / w, -v * (cosh - coshwdt) / (w * w) + v * sinhwdt * dt / w, 0, 0,
         1, -v*sin(h)*dt, 0, 0,
         sin(h), v*cos(h)*dt, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1,
         0, dt, 0, 0;

    Matrix<double, 4, 4> M;
    M << 10 * a1 * v * v + a2 * w * w, 0, 0, 0,
         0, a3 * v * v + a4 * w * w, 0, 0,
         0, 0, sigma_ax * sigma_ax, 0,
         0, 0, 0, sigma_ay * sigma_ay;

    v = get_vx(v);
    double vy_dt = get_vy(v, w, dt);

    Matrix<double, 7, 1> x_temp;
    x_temp << -r * sin(x(6)) + r * sin(x(6) + w * dt),
             r * cos(x(6)) - r * cos(x(6) + w * dt),
             v * dt,
             vy_dt,
             (v - v_x_prev) / dt,
             (vy_dt/dt - v_y_prev) / dt,
             w * dt;

    x += x_temp;
    P = G * P * G.transpose() + V * M * V.transpose();
}

void ekfloc(Matrix<double, 7, 1>& x, Matrix<double, 7, 7>& P, const Vector2d& u, const Vector2d& u_prev, const Matrix<double, 7, 1>& z, double dt) {
    double h = x(6);
    double v = u(0);
    double w = u(1);
    double v_x_prev = u_prev(0);
    double v_y_prev = u_prev(1);

    LowPassFilter lpf_ax(0.01, 0.1);
    LowPassFilter lpf_vx(0.1, 0.1);

    if (w == 0) {
        w = 1.e-30;
    }
    double r = v / w;

    double sinh = sin(h);
    double sinhwdt = sin(h + w * dt);
    double cosh = cos(h);
    double coshwdt = cos(h + w * dt);

    Matrix<double, 7, 7> F;
    F << 1, 0, dt, 0, 0.5*(dt*dt), 0, -r*cosh+r*coshwdt,
         0, 1, 0, dt, 0, 0.5*(dt*dt), -r*sinh+r*sinhwdt,
         0, 0, 1, 0, dt, 0, 0,
         0, 0, 0, get_vy(v, w, dt), 0, dt, 0,
         0, 0, 0, 0, 1, 0, 0,
         0, 0, 1/dt, 0, 0, 1, 0,
         0, 0, 0, 1/dt, 0, 0, 1;

    Matrix<double, 7, 4> V;
    V << (-sinh + sinhwdt) / w, v * (sin(h) - sinhwdt) / (w * w) + v * coshwdt * dt / w, 0, 0,
         (cosh - coshwdt) / w, -v * (cosh - coshwdt) / (w * w) + v * sinhwdt * dt / w, 0, 0,
         1, -v*sin(h)*dt, 0, 0,
         sin(h), v*cos(h)*dt, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1,
         0, dt, 0, 0;

    Matrix<double, 4, 4> M;
    M << a1 * v * v + a2 * w * w, 0, 0, 0,
         0, a3 * v * v + a4 * w * w, 0, 0,
         0, 0, 100 * sigma_ax * sigma_ax, 0,
         0, 0, 0, sigma_ay * sigma_ay;

    double vy_dt = get_vy(v, w, dt);
    v = lpf_vx.apply(v);
    
    Matrix<double, 7, 1> x_temp;
    if(dt >= 0.05) {
        x_temp << -r * sinh + r * sinhwdt,
         r * cosh - r * coshwdt,
         v * dt,
         get_vy(v, w, dt),
         lpf_ax.apply((v - x(2))/dt);
         (vy_dt - x(3)) / dt,
         w * dt;
    }
    x += x_temp;
    P = F * P * F.transpose() + V * M * V.transpose();

    Matrix<double, 7, 7> R = Matrix<double, 7, 7>::Zero();
    R.diagonal() << sigma_x * sigma_x, sigma_y * sigma_y, sigma_vx * 0.1,
                    sigma_vy * sigma_vy * 0.1, 2000 * sigma_ax*sigma_ax, sigma_ay*sigma_ay,
                    sigma_h * sigma_h;

    Matrix<double, 7, 1> z_est = x;
    Matrix<double, 7, 7> H = Matrix<double, 7, 7>::Identity();
    Matrix<double, 7, 7> S = H * P * H.transpose() + R;
    Matrix<double, 7, 7> K = P * H.transpose() * S.inverse();

    Matrix<double, 7, 1> y = z - z_est;
    normalize_angle(y, 6);
    x += K * y;

    // std::cout << "K\n" << K << std::endl;
    // std::cout << "y\n" << y << std::endl;

    Matrix<double, 7, 7> I = Matrix<double, 7, 7>::Identity();
    Matrix<double, 7, 7> I_KH = I - K * H;
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
        a_x = data->linear_acceleration.y;
        a_y = data->linear_acceleration.x;

        double current_time = data->header.stamp.toSec();

        if (prev_acq_time_ == 0.0) {
            dt_ = 0.0; // First measurement
        } else {
            double raw_dt_ = current_time - prev_acq_time_;
            dt_ = alpha * dt_ + (1 - alpha) * raw_dt_;
            // if(dt_ < 0.0035) {
            //     dt_ = 0.0052;
            // }
            dt_ = current_time - prev_acq_time_;
            // std::cout << dt_ << std::endl;
        }
        acq_time_ = current_time;
        prev_acq_time_ = current_time;
    }

    // Odom callback
    void callback_odom(const nav_msgs::Odometry::ConstPtr& data)
    {
        x_ = data->pose.pose.position.x;
        y_ = data->pose.pose.position.y;
        v_x = data->twist.twist.linear.x;
        v_y = data->twist.twist.linear.y;
        w_ = data->twist.twist.angular.z;
        v_ = v_x;
    }

    // Main update
    void pub_posteriori()
    {
        // Prepare state & control vectors
        X_ << x_, y_, v_x, v_y, a_x, a_y, theta_;
        U_ << v_x, w_;
        P_ = Matrix<double, 7, 7>::Identity();
        P_(2,2) = 0.01;

        synced_time = ros::Time::now();

        // Predict step
        Matrix<double, 7, 1> xp = X_;
        ekfloc_predict(xp, P_, U_, U_prev, dt_);

        // Update step
        z_ << x_, y_, v_x, v_y, a_x, a_y, theta_;
        ekfloc(X_, P_, U_, U_prev, z_, dt_);

        v_x_prev = v_x;
        v_y_prev = v_y;
        U_prev << v_x_prev, v_y_prev;
        t_start_ = acq_time_;  // store last IMU time

        // Publish corrected odom
        publishOdom();

        // Publish predicted IMU
        publishImu();
    }

private:
    void publishOdom()
    {
        odom_.header.stamp = synced_time;
        odom_.pose.pose.position.x = X_(0);
        odom_.pose.pose.position.y = X_(1);

        tf::Quaternion q = tf::createQuaternionFromYaw(X_(6));
        odom_.pose.pose.orientation.x = q.x();
        odom_.pose.pose.orientation.y = q.y();
        odom_.pose.pose.orientation.z = q.z();
        odom_.pose.pose.orientation.w = q.w();

        odom_.twist.twist.linear.x = v_x;
        odom_.twist.twist.angular.z = w_;

        pub_.publish(odom_);
    }

    void publishImu()
    {
        // Construct an Imu message with EKF-predicted orientation & acceleration
        sensor_msgs::Imu imu_msg;
        imu_msg.header.stamp = synced_time;

        // Orientation from X_(6) = yaw
        tf::Quaternion q = tf::createQuaternionFromYaw(X_(6));
        imu_msg.orientation.x = q.x();
        imu_msg.orientation.y = q.y();
        imu_msg.orientation.z = q.z();
        imu_msg.orientation.w = q.w();

        // We only have linear accel in X_(1) = a_
        // For a simple example, assume all accel in X direction
        imu_msg.linear_acceleration.x = X_(4);
        imu_msg.linear_acceleration.y = X_(5);
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

    ros::Time synced_time;

    double alpha = 0.9;

    double t_start_;
    bool use_ekf_odom_;

    // Timing
    double acq_time_ = 0.0;
    double prev_acq_time_ = 0.0;
    double dt_ = 0.0;

    // EKF state
    double x_ = 0.0;
    double y_ = 0.0;
    double theta_ = 0.0;
    double a_x = 0.0, a_y = 0.0;
    double v_x = 0.0, v_y = 0.0;
    double v_x_prev = 0.0, v_y_prev = 0.0;
    double v_ = 0.0, w_ = 0.0;

    nav_msgs::Odometry odom_;

    Matrix<double, 7, 1> X_;
    Matrix<double, 7, 1> z_;
    Vector2d U_, U_prev;
    Matrix<double, 7, 7> P_;
};

int main(int argc, char** argv)
{
    ros::init(argc, argv, "imu_bias");
    ros::NodeHandle nh;

    double t_now = ros::Time::now().toSec() * 1e9 + ros::Time::now().toNSec();
    ImuBias imu_bias_obj(t_now);

    ros::Rate rate(10);
    while (ros::ok()) {
        // Process callbacks (IMU + Odom)
        ros::spinOnce();

        // Perform the EKF update & publish
        imu_bias_obj.pub_posteriori();

        rate.sleep();
    }
    return 0;
}
