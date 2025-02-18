#!/usr/bin/env python
import roslib; roslib.load_manifest('gps_nav')
import rospy, math
from nav_msgs.msg import Odometry
from sensor_msgs.msg import Imu
from math import cos, sin, sqrt, atan2
import numpy as np
from numpy import array, dot
from numpy.linalg import pinv


def print_x(x):
    print(x[0, 0], x[1, 0], np.degrees(x[2, 0]))


def control_update(x, u, dt):
    """ x is [x, a, hdg], u is [vel, omega] """

    v = u[0]
    w = u[1]
    if w == 0:
        # approximate straight line with huge radius
        w = 1.e-30
    r = v/w # radius


    return x + np.array([[-r*sin(x[2]) + r*sin(x[2] + w*dt)],
                         [ r*cos(x[2]) - r*cos(x[2] + w*dt)],
                         [w*dt]])


a1 = 0.001
a2 = 0.001
a3 = 0.001
a4 = 0.001

sigma_r = 0.1
sigma_h = a_error = np.radians(1)
sigma_s = 0.00001

def normalize_angle(x, index):
    if x[index] > np.pi:
        x[index] -= 2*np.pi
    if x[index] < -np.pi:
        x[index] = 2*np.pi


def ekfloc_predict(x, P, u, dt):

    h = x[2]
    v = u[0]
    w = u[1]

    if w == 0:
        # approximate straight line with huge radius
        w = 1.e-30
    r = v/w # radius

    sinh = sin(h)
    sinhwdt = sin(h + w*dt)
    cosh = cos(h)
    coshwdt = cos(h + w*dt)

    G = array(
       [[1, 0.5*(dt**2), -r*cosh + r*coshwdt],
        [0, 1, -r*sinh + r*sinhwdt],
        [0, 0, 1]])

    V = array(
        [[(-sinh + sinhwdt)/w, v*(sin(h)-sinhwdt)/(w**2) + v*coshwdt*dt/w],
         [(cosh - coshwdt)/w, -v*(cosh-coshwdt)/(w**2) + v*sinhwdt*dt/w],
         [0, dt]])


    # covariance of motion noise in control space
    M = array([[a1*v**2 + a2*w**2, 0],
               [0, a3*v**2 + a4*w**2]])



    x = x + array([[-r*sinh + r*sinhwdt],
                   [r*cosh - r*coshwdt],
                   [w*dt]])

    P = dot(G, P).dot(G.T) + dot(V, M).dot(V.T)

    return x, P

def ekfloc(x, P, u, z, dt):

    h = x[2]
    v = u[0]
    w = u[1]

    if w == 0:
        # approximate straight line with huge radius
        w = 1.e-30
    r = v/w # radius

    sinh = sin(h)
    sinhwdt = sin(h + w*dt)
    cosh = cos(h)
    coshwdt = cos(h + w*dt)

    F = array(
       [[1, 0.5*(dt**2), -r*cosh + r*coshwdt],
        [0, 1, -r*sinh + r*sinhwdt],
        [0, 0, 1]])

    V = array(
        [[(-sinh + sinhwdt)/w, v*(sin(h)-sinhwdt)/(w**2) + v*coshwdt*dt/w],
         [(cosh - coshwdt)/w, -v*(cosh-coshwdt)/(w**2) + v*sinhwdt*dt/w],
         [0, dt]])


    # covariance of motion noise in control space
    M = array([[a1*v**2 + a2*w**2, 0],
               [0, a3*v**2 + a4*w**2]])


    x = x + array([[-r*sinh + r*sinhwdt],
                   [r*cosh - r*coshwdt],
                   [w*dt]])

    P = dot(F, P).dot(F.T) + dot(V, M).dot(V.T)

    R = np.diag([sigma_r**2, sigma_h**2, sigma_s**2])

    z_est = x

    H = array(
        [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 1]])

    S = dot(H, P).dot(H.T) + R

    #print('S', S)
    K = dot(P, H.T).dot(pinv(S))
    y = z - z_est
    normalize_angle(y, 2)
    y = array([y]).T
    #print('y', y)
    
    # Final Estimate
    x = x + dot(K, y)
    I = np.eye(P.shape[0])
    I_KH = I - dot(K, H)
    #print('i', I_KH)

    P = dot(I_KH, P).dot(I_KH.T) + dot(K, R).dot(K.T)

    return x, P


class imuBias(object):
    def __init__(self, t_start, use_ekf_odom=False):
        self.t_start = t_start
        self.acq_time = 0
        self.use_ekf_odom = use_ekf_odom

        self.x = 0
        self.y = 0
        self.theta = 0
        self.a = 0
        self.dt = 0
        self.v, self.w = 0., 0.
        self.odom_ = Odometry()

        self.pub = rospy.Publisher('bias_corrected_odom', Odometry, queue_size=1)
    
    def quaternion_to_euler(self, x, y, z, w):
        ysqr = y * y

        t3 = +2.0 * (w * z + x * y)
        t4 = +1.0 - 2.0 * (ysqr + z * z)
        Z = math.atan2(t3, t4)

        return Z
    
    def euler_to_quaternion(self, theta):
        w = math.cos(theta/2.0)
        x = 0.
        y = 0.
        z = math.sin(theta/2.0)

        return x, y, z, w
    
    def normalize_angle(self, x):
        if x > np.pi:
            x -= 2*np.pi
        if x < -np.pi:
            x = 2*np.pi

    def get_params(self):
        '''
        If no odometer available, use local EKF odometery if specified
        '''
        try:
            rospy.Subscriber('odom', Odometry, self.callback_odom)
            msg = rospy.wait_for_message('odom', Odometry, 1)
        except rospy.exceptions.ROSException:
            self.v = 0.
            self.w = 0.
            if self.use_ekf_odom:
                rospy.Subscriber('odometry/filtered', Odometry, self.callback_odom)
            else:
                self.x = 0
        rospy.Subscriber('imu/data', Imu, self.callback_imu)

    def callback_imu(self, data):
        self.theta = self.quaternion_to_euler(data.orientation.x,
                                         data.orientation.y,
                                         data.orientation.z,
                                         data.orientation.w)
                                         
        self.a = data.linear_acceleration.x
        self.acq_time = data.header.stamp.secs * 1e9 + data.header.stamp.nsecs
        self.dt = self.acq_time - self.t_start
        # print(self.acq_time)
    
    def callback_odom(self, data):
        # These are the measured x and y coorinates of the rover
        self.x = data.pose.pose.position.x
        self.y = data.pose.pose.position.y

        # These are linear and angular input velocities
        self.v = data.twist.twist.linear.x
        self.w = data.twist.twist.angular.z

    def pub_posteriori(self):
        self.get_params()
        self.X = np.array([[self.x, self.a, self.theta]]).T
        self.U = np.array([self.v, self.w])
        self.P = np.diag([1., 1., 1.])

        xp = self.X.copy()
        xp, _ = ekfloc_predict(xp, self.P, self.U, self.dt) # Prediction
        self.X, self.P = ekfloc(self.X, self.P, self.U, xp, self.dt) # Update
        self.t_start = self.acq_time

        self.odom_.header.stamp.secs = int(self.acq_time * 1e-9)
        self.odom_.header.stamp.nsecs = self.acq_time - self.odom_.header.stamp.secs * 1e9
        self.odom_.pose.pose.position.x = self.x
        self.odom_.pose.pose.position.y = self.y

        # print("dt: ", self.dt)
        # print("Estimated theta: ", self.X[2,1])
        x,y,z,w = self.euler_to_quaternion(self.X[2,1])
        self.odom_.pose.pose.orientation.x = x
        self.odom_.pose.pose.orientation.y = y
        self.odom_.pose.pose.orientation.z = z
        self.odom_.pose.pose.orientation.w = w

        self.odom_.twist.twist.linear.x = self.v
        self.odom_.twist.twist.angular.z = self.w

        self.pub.publish(self.odom_)

if __name__ == '__main__':
    rospy.init_node('imu_bias', anonymous=False)
    rate = rospy.Rate(30)

    t_now = rospy.get_rostime().secs * 1e9 + rospy.get_rostime().nsecs
    imu_bias_obj = imuBias(t_now)

    while not rospy.is_shutdown():
        imu_bias_obj.pub_posteriori()
        rate.sleep()
        # print(rospy.get_rostime().secs * 1e9 + rospy.get_rostime().nsecs)
