#include "ukf.h"
#include "Eigen/Dense"
#include<iostream>
#include<stdio.h>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state vector dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

   // Augmented state dimension
   n_aug_ = n_x_ + 2;

   // lambda calculation for sigma points
   lambda_ = 3 - n_aug_;

   //time for first step
   time_us_ = 0.0;

   // Initialize weights_
   weights_ = VectorXd(2 * n_aug_ + 1);
   weights_.fill(0.5/(lambda_ + n_aug_));
   weights_(0) = lambda_/(lambda_ + n_aug_);

   // sigma point prediction
   Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * This function initilize the state vector and covariance and switches between lidar and radar
   * measurements based on the boolean operation
   */
   // Initilization of initial state and covariance
   if(!is_initialized_)
   {
     if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
     { // For radar Sensor
      double rho = meas_package.raw_measurements_[0]; // range
      double phi = meas_package.raw_measurements_[1]; // bearing
      double rho_dot = meas_package.raw_measurements_[2]; //velocity in radial direction

      // to convert from radial coordinates to cartesian coordinates
      double x = rho * cos(phi); // position x
      double y = rho * sin(phi); // position y
      double vx = rho_dot * cos(phi); // velocity x-component
      double vy = rho_dot * sin(phi); // velocity y-component
      double v = sqrt(vx * vx + vy * vy); // velocity magnitude
      
      // initilized state vector
      x_ << x , y, v, 0, 0;
     }
     else // initilization by Lidar measurement
     { // lidar measurement directly include x and y position
       x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
     }

     // done initializing no need to predict or update
     is_initialized_ = true;
     return;
   }

   // Covariance matrix initialize as identity matrix
   P_ << 0.06, 0, 0, 0, 0,
        0, 0.06, 0, 0, 0,
        0, 0, 0.06, 0, 0,
        0, 0, 0, 0.06, 0,
        0, 0, 0, 0, 0.06;

   // calculate dt
   double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
   time_us_ = meas_package.timestamp_;

   // Prediction step
   Prediction(dt);

   // Update step
   if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
     UpdateRadar(meas_package);
   }

   if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
     UpdateLidar(meas_package);
   }
}

void UKF::Prediction(double delta_t) {
  /**
   * this function estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

   //-----------------------------------------Generate augmented sigma points---------------------------------------------//
   // create augmented matrix
   VectorXd x_aug = VectorXd(n_aug_);

   //create augmented covariance matrix
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

   //create augmented sigma points
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

   x_aug.head(5) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;

   P_aug.fill(0.0);
   P_aug.topLeftCorner(5,5) = P_;
   P_aug(5,5) = std_a_*std_a_;
   P_aug(6,6) = std_yawdd_*std_yawdd_;

   // create square root matrix
   MatrixXd L = P_aug.llt().matrixL();
   // create augmented sigma points
   Xsig_aug.col(0) = x_aug;
   for(int i = 0; i < n_aug_; i++){
     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
   }

   //---------------------------------------------sigma point prediction--------------------------------------------------//

   for(int i = 0; i<2*n_aug_+1; i++){
     double p_x = Xsig_aug(0, i);// x-position
     double p_y = Xsig_aug(1, i);// y-position
     double v = Xsig_aug(2, i); // velocity
     double yaw = Xsig_aug(3, i); // yaw angle
     double yawd = Xsig_aug(4,i); // yaw rate
     double acc = Xsig_aug(5,i); // acceleration
     double yawdd = Xsig_aug(6,i); //angular acceleraion

     // create predicted state values
     double px_p, py_p, v_p, yaw_p, yawd_p;

     // avoid division by zero
     if(fabs(yawd) < 0.001)
     {
       px_p = p_x + v * delta_t * cos(yaw);
       py_p = p_y + v * delta_t * sin(yaw);
     } 
     else
     {
       px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
       py_p = p_y + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
     }

     //process terms
     v_p = v;
     yaw_p = yaw + yawd * delta_t;
     yawd_p = yawd;

     // adding noise terms
     px_p = px_p + 0.5 * acc * delta_t * delta_t * cos(yaw);
     py_p = py_p + 0.5 * acc * delta_t * delta_t * sin(yaw);
     v_p = v_p + acc * delta_t;
     yaw_p = yaw_p + 0.5 * yawdd * delta_t * delta_t;
     yawd_p = yawd_p + yawdd * delta_t;

     // write predicted sigma points
     Xsig_pred_(0, i) = px_p;
     Xsig_pred_(1, i) = py_p;
     Xsig_pred_(2, i) = v_p;
     Xsig_pred_(3, i) = yaw_p;
     Xsig_pred_(4, i) = yawd_p;
   }

   // --------------------------------------------------Predict state mean------------------------------------------------//
   x_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++)
   {  
     x_ = x_ + weights_(i)*Xsig_pred_.col(i);
   }

   // Predict state covairance
   P_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++)
   {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

     // angle normalization
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

     P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
   }
  }

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * This function use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
   //---------------------------------------------sigma point from measurement matrix-------------------------------------// 
   // create measurement matrix
   VectorXd z_ = meas_package.raw_measurements_;

   //measurement matrix dimensions
   int n_z_ = z_.rows();

   //sigma points in measurement domain
   MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_+1);
   for(int i = 0; i < 2*n_aug_+1; i++){
     Zsig(0, i) = Xsig_pred_(0, i);
     Zsig(1, i) = Xsig_pred_(1, i);
   }

   //-----------------------------------------Predicted mean and covariance measurement-----------------------------------//
   //predicted mean 
   VectorXd z_pred_= VectorXd(n_z_);
   z_pred_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++)
   {
     z_pred_ = z_pred_ + weights_(i) * Zsig.col(i);
   }

   // predicted covariance from measurement
   MatrixXd S = MatrixXd(n_z_, n_z_);
   S.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++)
   {
     VectorXd z_diff = Zsig.col(i) - z_pred_;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //------------------------------------------------UKF update from measurement------------------------------------------//
   // measurement noise
   MatrixXd R = MatrixXd(n_z_,n_z_);
   R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
   S = S + R;

   // Cross correlation matix
   MatrixXd Tc = MatrixXd(n_x_, n_z_);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

     VectorXd z_diff = Zsig.col(i) - z_pred_;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   // calculate Kalman gain K
   MatrixXd K = Tc * S.inverse();

   // update state mean and covariance
   // residual
   VectorXd z_diff = z_ - z_pred_;

   x_ = x_ + K*z_diff;

   P_ = P_ - K*S*K.transpose();

   //----------------------------------------------------Consistency check -----------------------------------------------//
   double NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

   //write NIS values
   std::ofstream myfile; 
   myfile.open("NIS_lidar.txt",std::ios_base::app); 

   myfile << NIS_laser_ << std::endl;
   myfile.close();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * this function use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   // Measurement prediction
   // set measurement dimension, radar can measure r, phi, and r_dot

   //---------------------------------------------sigma point from measurement matrix-------------------------------------//
   //extract measurement as VectorXd
   VectorXd z_ = meas_package.raw_measurements_;

   int n_z_ = 3;
   MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_+1);

   for(int i = 0; i < 2 * n_aug_ + 1; i++){
     double p_x = Xsig_pred_(0, i);
     double p_y = Xsig_pred_(1, i);
     double v = Xsig_pred_(2, i);
     double yaw = Xsig_pred_(3, i);
     double yawd = Xsig_pred_(4, i);

     double vx = cos(yaw)*v;
     double vy = sin(yaw)*v;

     Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                      // range
     Zsig(1, i) = atan2(p_y, p_x);                              // bearing
     Zsig(2, i) = (p_x*vx + p_y*vy)/(sqrt(p_x*p_x + p_y*p_y));  // radial velocity
   }

   //-----------------------------------------Predicted mean and covariance measurement-----------------------------------//
   // calculate mean predicted measurement
   VectorXd z_pred_ = VectorXd(n_z_);
   z_pred_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     z_pred_ = z_pred_ + weights_(i)*Zsig.col(i);
   }
   // calculate covariance of predicted measurement
   MatrixXd S = MatrixXd(n_z_, n_z_);
   S.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd z_diff = Zsig.col(i) - z_pred_;

     while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
     while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }

   //------------------------------------------------UKF update from measurement------------------------------------------//
   // measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z_,n_z_);
    R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
    S = S + R;

   // UKF update
   // Cross correlation matrixc between sigma points in state space
   // and measurement space
   MatrixXd Tc = MatrixXd(n_x_, n_z_);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

     VectorXd z_diff = Zsig.col(i) - z_pred_;
     while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
     while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   // calculate Kalman gain K
   MatrixXd K = Tc * S.inverse();

   // update state mean and covariance
   // residual
   VectorXd z_diff = z_ - z_pred_;
   while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
   while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

   //updated mean
   x_ = x_ + K*z_diff;

   //updated covariance
   P_ = P_ - K*S*K.transpose();

   //----------------------------------------------------Consistency check -----------------------------------------------//
   //calculate NIS
   double NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

   std::ofstream myfile;
   myfile.open("NIS_radar.txt",std::ios_base::app);

   myfile << NIS_radar_ << std::endl;
   myfile.close();

}
