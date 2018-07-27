#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd::Zero(2 * n_aug_ + 1);
  x_ << VectorXd::Zero(n_x_);
  P_ << MatrixXd::Identity(n_x_, n_x_);
  is_initialized_ = false;
  time_us_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initialize state 
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //Convert radar from polar to cartesian coordinates
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      // initialize state.
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    }

    is_initialized_ = true;
    // no need to predict or update for the first measurement
    return;
  }

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd::Zero(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set lambda for non-augmented sigma points
  lambda_ = 3 - n_x_;

  double factor = sqrt(lambda_ + n_x_);

  Xsig.block(0,1,n_x_,n_x_)  << factor * A;
  Xsig.block(0,1+n_x_,n_x_,n_x_) << -factor * A;
  for( int i = 0; i < 2 * n_x_ + 1;++i)
    Xsig.col(i) += x_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  //set lambda for augmented sigma points
  lambda_ = 3 - n_aug_;

  //create augmented mean state
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  MatrixXd Q  = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
                   0, std_yawdd_ * std_yawdd_;
  P_aug.block(n_x_, n_x_, 2, 2) << Q;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  factor = sqrt(lambda_ + n_aug_);
  Xsig_aug.block(0, 1, n_aug_, n_aug_)  << factor * A_aug;
  Xsig_aug.block(0, 1+n_aug_, n_aug_, n_aug_) << -factor * A_aug;
  for( int i = 0; i < 2 * n_aug_ + 1;++i)
    Xsig_aug.col(i) += x_aug;

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    //extract values for better readability
    double p_x      = Xsig_aug(0, i);
    double p_y      = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else
    {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //create vector for predicted state
  VectorXd x = VectorXd::Zero(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);

  // set weights
  double weight_zero = lambda_ / (lambda_ + n_aug_);
  double weight_nonzero = 0.5 / (lambda_ + n_aug_);
  weights_(0) = weight_zero;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
    weights_(i) = weight_nonzero;

  //predict state mean
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    x = x + weights_(i) * Xsig_pred_.col(i);

  //predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // diff
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    P = P + weights_(i) * x_diff * x_diff.transpose();
  }
  P_ = P;
  x_ = x;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;

  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    z_pred = z_pred + weights_(i) * Zsig.col(i);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    //diff
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    //diff
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //diff
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z = 3;
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1,i) = atan2(p_y,p_x);
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/sqrt(p_x*p_x + p_y*p_y);
  }
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
    z_pred = z_pred + weights_(i) * Zsig.col(i);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    //diff
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    //diff
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
