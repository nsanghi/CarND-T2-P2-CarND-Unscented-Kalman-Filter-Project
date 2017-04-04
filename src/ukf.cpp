#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  //std_a_ = 30;
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.6;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  lambda_ = 3 - n_x_;
  n_aug_ = 7;
  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  double weight = 0.5/(n_aug_+lambda_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = weight;
  }

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
  /*****************************************************************************
 *  Initialization
 ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    P_ << 1.0, 0, 0, 0, 0,
          0, 1.0, 0, 0, 0,
          0, 0, 1.0, 0, 0,
          0, 0, 0, 1.0, 0,
          0, 0, 0, 0, 1.0;


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      x_ << rho * cos(phi),
            rho * sin(phi),
            rho_dot,
            phi,
            0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ <<  meas_package.raw_measurements_[0],
             meas_package.raw_measurements_[1],
             0,
             0,
             0;
    }

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    cout <<"initializing..." << endl;
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
    cout << "**********************************" << endl;

    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  if (fabs(dt > 0.00001)) {

    // process prediction in small steps as the nonlinearity will work well
    // with small dt values
    while (dt > 0.1) {
      const double delta_t = 0.05;
      Prediction(delta_t);
      dt -= delta_t;
    }
    Prediction(dt);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  cout << "**********************************" << endl;
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

  /***************************************************************************/
  // 1.0 Gnerate Augmented Sigma points
  /***************************************************************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  /***************************************************************************/
  // 2.0 Predict Sigma points
  /***************************************************************************/

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
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
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  /***************************************************************************/
  // 3.0 Predict Mean and Covariance
  /***************************************************************************/

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);


  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //assign predicated Mean and Covariance to Class state and Covariance var
  x_ = x;
  P_ = P;
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

  //actual measurement from Lidar
  VectorXd z = meas_package.raw_measurements_;

  //Using regular kalman filter equations as Laser measurement update has
  // no non-linearity

  MatrixXd H_ = MatrixXd(2, 5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  MatrixXd R_ = MatrixXd(2, 2);
  R_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  //calculate NIS of Laaser
  NIS_laser_ = y.transpose() * Si * y;

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

  /***************************************************************************/
  // 0.0 Gnerate Sigma points
  /***************************************************************************/

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  //size-5x11

  //create square root matrix
  MatrixXd L = P_.llt().matrixL();
  //size-5x5
  //cout << "L = " << L << endl;


  //create sigma points
  Xsig.col(0)  = x_;
  for (int i = 0; i< n_x_; i++)
  {
    Xsig.col(i+1)       = x_ + sqrt(lambda_+n_x_) * L.col(i);
    Xsig.col(i+1+n_x_)  = x_ - sqrt(lambda_+n_x_) * L.col(i);
  }

  //cout << "Xsig = " << Xsig << endl;
  /***************************************************************************/
  // 1.0 Predict Measurement
  /***************************************************************************/

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_x_ + 1);
  //size-3x11

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_x_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig(0,i);
    double p_y = Xsig(1,i);
    double v  = Xsig(2,i);
    double yaw = Xsig(3,i);

    if (fabs(sqrt(p_x*p_x+p_y*p_y)) < 0.0001) {
      p_x = 0.00001;
      p_y = 0.00001;
      v = 0.0;
      yaw = 0.0;
    }

    double v_x = cos(yaw)*v;
    double v_y = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v_x + p_y*v_y ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //cout << "Zsig = " << Zsig << endl;

  // weights for n_x_
  VectorXd weights = VectorXd(2*n_x_+1);
  //size-11
  weights(0) = lambda_/(lambda_+n_x_);
  double weight = 0.5/(n_x_+lambda_);
  for (int i=1; i<2*n_x_+1; i++) {  //2n+1 weights
    weights(i) = weight;
  }


  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  //size-3
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_x_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //cout << "z_pred= " << z_pred << endl;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  //size-3x3
  S.fill(0.0);
  for (int i = 0; i < 2 * n_x_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //size-3
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //cout << "S = " << S << endl;

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  //size-3x3
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //cout << "S = " << S << endl;

  /***************************************************************************/
  // 2.0 Update State
  /***************************************************************************/

  //actual measurement from Radar
  VectorXd z = meas_package.raw_measurements_;
  //size-3

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //size-5x3

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_x_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //size-3
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig.col(i) - x_;
    //size-5
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  MatrixXd Si = S.inverse();

  //cout << "Si = " << Si << endl;

  //Kalman gain K;
  MatrixXd K = Tc * Si;
  //size-5x3

  //cout << "K = " << K << endl;

  //residual
  VectorXd z_diff = z - z_pred;
  //size-3

  //cout << "z_diff = " << z_diff << endl;


  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  //size -- 5 = 5 + 5x3 * 3
  P_ = P_ - K*S*K.transpose();
  //size -- 5x5 = 5x5 - 5x3*3x3*3x5

  //Calculate NIS of Radar
  NIS_radar_ = z_diff.transpose() * Si * z_diff;

  //cout << "NIS_radar_ = " << NIS_radar_ << endl;

}
