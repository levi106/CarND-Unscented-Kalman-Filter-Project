#include "ukf.h"
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  // State dimension
  n_x_ = 5;

  // Augumented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = (double)lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  // Previous timestamp
  previous_timestamp_ = 0LL;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
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
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float x = rho * cos(phi);
      float y = rho * sin(phi);
      x_ << x, y, 0, 0, 0;
    } else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
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
  MatrixXd Xsig_aug{n_aug_, 2*n_aug_+1};
  AugmentedSigmaPoints(&Xsig_aug);

  SigmaPointPrediction(Xsig_aug, delta_t);

  PredictMeanAndCovariance();
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
  MatrixXd Zsig_pred;
  TransformIntoLidarMeasurementSpace(Xsig_pred_, &Zsig_pred);

  VectorXd z_pred;
  MatrixXd S;
  PredictLidarMeasurement(Zsig_pred, &z_pred, &S);

  UpdateState(Zsig_pred, z_pred, S, meas_package.raw_measurements_, 2);
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
  MatrixXd Zsig_pred;
  TransformIntoRadarMeasurementSpace(Xsig_pred_, &Zsig_pred);

  VectorXd z_pred;
  MatrixXd S;
  PredictRadarMeasurement(Zsig_pred, &z_pred, &S);

  UpdateState(Zsig_pred, z_pred, S, meas_package.raw_measurements_, 3);
}

/**
 * Calculates augmented sigma points
 * @param Xsig_out Augmented sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  VectorXd x_aug{7};
  MatrixXd P_aug{7,7};
  MatrixXd Xsig_aug{n_aug_, 2*n_aug_+1};

  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  // write result;
  *Xsig_out = Xsig_aug;
}

/**
 * Predicts sigma points
 * @param Xsig_aug The augmented sigma points
 * @param delta_t Time between k and k+1 in s
 */
void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug, double delta_t) {
  for (int i = 0; i < 2*n_aug_+1; i++) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
* Predicts mean and covariance
*/
void UKF::PredictMeanAndCovariance() {
  // predicted state mean
  x_ = Xsig_pred_ * weights_;

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
* Transform sigma points into Lidar measurement space
* @param Xsig_pred The predicted sigma points
* @param Zsig_out The predicted sigma points in measurement space
*/
void UKF::TransformIntoLidarMeasurementSpace(MatrixXd& Xsig_pred, MatrixXd* Zsig_out) {
  MatrixXd Zsig{2, 2 * n_aug_ + 1};

  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i ++) {
    Zsig(0,i) = Xsig_pred(0,i);
    Zsig(1,i) = Xsig_pred(1,i);
  }

  *Zsig_out = Zsig;
}

/**
* Transform sigma points into Radar measurement space
* @param Xsig_pred The predicted sigma points
* @param Zsig_out The predicted sigma points in measurement space
*/
void UKF::TransformIntoRadarMeasurementSpace(MatrixXd& Xsig_pred, MatrixXd* Zsig_out) {
  MatrixXd Zsig{3, 2 * n_aug_ + 1};

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  *Zsig_out = Zsig;
}

/**
* Calculate mean predicted measurement and measurement covariance matrix
* @param Zsig_pred The predicted sigma points in measurement space
* @param z_out The mean predicted measurement
* @param S_out The measurement covariance matrix
*/
void UKF::PredictLidarMeasurement(MatrixXd& Zsig_pred, VectorXd* z_out, MatrixXd* S_out) {
  // mean predicted measurement
  VectorXd z_pred{2};

  // measurement covariance matrix S
  MatrixXd S{2, 2};

  // calculate mean predicted measurement
  z_pred = Zsig_pred * weights_;

  // calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R{2, 2};
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  S = S + R;

  *z_out = z_pred;
  *S_out = S;
}

/**
* Calculate mean predicted measurement and measurement covariance matrix
* @param Zsig_pred The predicted sigma points in measurement space
* @param z_out The mean predicted measurement
* @param S_out The measurement covariance matrix
*/
void UKF::PredictRadarMeasurement(MatrixXd& Zsig_pred, VectorXd* z_out, MatrixXd* S_out) {
  // mean predicted measurement
  VectorXd z_pred{3};

  // measurement covariance matrix S
  MatrixXd S{3, 3};

  // calculate mean predicted measurement
  z_pred = Zsig_pred * weights_;

  // calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R{3, 3};
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  S = S + R;

  *z_out = z_pred;
  *S_out = S;
}

/**
* Update state mean and covariance matrix
* @param Zsig_pred The predicted sigma points in measurement space
* @param z_pred The mean predicted measurement
* @param S The matrix for predicted measurement covariance
* @param z The incoming measurement
* @param n_z The dimension of incoming measurement
*/
void UKF::UpdateState(MatrixXd& Zsig_pred, VectorXd& z_pred, MatrixXd& S, VectorXd& z, int n_z) {
  // set measurement dimen
  MatrixXd Tc{n_x_, n_z};

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
