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
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  /**  TODO:Complete the initialization. See ukf.h for other member properties. */
  //Set is_initialized to false.
  is_initialized_ = false;

  //Set time_us to zero.
  time_us_ = 0;
  //Size of state vector 5
  n_x_ = 5;
  //Size of state augmented vector 7
  n_aug_= 7;
  //Spread parameter lambda.
  lambda_ = (3 - n_aug_);
  // Number of columns for augmented sigma points
  int aug_cols = (2* n_aug_ + 1);
  // Define weights values
  weights_ = VectorXd(aug_cols);
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < aug_cols; i++)
  {
    weights_[i] = 1/(2*(lambda_ + n_aug_));
  }
  //Define Xsig_pred_ matrix size
  Xsig_pred_ = MatrixXd(n_x_, aug_cols);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //Define delta_time
  double delta_time;
  //std::cout <<"Sensor Type: "<< meas_package.sensor_type_<<std::endl;
  //std::cout << "Data received from sensor: " << meas_package.raw_measurements_<< std::endl;

  // Initialize state vector   LASER / RADAR
  if(!is_initialized_){
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER){
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
    }
    else {
      x_[0] = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      x_[1] = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  //Calculate Delta Time for each follow measurements.
  delta_time = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  //std::cout<<"Delta time: \n"<<delta_time<<std::endl;

  //Prediction step, function call .
  UKF::Prediction(delta_time);

  // Update step LASER/RADAR.
  if(meas_package.sensor_type_== MeasurementPackage::SensorType::LASER)
  {
    UKF::UpdateLidar(meas_package);
  } else {
    UKF::UpdateRadar(meas_package);
  }


}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //Generate Augmented Sigma Points.
  VectorXd x_aug = VectorXd(7);
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  //std::cout << " x aumentada: \n " << x_aug << std::endl;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //std::cout << "P Augmentada:  \n" << P_aug << std::endl;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;

  P_aug = P_aug * (lambda_ + n_aug_);

  // calculate square root P_aug.
  MatrixXd A = MatrixXd(n_aug_,n_aug_);
  A.fill(0.0);
  A = P_aug.llt().matrixL();
  //std::cout << "L_Pa:  \n" << A << std::endl;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - A.col(i);
  }
  //std::cout << "Xsig_aug  \n" << Xsig_aug << std::endl;

  // Calculate Xsig_pred_
  VectorXd coln = VectorXd(5);
  coln.fill(0.0);
  Xsig_pred_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    double Vk = Xsig_aug(2, i);
    double Yk = Xsig_aug(3, i);
    double YDk = Xsig_aug(4, i);
    double Vak = Xsig_aug(5, i);
    double YDDk = Xsig_aug(6, i);
    double Pxp;
    double Pyp;
    double vp;
    double Yp;
    double YDp;


    if (fabs(YDk) > 0.0001) {
      Pxp = Xsig_aug(0, i) + ((Vk / YDk) * (sin(Yk + (YDk * delta_t)) - sin(Yk))) +
            (0.5 * (delta_t * delta_t) * cos(Yk) * Vak);
      Pyp = Xsig_aug(1, i) + (Vk / YDk) * (cos(Yk) - cos(Yk + (YDk * delta_t))) +
            (0.5 * (delta_t * delta_t) * sin(Yk) * Vak);
      vp = Xsig_aug(2, i) + delta_t * Vak;
      Yp = Xsig_aug(3, i) + YDk * delta_t + (0.5 * (delta_t * delta_t)) * YDDk;
      YDp = Xsig_aug(4, i) + delta_t * YDDk;

    } else {
      Pxp = Xsig_aug(0, i) + (Vk * cos(Yk) * delta_t) + (0.5 * (delta_t * delta_t) * cos(Yk) * Vak);
      Pyp = Xsig_aug(1, i) + (Vk * (sin(Yk) * delta_t)) + (0.5 * (delta_t * delta_t) * sin(Yk) * Vak);
      vp = Xsig_aug(2, i) + delta_t * Vak;
      Yp = Xsig_aug(3, i) + (0.5 * (delta_t * delta_t)) * YDDk;
      YDp = Xsig_aug(4, i) + delta_t * YDDk;
    }

    coln << Pxp, Pyp, vp, Yp, YDp;

    Xsig_pred_.col(i) = coln;
  }
  //std::cout << "Xsig_pred_  \n" << Xsig_pred_ << std::endl;
  //std::cout << "weights \n" << weights_ << std::endl;
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //iterate over sigma points
    x_ +=  weights_(i) * Xsig_pred_.col(i);
    //std::cout << "x_ predicted \n" << x_<< std::endl;
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //std::cout << "x_ predicted \n" << x_<< std::endl;
  //std::cout << "P_ predicted  \n" << P_ << std::endl;

}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  int n_z = 2;
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  MatrixXd z_sig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_sig.fill(0.0);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  VectorXd coln = VectorXd(5);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    coln = Xsig_pred_.col(i);
    z_sig.col(i)[0] = coln[0];
    z_sig.col(i)[1] = coln[1];
    z_pred += weights_[i] * z_sig.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    S += weights_[i] * ((z_sig.col(i) - z_pred) * (z_sig.col(i) - z_pred).transpose());
  }
  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  VectorXd x_diff = VectorXd(2);
  for (int i = 0; i < n_aug_ * 2 + 1; i++) {
    x_diff = (Xsig_pred_.col(i) - x_);

    Tc += weights_[i] * (x_diff * (z_sig.col(i) - z_pred).transpose());
  }
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  VectorXd z_diff = VectorXd(n_z);
  z_diff = z - z_pred;
  double E;
  E = z_diff.transpose() * S.inverse() * z_diff;
  //std::cout<<"EL "<<E <<"\n";
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
  int n_z = 3;
  int n_x = 5;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  VectorXd z = VectorXd(n_z);
  MatrixXd Tc = MatrixXd(n_x, n_z);

  z<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  //std::cout<< "Z measurement " << z<<"\n";
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double Px = Xsig_pred_(0, i);
    double Py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double Vp1 = sqrt(Px * Px + Py * Py);
    double ph1 = atan2(Py, Px);///formula
    if(fabs(Vp1) < 0.00001){
      std::cout<<"\n Zero Division  !!!!\n";//Vp1 = 0.0001;
    }
    double phD = (Px * cos(yaw) * v + Py * sin(yaw) * v) / Vp1;

    Zsig.col(i)<< Vp1, ph1, phD;
  }
  //std::cout<<"Zsig"<< Zsig<<"\n";
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i =0; i < 15; i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
  //std::cout<< "z_pred \n"<<z_pred<<"\n";

  //calculate measurement covariance matrix S

  S.fill(0.0);
  S(0,0) = (std_radr_ * std_radr_);
  S(1,1) = (std_radphi_ * std_radphi_);
  S(2,2) = (std_radrd_ * std_radrd_);

  VectorXd resta = VectorXd(3);
  resta.fill(0.0);
  for(int i = 0; i < 15; i++)
  {
    resta = Zsig.col(i) - z_pred;
    S += weights_(i) * resta * resta.transpose();
  }
  Tc.fill(0.0);
  //calculate cross correlation matrix
  VectorXd Xdiff = VectorXd(5);
  Xdiff.fill(0.0);
  VectorXd Zdiff = VectorXd(3);
  Zdiff.fill(0.0);
  for (int i = 0; i < 15; i++)
  {
    Xdiff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (Xdiff(3)> M_PI) Xdiff(3)-=2.*M_PI;
    while (Xdiff(3)<-M_PI) Xdiff(3)+=2.*M_PI;

    Zdiff = Zsig.col(i) - z_pred;
    //angle normalization
    while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
    while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;

    Tc = Tc + weights_(i) * Xdiff * Zdiff.transpose();
  }
  //std::cout<< Tc;
  //calculate Kalman gain K;
  MatrixXd K(5,3);
  K = Tc * S.inverse();
  //std::cout<< "K  \n"<<K;


  //update state mean and covariance matrix
  VectorXd z_diff;
  z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  double E;
  E = z_diff.transpose() * S.inverse() * z_diff;
  //std::cout<<"ER "<<E <<"\n";

  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();
  //std::cout<< "x_  \n"<<x_;
  //std::cout<< "P_  \n"<<P_;
}