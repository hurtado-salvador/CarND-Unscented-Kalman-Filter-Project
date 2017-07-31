#include "ukfback.h"
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
  // State Vector size 5
  n_x_ = 5;  //ok
  // Augmented State Vector size 7
  n_aug_ = 7; // ok
  //Define Lambda
  lambda_ = 3 - n_aug_; //ok
  // Time variable to calculate delta time.
  time_us_ = 0;//ok
  // Predicted X Sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1 );//ok
  Xsig_pred_.fill(0.0);
  //Weights calculation.
  weights_ = VectorXd(2 * n_aug_ + 1); //ok
  weights_(0) = (lambda_ / (lambda_ + n_aug_));
  for (int i = 1; i < 15; i++) {
    weights_(i) = (1 / (2 * (lambda_ + n_aug_)));

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  // initial state vector
  is_initialized_ = false;

  use_laser_ = use_radar_ = false; //what is this for??
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_) * 0.1; //Why *0.1??
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;  //Parameter to tune.
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; //Parameter to tune.
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
  std::cout << "Incoming measurement: " << meas_package.raw_measurements_
            << std::endl;
  // Calculate delta_t
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

     // Initialization
     if (!is_initialized_)
     {
       UKF::Inicializar(meas_package);
       return;;
     }


     UKF::Prediction(delta_t);

     if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
     {
        UKF::UpdateLidar(meas_package);
     }
     else if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
     {
         UKF::UpdateRadar(meas_package);
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

     //Generating augmented Sigma Points.
     MatrixXd SigmaAug = UKF::GenerateAugmentedSigmaPoints();

     // Predict Sigma Points.
     UKF::PredictSigmaPoints(SigmaAug, delta_t);

     //Predict x_ and P_
     UKF::PredictXandP();


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
     int n_z = 2;
     MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

     // mean predicted measurement
     VectorXd z_pred = VectorXd::Zero(n_z);

     // measurement covariance matrix S
     MatrixXd S = MatrixXd::Zero(n_z, n_z);

     /*******************************************************************************
      * Student part begin
      ******************************************************************************/

     for (int i = 0; i < 2 * n_aug_ + 1; i++) {
         VectorXd iter_col = Xsig_pred_.col(i);
         Zsig.col(i)[0] = iter_col[0];
         Zsig.col(i)[1] = iter_col[1];
         z_pred += weights_[i] * Zsig.col(i);
     }

     for (int i = 0; i < 2 * n_aug_ + 1; i++) {
         S += weights_[i] *
              ((Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose());
     }
     S(0, 0) += std_laspx_ * std_laspx_;
     S(1, 1) += std_laspy_ * std_laspy_;

     MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

     for (int i = 0; i < n_aug_ * 2 + 1; i++) {
         VectorXd Xdiff = (Xsig_pred_.col(i) - x_);
         while (Xdiff(3)> M_PI) Xdiff(3)-=2.*M_PI;
         while (Xdiff(3)<-M_PI) Xdiff(3)+=2.*M_PI;
         //NormAng(&x_diff[3]);
         Tc += weights_[i] * (Xdiff* (Zsig.col(i) - z_pred).transpose());
     }
     MatrixXd K = MatrixXd(n_x_, n_z);
     K = Tc * S.inverse();
  double E;
  VectorXd  restaX = meas_package.raw_measurements_ - z_pred;
  E = restaX.transpose() * S.inverse() * restaX;
  //std::cout<<"EL "<<E <<"\n";


     x_ = x_ + K * (meas_package.raw_measurements_ - z_pred);
     // calculate cross correlation matrix
     // calculate Kalman gain K;
     // update state mean and covariance matrix
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
  Zsig.fill(0.0);
     VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
     MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
     VectorXd z = VectorXd(n_z);
  z.fill(0.0);
     MatrixXd Tc = MatrixXd(n_x, n_z);
  Tc.fill(0.0);
  VectorXd resta = VectorXd(3);
  resta.fill(0.0);
  VectorXd Xdiff = VectorXd(5);
  Xdiff.fill(0.0);
  VectorXd Zdiff = VectorXd(3);
  Zdiff.fill(0.0);
  VectorXd z_diff;
  z_diff.fill(0.0);
  MatrixXd K(5,3);
  K.fill(0.0);

     z<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

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
     //calculate mean predicted measurement

     for(int i =0; i < 15; i++)
     {
         z_pred += weights_(i) * Zsig.col(i);
     }

     //calculate measurement covariance matrix S





     for(int i = 0; i < 15; i++)
     {
         resta = Zsig.col(i) - z_pred;
       while (resta[1] > M_PI) resta[1] -= 2.0 * M_PI;
       while (resta[1] < -M_PI) resta[1] += 2.0 * M_PI;

         S += weights_(i) * resta * resta.transpose();
     }
  S(0,0) += (std_radr_ * std_radr_);
  S(1,1) += (std_radphi_ * std_radphi_);
  S(2,2) += (std_radrd_ * std_radrd_);


     //calculate cross correlation matrix

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
     //calculate Kalman gain K;

     K = Tc * S.inverse();

     //update state mean and covariance matrix

     z_diff = z - z_pred;
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
     double E;
     E = z_diff.transpose() * S.inverse() * z_diff;
     //std::cout<<"ER "<<E <<"\n";

     x_ = x_ + K * (z_diff);

     P_ = P_ - K * S * K.transpose();

}

 void UKF::Inicializar(MeasurementPackage meas_package) {


     if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) //First measurement = LASER
     {
         x_[0] = meas_package.raw_measurements_[0];
         x_[1] = meas_package.raw_measurements_[1];

     } else if(meas_package.sensor_type_== MeasurementPackage::SensorType::RADAR) //First measurement = RADAR
     {
         double rho = meas_package.raw_measurements_[0];
         double phi = meas_package.raw_measurements_[1];
         x_[0] = rho * sin(phi);  // inverted sin and cos.
         x_[1] = rho * cos(phi);
     }
     is_initialized_ = true;
     return;
 }

 MatrixXd UKF::GenerateAugmentedSigmaPoints() {

     VectorXd x_aug = VectorXd(7);
   x_aug.fill(0.0);


     for (int i = 0; i < 5; i++) {
         x_aug[i] = x_[i];
     }
     x_aug[n_x_] = 0.0;
     x_aug[n_x_ + 1] = 0.0;


     MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
     P_aug.fill(0);

     P_aug.topLeftCorner(n_x_, n_x_) = P_;
     P_aug(n_x_, n_x_) = std_a_ * std_a_;
     P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
   P_aug *= lambda_ + n_aug_;

     MatrixXd SigmaAug(n_aug_, 2* n_aug_  + 1);
     SigmaAug.fill(0.0);


     MatrixXd A = MatrixXd(n_aug_, n_aug_);
     A.fill(0.0);
     A = P_aug.llt().matrixL();
     //A = A * sqrt(lambda_ + n_aug_);
     SigmaAug.col(0) = x_aug;

     for (int i = 0; i < n_aug_; i++) {
         SigmaAug.col(i + 1) = x_aug + A.col(i);
         SigmaAug.col(i + n_aug_+1) = x_aug -  A.col(i);
     }
     return SigmaAug;
 }
 void UKF::PredictSigmaPoints(MatrixXd SigmaAug, double delta_t) {
     Xsig_pred_.fill(0.0);
   for (int i = 0; i< 2*n_aug_+1; i++)
   {
     //extract values for better readability
     double p_x = SigmaAug(0,i);
     double p_y = SigmaAug(1,i);
     double v = SigmaAug(2,i);
     double yaw = SigmaAug(3,i);
     double yawd = SigmaAug(4,i);
     double nu_a = SigmaAug(5,i);
     double nu_yawdd = SigmaAug(6,i);

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
     Xsig_pred_(0,i) = px_p;
     Xsig_pred_(1,i) = py_p;
     Xsig_pred_(2,i) = v_p;
     Xsig_pred_(3,i) = yaw_p;
     Xsig_pred_(4,i) = yawd_p;
   }
     }


 void UKF::PredictXandP() {

     x_.fill(0.0);
   x_ = weights_[0] * Xsig_pred_.col(0);
         for (int i = 1; i < 15; i++) {
             x_ += weights_(i) * Xsig_pred_.col(i);
         }

     P_.setIdentity();
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

         // state difference
         VectorXd x_diff = Xsig_pred_.col(i) - x_;

         //angle normalization
         while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
         while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

         P_ = P_ + (weights_(i) *( x_diff * x_diff.transpose()));
     }
 }