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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
     P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
   n_x_ = 5;
   n_aug_ = 7;
   lambda_ = 3 - n_aug_;
   time_us_ = 0.0;

   Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1 );
   weights_ = VectorXd(2 * n_aug_ + 1);
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
     int St = meas_package.sensor_type_;
     // Initialization
     if (!is_initialized_)
     {
         UKF::Inicializar(St, meas_package);
     }
     // Calculate delta_t
     double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

     time_us_ = meas_package.timestamp_;



     if(St == 0)
     {
        UKF::UpdateLidar(meas_package);
     }
     else if(St == 1)
     {
         UKF::Prediction(delta_t);
         UKF::UpdateRadar(meas_package);
     }
    // else{std::cout<< "\n No detecta RADAR o LIDAR!!!!!!!!!!\n";}

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

     //std::cout<<"X_ predicted \n"<< x_<<"\n";
     //std::cout<<"P_ Predicted \n"<< P_<<"\n";

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
     std::cout<<"E "<<E <<"\n";

     x_ = x_ + K * (z_diff);

     P_ = P_ - K * S * K.transpose();
     //std::cout<< "x_  \n"<<x_;
     //std::cout<< "P_  \n"<<P_;



}

 void UKF::Inicializar(int St, MeasurementPackage meas_package) {


     if (St == 0) //First measurement = LASER
     {
         double px = meas_package.raw_measurements_[0];
         double py = meas_package.raw_measurements_[1];
         x_ << px, py, 0, 0, 0;
         P_.fill(0.05);
         time_us_ = meas_package.timestamp_;

     } else if (St == 1) //First measurement = RADAR
     {
         double rho = meas_package.raw_measurements_[0];
         double phi = meas_package.raw_measurements_[1];
         double px = rho * cos(phi);
         double py = rho * sin(phi);
         x_ << px, py, 0, 0, 0;
         P_.fill(0.05);
         time_us_ = meas_package.timestamp_;
     }
     is_initialized_ = true;
     return;
 }

 MatrixXd UKF::GenerateAugmentedSigmaPoints() {

     VectorXd x_aug(7);
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

     //std::cout<<"\n P_aug \n"<<P_aug<<"\n";

     MatrixXd SigmaAug(n_aug_, 2* n_aug_  + 1);
     SigmaAug.fill(0.0);
     MatrixXd A = P_aug.llt().matrixL();
     A = A * sqrt(lambda_ + n_aug_);
     SigmaAug.col(0) = x_aug;
     VectorXd tem1(7);
     VectorXd tem2(7);

     for (int i = 0; i < n_aug_; i++) {
         tem1 = x_aug + A.col(i);
         tem2 = x_aug - A.col(i);
         SigmaAug.col(i + 1) = tem1;
         SigmaAug.col(i + n_aug_+1) = tem2;
     }

     //std::cout << "Sigma Augmented  \n" << SigmaAug;
     return SigmaAug;
 }
 void UKF::PredictSigmaPoints(MatrixXd SigmaAug, double delta_t) {
     Xsig_pred_.fill(0.0);
     VectorXd coln = VectorXd(5);
     for (int i = 0; i < 15; i++) {
         double Vk = SigmaAug(2, i);
         double Yk = SigmaAug(3, i);
         double YDk = SigmaAug(4, i);
         double Vak = SigmaAug(5, i);
         double YDDk = SigmaAug(6, i);
         double Pxp;
         double Pyp;
         double vp;
         double Yp;
         double YDp;


         if (fabs(YDk) > 0.000001) {
             Pxp = SigmaAug(0, i) + ((Vk / YDk) * (sin(Yk + (YDk * delta_t)) - sin(Yk))) +
                   (0.5 * (delta_t * delta_t) * cos(Yk) * Vak);
             Pyp = SigmaAug(1, i) + (Vk / YDk) * (cos(Yk) - cos(Yk + (YDk * delta_t))) +
                   (0.5 * (delta_t * delta_t) * sin(Yk) * Vak);
             vp = SigmaAug(2, i) + delta_t * Vak;
             Yp = SigmaAug(3, i) + YDk * delta_t + (0.5 * (delta_t * delta_t)) * YDDk;
             YDp = SigmaAug(4, i) + delta_t * YDDk;

         } else {
             Pxp = SigmaAug(0, i) + (Vk * cos(Yk) * delta_t) + (0.5 * (delta_t * delta_t) * cos(Yk) * Vak);
             Pyp = SigmaAug(1, i) + (Vk * (sin(Yk) * delta_t)) + (0.5 * (delta_t * delta_t) * sin(Yk) * Vak);
             vp = SigmaAug(2, i) + delta_t * Vak;
             Yp = SigmaAug(3, i) + (0.5 * (delta_t * delta_t)) * YDDk;
             YDp = SigmaAug(4, i) + delta_t * YDDk;
         }

         coln << Pxp, Pyp, vp, Yp, YDp;

         Xsig_pred_.col(i) = coln;
     }


     //std::cout << "Sigma Predicted \n" << Xsig_pred_ << "\n";

 }

 void UKF::PredictXandP() {
     weights_(0) = (lambda_ / (lambda_ + n_aug_));
     double weight2 = (1 / (2 * (lambda_ + n_aug_)));
     for (int i = 1; i < 15; i++) {
         weights_(i) = weight2;
     }
     //std::cout << weights_ << "\n";
     x_.fill(0.0);
     for (int i = 0; i < 5; i++) {
        // x_(i) = 0.0;
         for (int j = 0; j < 15; j++) {
             x_(i) += weights_(j) * Xsig_pred_(i, j);
         }
         std::cout << "x_ predicted \n" << x_ << "\n";
     }
     P_.setIdentity();
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

         // state difference
         VectorXd x_diff = Xsig_pred_.col(i) - x_;
         std::cout <<"i "<<i<< " Xdiff  "<<x_diff<<"Xsig_pred_ "<<Xsig_pred_<<" -x "<< x_<<"\n";
         //std::cout << "angle norm before \n" << x_diff << " x_diff \n";
         //angle normalization
         //std::cout << "angle norm before \n" << x_diff(3) << " x_diff \n";
         while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
         //std::cout << "angle norm middle \n" << x_diff(3) << " x_diff \n";
         while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
         //std::cout << "angle norm end \n" << x_diff(3) << " x_diff \n";

         P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
         std::cout<<"P_ "<<P_<<"\n";
     }
 }