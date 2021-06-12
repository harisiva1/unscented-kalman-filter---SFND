#include "ukf.h"
#include "Eigen/Dense"
#include"iostream"
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
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
  previous_timestamp_ = 0;
  is_initialized_ = false;
  n_x_ = 5;
  Xsig_ = MatrixXd(n_x_,2*n_x_+1);
  n_aug_ = 7;
  //lambda_ = 3-n_x_;
  lambda_aug = 3-n_aug_;
  x_aug = VectorXd(n_aug_);
  P_aug = MatrixXd(n_aug_, n_aug_);

  Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);

  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_aug+n_aug_)));
  weights_(0) = lambda_aug/(lambda_aug + n_aug_);


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {

      x_<<meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          0,
          0,
          0;
      P_<<std_laspx_*std_laspx_,0,0,0,0,
          0,std_laspy_*std_laspy_,0,0,0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0] ;
      double psi = meas_package.raw_measurements_[1];
      double rhod = meas_package.raw_measurements_[2];
      double vx = rhod*cos(psi);
      double vy = rhod*sin(psi);
      double v = sqrt(vx*vx + vy*vy);
      x_<<rho * cos(psi),
          rho * sin(psi),
          v,
          0,
          0;
      P_<<std_radr_*std_radr_,0,0,0,0,
          0,std_radr_*std_radr_,0,0,0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    }
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    std::cout<<"Init done"<<std::endl;
  }
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  std::cout<<"predict done"<<std::endl;

  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    std::cout<<"update lidar"<<std::endl;
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    std::cout<<"update radar"<<std::endl;
    UpdateRadar(meas_package);
  }
  
}

void UKF::Prediction(double delta_t) 
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  

  MatrixXd n = MatrixXd(2,2);
  n<<std_a_*std_a_,0,
     0,std_yawdd_*std_yawdd_;
  P_aug.fill(0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = n;
  

  // create square root matrix
  MatrixXd a = P_aug.llt().matrixL();
  // create augmented sigma points
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  Xsig_aug_.fill(0);
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) 
  {
    Xsig_aug_.col(i+1)     = x_aug + sqrt(lambda_aug+n_aug_) * a.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug+n_aug_) * a.col(i);
  }

  double t2  = delta_t*delta_t;

  // predict sigma points
  MatrixXd xk = MatrixXd(n_x_,2 * n_aug_ + 1);
  xk = Xsig_aug_.topLeftCorner(5,15);

  for(int i=0; i<(2*n_aug_+1);i++)
    {
        if(Xsig_aug_(4,i) !=0)
        {
            
            float a = (xk(2,i)/xk(4,i)) *(sin(xk(3,i) + (xk(4,i)*delta_t)) - sin(xk(3,i)));
            float b = 0.5*t2*cos(xk(3,i))*Xsig_aug_(5,i);
            float c = (xk(2,i)/xk(4,i)) *(-cos(xk(3,i) + (xk(4,i)*delta_t)) + cos(xk(3,i)));
            float d = 0.5*t2*sin(xk(3,i))*Xsig_aug_(5,i);
            Xsig_pred_(0,i) = xk(0,i) + a + b;
            Xsig_pred_(1,i) = xk(1,i) + c + d;
            Xsig_pred_(2,i) = xk(2,i) + delta_t * Xsig_aug_(5,i);
            Xsig_pred_(3,i) = xk(3,i) + (xk(4,i)*delta_t) + (0.5*t2*Xsig_aug_(6,i));
            Xsig_pred_(4,i) = xk(4,i) + delta_t * Xsig_aug_(6,i);

        }
        else
        {
            float a = xk(2,i) * cos(xk(3,i))*delta_t;
            float b = 0.5*t2*cos(xk(3,i))*Xsig_aug_(5,i);
            float c = xk(2,i) * sin(xk(3,i))*delta_t;
            float d = 0.5*t2*sin(xk(3,i))*Xsig_aug_(5,i);
            Xsig_pred_(0,i) = xk(0,i) + a + b;
            Xsig_pred_(1,i) = xk(1,i) + c + d;
            Xsig_pred_(2,i) = xk(2,i) + delta_t * Xsig_aug_(5,i);
            Xsig_pred_(3,i) = xk(3,i) + (xk(4,i)*delta_t) + (0.5*t2*Xsig_aug_(6,i));
            Xsig_pred_(4,i) = xk(4,i) + delta_t * Xsig_aug_(6,i);

        }    
    }


  // predict state mean
  x_.fill(0);
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
      x_(0) += weights_(i) * Xsig_pred_(0,i);
      x_(1) += weights_(i) * Xsig_pred_(1,i);
      x_(2) += weights_(i) * Xsig_pred_(2,i);
      x_(3) += weights_(i) * Xsig_pred_(3,i);
      x_(4) += weights_(i) * Xsig_pred_(4,i);
  }
  // predict state covariance matrix
  MatrixXd q = MatrixXd(5,1);
  P_.fill(0);
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
      q = Xsig_pred_.col(i) - x_;
      while (q(3)> M_PI) q(3)-=2.*M_PI;
      while (q(3)<-M_PI) q(3)+=2.*M_PI;
      P_ += weights_(i) * q * q.transpose();
  }
  std::cout<<"pd"<<std::endl; 
  //std::cout<<x_<<std::endl;
  //std::cout<<P_<<std::endl; 

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  n_z = 2;
  Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred = VectorXd(n_z);
  z_pred.fill(0);
  S = MatrixXd(n_z,n_z);
  S.fill(0);
  R = MatrixXd(n_z,n_z);
  R<<std_laspx_*std_laspx_,0,
     0,std_laspy_*std_laspy_;
     

  z = VectorXd(n_z);
  Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  kal = MatrixXd(n_x_, n_z);

  for(int i=0; i< (2 * n_aug_ + 1) ; i++)
  {

    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);

  }
  // calculate mean predicted measurement
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
    z_pred(0) += weights_(i) * Zsig(0,i);
    z_pred(1) += weights_(i) * Zsig(1,i);

  }

  
  // calculate innovation covariance matrix S
  MatrixXd q = MatrixXd(3,1);
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
    q = Zsig.col(i) - z_pred;
    while (q(1)> M_PI) q(1)-=2.*M_PI;
    while (q(1)<-M_PI) q(1)+=2.*M_PI;
    S += weights_(i) * q * q.transpose();
  }
  S += R;
  z<<meas_package.raw_measurements_[0],
     meas_package.raw_measurements_[1];


  // calculate cross correlation matrix
  for(int i=0 ; i < (2 * n_aug_ + 1) ; i++)
  {
    
   // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd kal = MatrixXd(n_x_, n_z);
  kal = Tc * S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + kal * z_diff;
  //x = x + kal*(z - z_pred);
  
  P_ = P_ - (kal*S*kal.transpose());

  std::cout<<"x and p updated from lidar"<<std::endl;
  std::cout<<x_<<std::endl;
  std::cout<<P_<<std::endl; 

}

void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  n_z = 3;
  Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred = VectorXd(n_z);
  z_pred.fill(0);
  S = MatrixXd(n_z,n_z);
  S.fill(0);
  R = MatrixXd(n_z,n_z);
  R<<std_radr_*std_radr_,0,0,
     0,std_radphi_*std_radphi_,0,
     0,0,std_radrd_*std_radrd_;

  z = VectorXd(n_z);
  Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  kal = MatrixXd(n_x_, n_z);

  for(int i=0; i< (2 * n_aug_ + 1) ; i++)
  {
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double psi = Xsig_pred_(3,i);
      double vel = Xsig_pred_(2,i);
      
      Zsig(0,i) = sqrt(px*px + py*py);
      Zsig(1,i) = atan2(py,px);
      Zsig(2,i) = ((px*cos(psi)*vel) + (py*sin(psi)*vel))/sqrt(px*px + py*py);
  }
  // calculate mean predicted measurement
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
      z_pred(0) += weights_(i) * Zsig(0,i);
      z_pred(1) += weights_(i) * Zsig(1,i);
      z_pred(2) += weights_(i) * Zsig(2,i);

  }

  
  // calculate innovation covariance matrix S
  MatrixXd q = MatrixXd(3,1);
  for(int i = 0 ; i < (2*n_aug_+1) ; i++)
  {
      q = Zsig.col(i) - z_pred;
      while (q(1)> M_PI) q(1)-=2.*M_PI;
      while (q(1)<-M_PI) q(1)+=2.*M_PI;
      S += weights_(i) * q * q.transpose();
  }
  S += R;
  z<<meas_package.raw_measurements_[0],
     meas_package.raw_measurements_[1],
     meas_package.raw_measurements_[2];

  // calculate cross correlation matrix
  for(int i=0 ; i < (2 * n_aug_ + 1) ; i++)
  {
    
   // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd kal = MatrixXd(n_x_, n_z);
  kal = Tc * S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + kal * z_diff;
  //x = x + kal*(z - z_pred);
  
  P_ = P_ - (kal*S*kal.transpose());

  std::cout<<"x and p updated from radar"<<std::endl;

  std::cout<<x_<<std::endl;
  std::cout<<P_<<std::endl; 
  
}
