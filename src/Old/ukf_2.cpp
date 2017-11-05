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

using namespace std;
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5;

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

  //long long previous_timestamp_;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_= 5;
  n_aug_= 7;
  is_initialized_= false;
  x_=VectorXd(n_x_);
  x_<< 1,1,0,0,0;
  P_ << 1,0,0,0,0,
		  0,1,0,0,0,
		  0,0,1,0,0,
		  0,0,0,1,0,
		  0,0,0,0,1;

  lambda_= 3-n_aug_;
  Xsig_pred_=MatrixXd(n_x_,2*n_aug_+1);
  weights_=VectorXd(2*n_aug_+1);
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
	long long previous_timestamp_;
	
	if (!is_initialized_){
		if (meas_package.sensor_type_== MeasurementPackage::LASER){
			x_<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;
			previous_timestamp_= meas_package.timestamp_;
		}
		if (meas_package.sensor_type_== MeasurementPackage::RADAR){
			float rho= meas_package.raw_measurements_[0];
			float theta= meas_package.raw_measurements_[1];
			x_ << (rho*std::cos(theta)), (rho*std::sin(theta)), 0,0, 0;
			previous_timestamp_= meas_package.timestamp_;
		}
		is_initialized_=true;

		return;
	}
	std::cout << "start prediction" << std::endl;



//	use_laser_=false;
	if (meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_){
		double dt = (meas_package.timestamp_-previous_timestamp_)/1000000.0;
		previous_timestamp_= meas_package.timestamp_;
		Prediction(dt);
		std::cout << "end prediction" << std::endl;
		UpdateLidar(meas_package);
	}
	else if(meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_) {
		double dt = (meas_package.timestamp_-previous_timestamp_)/1000000.0;
		previous_timestamp_= meas_package.timestamp_;
		Prediction(dt);
		std::cout << "end prediction" << std::endl;
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

//	lambda_= 3-n_x_;
	VectorXd x_aug= VectorXd(n_aug_);
	MatrixXd P_aug= MatrixXd(n_aug_,n_aug_);
	lambda_= 3-n_aug_;
	MatrixXd Xsig_aug= MatrixXd(n_aug_, 2*n_aug_+1);
	x_aug.head(5)= x_;
	x_aug(5)=0;
	x_aug(6)=0;
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5)=P_;
	P_aug(5,5)= pow(std_a_,2);
	P_aug(6,6)= pow(std_yawdd_,2);
	MatrixXd A_aug= P_aug.llt().matrixL();


	Xsig_aug.col(0)= x_aug;
	for (int i=0; i<n_aug_; i++){
		Xsig_aug.col(i+1)= x_aug+sqrt(lambda_+n_aug_)*A_aug.col(i);
		Xsig_aug.col(i+1+n_aug_)= x_aug- sqrt(lambda_+ n_aug_)*A_aug.col(i);

		}
	cout<< "sig moid x  "<<Xsig_aug<<endl;


    Xsig_pred_.fill(0.0);
	for(int i=0; i<2*n_aug_+1; i++){
		double p_x= Xsig_aug(0,i);
		double p_y= Xsig_aug(1,i);
		double v= Xsig_aug(2,i);
		double yaw= Xsig_aug(3,i);
		double yawd= Xsig_aug(4,i);
		double nu_a= Xsig_aug(5,i);
		double nu_yawd=Xsig_aug(6,i);

		double p_px, p_py, p_v, p_yaw, p_yawd;

		p_yaw= yaw+ yawd*delta_t;
		p_yawd= yawd;
		if (fabs(yawd)>.001){
			p_px= p_x+ (v/yawd)*(sin(yaw+yawd*delta_t)-sin(yaw));
			p_py= p_y+ (v/yawd)*(cos(yaw)-cos(yaw+yawd*delta_t));
		}
		else{
			p_px= p_x+ v*cos(yaw)*delta_t;
			p_py= p_y+ v*sin(yaw)*delta_t;
		}
		p_px += .5*pow(delta_t,2)*cos(yaw)*nu_a;
		p_py += .5*pow(delta_t,2)*sin(yaw)*nu_a;

		p_v=v+ delta_t*nu_a;
		p_yaw= p_yaw+ .5*pow(delta_t,2)*nu_yawd;
		p_yawd= p_yawd+0+delta_t*nu_yawd;

		Xsig_pred_(0,i)=p_px;
		Xsig_pred_(1,i)= p_py;
		Xsig_pred_(2,i)= p_v;
		Xsig_pred_(3,i)= p_yaw;
		Xsig_pred_(4,i)= p_yawd;




	}
	cout<< "sig pred x \n "<<Xsig_pred_<<endl;
	lambda_ = 3 - n_aug_;

	weights_(0)= lambda_/(lambda_+n_aug_);
	for (int i=1; i<2*n_aug_+1; i++)	{
		weights_(i)= .5/(lambda_+n_aug_);
	}
	cout<< "weights  "<<weights_<<endl;

	x_.fill(0.0);
	for (int i=0; i<2*n_aug_+1; i++){

		x_ += weights_(i)*Xsig_pred_.col(i);
//		cout<< " x=  \n"<<x_<<endl;
	}

	cout<< "final predicted x=  \n"<<x_<<endl;



	P_.fill(0.0);
	 for (int i=0; i<2*n_aug_+1; i++){
//		 cout<<"X sig= "<<Xsig_pred_.col(i)<<endl;
//		 cout<<"state= "<< x_<<endl;

		  VectorXd x_diff= Xsig_pred_.col(i)-x_;
//		  cout<< "diff x  "<<x_<<endl;
		  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
//		  cout<< "diff x2  "<<x_<<endl;
		  P_+=  weights_(i)* x_diff*x_diff.transpose();
//		  cout<< "P  "<<P_<<endl;
//		  cout << i<< endl;
	  }
	 cout<< "pred P  "<<P_<<endl;
	 cout<< "end of predection  "<<endl;

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
	cout<<"start update lidar"<<endl;
	int n_z=2;
	VectorXd z_pred= VectorXd(n_z);
	VectorXd z= VectorXd(2);
	z<< meas_package.raw_measurements_(0),
			meas_package.raw_measurements_(1);

	MatrixXd Zsig= MatrixXd(n_z, 2*n_aug_+1);
	Zsig.fill(0.0);
	z_pred.fill(0.0);

	for (int i=0; i<2*n_aug_+1; i++){
		double p_x= Xsig_pred_(0,i);
		double p_y= Xsig_pred_(1,i);


		Zsig(0,i)= p_x;
		Zsig(1,i)= p_y;
		z_pred+= weights_(i)* Zsig.col(i);

	}



	MatrixXd R=MatrixXd(n_z,n_z);
	R.fill(0.0);
	R(0,0)=pow(std_laspx_,2);
	R(1,1)=pow(std_laspy_,2);


	MatrixXd S= MatrixXd(n_z,n_z);
	S.fill(0.0);

	for (int i=0; i<2*n_aug_+1;i++){
		VectorXd diff_z= Zsig.col(i)- z_pred;
		S+=weights_(i)* diff_z*diff_z.transpose();
	}
	S+= R;

	MatrixXd Tc= MatrixXd(n_x_, n_z);
	Tc.fill(0.0);


	for (int i=0; i<2*n_aug_+1; i++){

		VectorXd z_diff= Zsig.col(i)-z_pred;
		VectorXd x_diff= Xsig_pred_.col(i)-x_;
//		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
//		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		Tc+= weights_(i)*x_diff*z_diff.transpose();
	}

	MatrixXd K= Tc*S.inverse();
	x_= x_+ K*(z-z_pred);
	P_= P_- K*S*K.transpose();
	std::cout << "Updated state x: " << std::endl << x_ << std::endl;



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
	cout<<"start update radar"<<endl;
	int n_z=3;
	VectorXd z_pred= VectorXd(n_z);
	VectorXd z=  meas_package.raw_measurements_;
	cout<<"z= "<<z<<endl;

	MatrixXd Zsig= MatrixXd(n_z, 2*n_aug_+1);

	for (int i=0; i<2*n_aug_+1; i++){
		double p_x= Xsig_pred_(0,i);
		double p_y= Xsig_pred_(1,i);
		double v= Xsig_pred_(2,i);
		double yaw=Xsig_pred_(3,i);

		Zsig(0,i)= sqrt(pow(p_x,2)+pow(p_y,2));
		Zsig(1,i)= atan2(p_y,p_x);
		double vx= cos(yaw)*v;
		double vy= sin(yaw)*v;
		Zsig(2,i)= (p_x*vx+p_y*vy)/sqrt(pow(p_x,2)+pow(p_y,2));


	}
	cout<<"Zsig= "<<Zsig <<endl;
	z_pred.fill(0.0);
	for (int i=0; i<2*n_aug_+1;i++){
		z_pred+= weights_(i)* Zsig.col(i);
	}
	cout<<"z_pred= "<<z_pred <<endl;

	MatrixXd R=MatrixXd::Zero(n_z,n_z);
//	R.fill(0.0);
	R(0,0)=pow(std_radr_,2);
	R(1,1)=pow(std_radphi_,2);
	R(2,2)= pow(std_radrd_,2);

	MatrixXd S= MatrixXd::Zero(n_z,n_z);
//	S.fill(0.0);

	for (int i=0; i<2*n_aug_+1;i++){
		VectorXd diff_z= Zsig.col(i)- z_pred;
		while (diff_z(1)> M_PI) diff_z(1)-=2.*M_PI;
		while (diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;
		S+=weights_(i)* diff_z*diff_z.transpose();
	}
	S+= R;
	cout<<"S= "<<S <<endl;

	MatrixXd Tc= MatrixXd(n_x_, n_z);
	Tc.fill(0.0);


	for (int i=0; i<2*n_aug_+1; i++){

		VectorXd z_diff= Zsig.col(i)-z_pred;
	    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
	    VectorXd x_diff= Xsig_pred_.col(i)-x_;
	    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
	    Tc+= weights_(i)*x_diff*z_diff.transpose();
	}
	cout<<"Tc= "<<Tc <<endl;

	MatrixXd K= Tc*S.inverse();
	VectorXd z_diff = z - z_pred;
	 while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	 while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
	x_= x_+ K*z_diff;
	P_= P_- K*S*K.transpose();

	std::cout << "Updated state x: " << std::endl << x_ << std::endl;







}
