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

  //Set the initiaization toggle variable
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
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

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set spreading parameter lambda
  lambda_ = 3 - n_aug_;
  
  //set weight vector
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);

  //Covariance matrix P
  /*P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;*/
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
	long long previous_timestamp;
	std::cout << "Starting a new cycle" << endl;

	if (!is_initialized_)
	{
		cout << "Very first initialization" << endl;
		x_ << 0, 0, 0, 0, 0;

		P_ << 1, 0, 0, 0, 0,
	  		0, 1, 0, 0, 0,
	  		0, 0, 1, 0, 0,
	  		0, 0, 0, 1, 0,
	  		0, 0, 0, 0, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/

			double x_pol_to_car, y_pol_to_car;
			double theta_ = meas_package.raw_measurements_[1];


			x_pol_to_car = meas_package.raw_measurements_[0] * cos(theta_);
			y_pol_to_car = meas_package.raw_measurements_[0] * sin(theta_);

			x_ << x_pol_to_car, y_pol_to_car, 0, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			/**
			Initialize state.
			*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

		}

		// done initializing, no need to predict or update
		previous_timestamp = meas_package.timestamp_;
		is_initialized_ = true;

		std::cout << "Initialization succesful" << endl;

		return;
	}
	std::cout << "Updating timestamp" << endl;
	
	double dt = (meas_package.timestamp_ - previous_timestamp) / 1000000.0;
	previous_timestamp = meas_package.timestamp_;

	//Call predict function and predict state vector x_ and covariance matrix P_
	std::cout << "Calling prediction" << endl;
	Prediction(dt);
	std::cout << "Prediction completed" << endl;

	//Calls either radar or lidar update functions
	

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		std::cout << "Calling radar update" << endl;
		UpdateRadar(meas_package);
	}
	

	else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		std::cout << "Calling laser update" << endl;

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
 

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);
	x_aug.fill(0.0);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);
	P_aug.fill(0.0);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.fill(0.0);

	//Defining a 2x2 Q matrix tp augment covariance P_aug
	MatrixXd Q =  MatrixXd(2,2);
	Q.fill(0.0);

	Q << std_a_ * std_a_, 0,
		0, std_yawdd_ * std_yawdd_;

	//Adding the first 5 known values of x to x_aug
	x_aug.head(n_x_) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//Augmenting P_aug by adding known 5x5 P values plus Q 2x2 matrix defined before
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	//P_aug.bottomLeftCorner(2, 2) = Q;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	std::cout << "x and P augmented succesfull" << endl;

	//Create a matrix to store squared P_aug
	MatrixXd L = MatrixXd(7, 7);
	L.fill(0.0);

	L = P_aug.llt().matrixL();

	/*Eigen::LLT<MatrixXd> lltOfPaug(P_aug);

	if (lltOfPaug.info() == Eigen::NumericalIssue) {

		cout << "LLT failed!" << endl;
		throw range_error("LLT failed");
	}

	MatrixXd P_sqr = lltOfPaug.matrixL();*/


	//Start sigma point creation
	Xsig_aug.col(0) = x_aug;

	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	std::cout << "Sigma points generation completed" << endl;

	//create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0.0);

	//write predicted sigma points into right column
	MatrixXd x_k = MatrixXd(5, 15);
	x_k.fill(0.0);

	x_k.topLeftCorner(5, 15) = Xsig_aug.topLeftCorner(5, 15);

	//Matrix for deterministic part
	MatrixXd det_vector = MatrixXd(5, 1);

	//Matrix for stochastic part
	MatrixXd stoch_vector = MatrixXd(5, 1);
	

	for (int i = 0; i < 15; i++)
	{
		//Get values into variable for readability
		/*double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		if (!(fabs(yawd) > 0.001))
		{
			//Calculating deterministic part when yawd equals 0
			det_vector(0, 0) = v * cos(yaw) * delta_t;
			det_vector(1, 0) = v * sin(yaw) * delta_t;
			det_vector(2, 0) = 0;
			det_vector(3, 0) = yawd * delta_t;
			det_vector(4, 0) = 0;
		}

		else
		{
			//Calculating deterministic part when yawd equals 1
			det_vector(0, 0) = (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
			det_vector(1, 0) = (v / yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw));
			det_vector(2, 0) = 0;
			det_vector(3, 0) = yawd * delta_t;
			det_vector(4, 0) = 0;
		}

		//Calculating stochastic part
		stoch_vector(0, 0) = (0.5) * (delta_t * delta_t) * cos(yaw) * nu_a;
		stoch_vector(1, 0) = (0.5) * (delta_t * delta_t) * sin(yaw) * nu_a;
		stoch_vector(2, 0) = delta_t * nu_a;
		stoch_vector(3, 0) = (0.5) * (delta_t * delta_t) * nu_yawdd;
		stoch_vector(4, 0) = delta_t * nu_yawdd;

		//Building predicted sigma values by putting deterministic and stochastic calculationa all together
		Xsig_pred_.col(i) = x_k.col(i) + det_vector + stoch_vector;*/

		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
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
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
	std::cout << "Sigma values predicted" << endl;
	std::cout << Xsig_pred_ << endl;

	//Setting weights
	/*for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		if (i + 1 == 1)
		{
			weights_(i) = lambda_ / (lambda_ + n_aug_);
		}
		else
		{
			weights_(i) = 0.5 * 1 / (lambda_ + n_aug_);
		}

	}*/

	weights_(0) = lambda_ / (lambda_ + n_aug_);
		for (int i = 1; i < 2 * n_aug_ + 1; ++i)
		{
			weights_(i) = 0.5/(lambda_ + n_aug_);
		}


	//Predicting x_
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	std::cout << "x_ predicted" << endl;
	std::cout << x_ << endl;

	//Predicting P_
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		std::cout << "x_diff calculated" << endl;
		std::cout << x_diff << endl;

		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
		//x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}

	std::cout << "P_ predicted" << endl;
	std::cout << P_ << endl;
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

	//Setting matrix H for lidar updates
	MatrixXd H = MatrixXd(2, 5);
	H.fill(0.0);

	H << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;

	std::cout << "H matrix generated"<<endl;

	//Noise matrix R
	MatrixXd R = MatrixXd(2, 2);
	R.fill(0.0);

	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_ * std_laspy_;

	std::cout << "R matrix generated" << endl;

	//Setting up the variable requirements
	VectorXd z_pred = H * x_;
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	std::cout << "Variable requirements done" <<endl;

	// Calculating the new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H) * P_;

	std::cout << "New estimation done" << endl;

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

    //create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
	Zsig.fill(0.0);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(3);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(3, 3);
	S.fill(0.0);

	//Calculate sigma points in measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//Extract parameters for code readability
		float x_pred = (Xsig_pred_(0, i));
		float y_pred = (Xsig_pred_(1, i));
		float v_pred = (Xsig_pred_(2, i));
		float yaw_pred = (Xsig_pred_(3, i));

		//Calculate polar coordinate parameters
		float rho = sqrt((x_pred * x_pred) + (y_pred * y_pred));
		float phi = atan2(y_pred, x_pred);
		float rho_dot = ((x_pred * cos(yaw_pred) * v_pred) + (y_pred * sin(yaw_pred) * v_pred)) / rho;

		//Assign calcualated values to corresponding column
		Zsig(0, i) = rho;
		Zsig(1, i) = phi;
		Zsig(2, i) = rho_dot;
	}
	std::cout << "Sigma points in measure space calculated" << endl;
	std::cout << Zsig << endl;

	//Calculate mean z measurements
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}
	std::cout << "Mean z finished" << endl;
	std::cout << z_pred << endl;

	//Calculate matrix measurement noise R
	MatrixXd R = MatrixXd(3, 3);
	R.fill(0.0);

	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_ * std_radphi_, 0,
		0, 0, std_radrd_ * std_radrd_;

	//Calculate covariance measurement matrix S

	VectorXd z_diff = VectorXd(3);
	z_diff.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_diff = Zsig.col(i) - z_pred;

		while (z_diff(1) > M_PI) z_diff(1) -= 2. *M_PI;
		while (z_diff(1) < -M_PI) z_diff(1) += 2. *M_PI;
		//z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	//Add noise to S
	S = S + R;
	std::cout << "Covariance S done" << endl;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, 3);
	Tc.fill(0.0);

	//Create vector x_diff to calculate difference between pred and measurement
	VectorXd x_diff = VectorXd(5);
	x_diff.fill(0.0);

	//Calculate cross correlation matrix Tc
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_diff = Xsig_pred_.col(i) - x_;
		z_diff = Zsig.col(i) - z_pred;
		std::cout << "x_ diff is:" << endl << x_diff << endl;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
		//z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));
		std::cout << "Angle of z normalized" << endl;

		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
		//x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));

		std::cout << "Angle of x normalized" << endl;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
		std::cout << " cycle of Tc done" << endl;
	}
	std::cout << "Cross correlation Tc done" << endl;

	//Declare and calculate kalman gain K_g
	MatrixXd K_g = MatrixXd::Zero(3,3);
	K_g = Tc * S.inverse();

	//Get measured values for readability
	VectorXd z_meas = VectorXd(3);
	z_meas(0) = meas_package.raw_measurements_(0);
	z_meas(1) = meas_package.raw_measurements_(1);
	z_meas(2) = meas_package.raw_measurements_(2);

	//Update state vector x_ and covariance matrix P_
	x_ = x_ + K_g * (z_meas - z_pred);
	P_ = P_ - K_g * S * K_g.transpose();
}
