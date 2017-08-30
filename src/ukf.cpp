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
    // x_ = VectorXd(5);

    // initial covariance matrix
    // P_ = MatrixXd(5, 5);

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

    // number of features in basic input vector
    n_x_ = x_.size();
    // 7 points for augmented vector
    n_aug_ = n_x_ + 2;
    // sigma point spreading factor
    lambda_ = 3 - n_aug_;

    // the number of sigma points
    n_sig = 2 * n_aug_ + 1;

    // weights to invert the spreading factor lambda to get the predicted covariance
    weights_ = VectorXd::Zero(n_sig);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for(int a = 0; a < n_sig; a++){ // 2n + 1
        weights_(a) = 0.5 / (lambda_ + n_aug_);
    }

    // initialize the state vector x
    x_ = VectorXd::Ones(5);

    // initialize state covariance matrix P
    P_ = MatrixXd::Identity(5, 5);
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

    if(!is_initialized_){
        // this is the first measurement
        if(meas_package.sensor_type_ == MeasurementPackage::LASER){
            // proceed with Laser data
            // set values in x_ state vector
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            // proceed with Radar data (rho, phi, rho_dot)
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            float dot = meas_package.raw_measurements_(2);

            // set values in x_ state vector
            x_(0) = rho * cos(phi);
            x_(1) = rho * sin(phi);
            x_(2) = dot * cos(phi);
            x_(3) = dot * sin(phi);

        }

        // initialized the x_ vector for the first measurement
        is_initialized_ = true;


    }
    else{
        // this is the 2nd or subsequent measurement

        /*****************************************************************************
        *                                Prediction                                  *
        *****************************************************************************/

        // calculate the change in time
        // assuming that delta_t is greater than 0 for each new measurement
        double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;// dt expressed in seconds
        // saving current timestamp as previous

        // calling predict function
        Prediction(delta_t);

        /*****************************************************************************
        *                                  Update                                    *
        *****************************************************************************/
        // if you want to turn off a sensor, set the use_laser_ or use_radar_ values to false in the constructor
        if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
            // updating for Laser
            UpdateLidar(meas_package);
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
            // updating for Radar
            UpdateRadar(meas_package);
        }
    }

    // setting the previous timestamp
    time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */


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
}
