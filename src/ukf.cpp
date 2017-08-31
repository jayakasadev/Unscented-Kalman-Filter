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

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, n_sig);

    // number of radar input values
    n_radar = 3;
    // number of laser input values
    n_laser = 2;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if(!is_initialized_){
        // this is the first measurement
        if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
            // proceed with Laser data
            // set values in x_ state vector
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);

        }
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
            // proceed with Radar data (rho, phi, rho_dot)
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double dot = meas_package.raw_measurements_(2);

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
    // generating augmented sigma points
    // matrix for generated sigma points
    MatrixXd Xsig_gen_ = MatrixXd(n_aug_, n_sig);
    // generating augmented sigma points
    // returns a 7 dimentional vector
    GenerateAugmentedSigma(Xsig_gen_);

    // predict sigma points by applying process model
    PredictSigmaPoints(Xsig_gen_, delta_t);

    // predict the mean and covariance
    CalculateMeanCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // measurement generated sigma points
    MatrixXd Zsig = MatrixXd(n_laser, n_sig);

    // predicted measurement mean
    VectorXd zpred = VectorXd::Zero(n_laser);

    // measurement covariance matrix
    MatrixXd S = MatrixXd::Zero(n_laser, n_laser);

    // cross correllation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_laser);

    // predict laser measurements
    PredictLaserMeasurements(Zsig, zpred, S, Tc);

    //update state
    NIS_laser_ = UpdateState(meas_package, Zsig, zpred, S, Tc);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    // measurement generated sigma points
    MatrixXd Zsig = MatrixXd(n_radar, n_sig);

    // predicted measurement mean
    VectorXd zpred = VectorXd::Zero(n_radar);

    // measurement covariance matrix
    MatrixXd S = MatrixXd::Zero(n_radar, n_radar);

    // cross correllation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_radar);

    // predict radar measurements
    PredictRadarMeasurements(Zsig, zpred, S, Tc);

    //update state
    NIS_radar_ = UpdateState(meas_package, Zsig, zpred, S, Tc);
}

/**
 * Private Methods
 */

/**
 * Method for generating the Augmented Sigma Points
 *
 * @param Xsig_gen_
 */
void UKF::GenerateAugmentedSigma(MatrixXd &Xsig_gen_) {
    //generate augmented state vector
    VectorXd x_aug = VectorXd::Zero(n_aug_);
    x_aug.head(n_x_) = x_; // TODO possible source of bug

    // generate augmented state covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = pow(std_a_, 2);
    P_aug(6, 6) = pow(std_yawdd_, 2);

    // retrieves A matrix of A x A^T
    MatrixXd A = P_aug.llt().matrixL();


    // set the columns of sigma points matrix
    Xsig_gen_.col(0) = x_aug;
    // get the remaining sigma points
    for(int a = 0; a < n_sig; a++){
        // common value
        VectorXd var = sqrt(lambda_ - n_aug_) * A.col(a);
        Xsig_gen_.col(a + 1) = x_ + var;
        Xsig_gen_.col(a + 1 + n_x_) = x_ - var;
    }
}

/**
 * Method for predicting sigma points by applying process model to generated sigma points
 * We are predicting where the sigma points will be at time k+1 given the generated sigma
 * points at time k
 *
 * Takes a 7x15 matrix in and outputs a 5x15 matrix in Xsig_pred_
 *
 * @param Xsig_gen_
 * @param delta_t
 */
void UKF::PredictSigmaPoints(MatrixXd &Xsig_gen_, const double &delta_t) {
    // could re/initialize or empty Xsig_pred_ here
    // but it is better to just overwrite the values

    //predict sigma values
    for(int a = 0; a < n_sig; a++){
        // extract all 7 features from generated sigma points
        double p_x = Xsig_gen_(0, a);
        double p_y = Xsig_gen_(1, a);
        double v = Xsig_gen_(2, a);
        double yaw = Xsig_gen_(3, a);
        double yaw_d = Xsig_gen_(4, a);
        // process noise values
        double nu_a = Xsig_gen_(5, a);
        double nu_yaw_dd = Xsig_gen_(6, a);

        // predicted state(x and y), speed, yaw angle, yaw rate
        double px_p, py_p, v_p = v, yaw_p = yaw + yaw_d * delta_t, yaw_d_p = yaw_d;

        //avoid division by zero
        if(fabs(yaw_d) > 0.001){ // TODO consider making this a global for more control
            // turning

            //shared value
            double inner = yaw + yaw_d * delta_t;

            px_p = p_x + v / yaw_d * (sin(inner) - sin(yaw));
            py_p = p_y + v / yaw_d * (cos(yaw) - cos(inner));
        }
        else{
            // going straight
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        // adding noise
        double shared = 0.5 * pow(delta_t, 2);

        px_p += shared * nu_a * cos(yaw);
        py_p += shared * nu_a * sin(yaw);

        v_p += nu_a * delta_t;

        yaw_p += shared * nu_yaw_dd;

        yaw_d_p += nu_yaw_dd * delta_t;

        //save predicted sigma points
        Xsig_pred_(0, a) = px_p;
        Xsig_pred_(1, a) = py_p;
        Xsig_pred_(2, a) = v_p;
        Xsig_pred_(3, a) = yaw_p;
        Xsig_pred_(4, a) = yaw_d_p;
    }
}

/**
 * Method for calculating the predicted state mean and covariance
 */
void UKF::CalculateMeanCovariance(){
    // resetting the values in the state vector and the covariance vector because we no longer need those values
    x_.fill(0);
    P_.fill(0);

    // predicted state mean
    for(int a = 0; a < n_sig; a++){
        x_ += weights_(a) * Xsig_pred_.col(a);
    }

    // predicted state covariance
    for(int a = 0; a < n_sig; a++){
        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;

        // angle normalization
        while(x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while(x_diff(3) < M_PI) x_diff(3) += 2. * M_PI;

        P_ += weights_(a) * x_diff * x_diff.transpose();
    }
}

/**
 * Method for transforming predicted state into measurement state for Laser and calculating cross correlation matrix
 *
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 */
void UKF::PredictLaserMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // transform sigma points into measurement space
    for(int a = 0; a < n_sig; a++){
        // extract values
        Zsig(0, a) = Xsig_pred_(0, a);
        Zsig(1, a) = Xsig_pred_(1, a);
    }

    // mean predicted measurement
    for(int a = 0; a < n_sig; a++){
        zpred += weights_(a) * Zsig.col(a);
    }

    // measurement covariance matrix
    for(int a = 0; a < n_sig; a++){
        // residual
        VectorXd z_diff = Zsig.col(a) - zpred;

        S += weights_(a) * z_diff * z_diff.transpose();

        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;

        // calculate cross correlation matrix
        Tc = Tc + weights_(a) * x_diff * z_diff.transpose();
    }

    MatrixXd R = MatrixXd::Zero(n_laser, n_laser);
    R(0, 0) = pow(std_laspx_, 2);
    R(1, 1) = pow(std_laspy_, 2);

    S += R;
}

/**
 * Method for transforming predicted state into measurement state for Radar and calculating cross correlation matrix
 *
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 */
void UKF::PredictRadarMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // transform sigma points into measurement space
    for(int a = 0; a < n_sig; a++){
        // extract values
        double p_x = Xsig_pred_(0, a);
        double p_y = Xsig_pred_(1, a);
        double v = Xsig_pred_(2, a);
        double yaw = Xsig_pred_(3, a);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        double shared = sqrt(pow(p_x, 2) + pow(p_y, 2));

        // measurement model
        Zsig(0, a) = shared;
        Zsig(1, a) = atan2(p_y, p_x);
        Zsig(2, a) = (p_x * v1 + p_y * v2) / shared;
    }

    // mean predicted measurement
    for(int a = 0; a < n_sig; a++){
        zpred += weights_(a) * Zsig.col(a);
    }

    // measurement covariance matrix
    for(int a = 0; a < n_sig; a++){
        // residual
        VectorXd z_diff = Zsig.col(a) - zpred;

        // angle normalization
        while(z_diff(1) < M_PI) z_diff(1) -= 2. * M_PI;
        while(z_diff(1) > M_PI) z_diff(1) += 2. * M_PI;

        S += weights_(a) * z_diff * z_diff.transpose();

        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;

        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        // calculate cross correlation matrix
        Tc = Tc + weights_(a) * x_diff * z_diff.transpose();
    }

    // add measurement noise
    MatrixXd R = MatrixXd::Zero(n_radar, n_radar);
    R(0, 0) = pow(std_radr_, 2);
    R(1, 1) = pow(std_radphi_, 2);
    R(2, 2) = pow(std_radr_, 2);

    S += R;
}

/**
 * Calculates Kalman Gain and NIS
 *
 * @param meas_package
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 * @return NIS value as double
 */
double UKF::UpdateState(MeasurementPackage meas_package, MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // Kalman Gain
    MatrixXd K = Tc * S.inverse();

    VectorXd z_diff = meas_package.raw_measurements_ - zpred;

    // update state mean and covariance
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

     // calculate NIS
    return z_diff.transpose() * S.inverse() * z_diff;
}