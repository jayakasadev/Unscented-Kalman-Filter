#include "ukf.h"
#include "Eigen/Dense"

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

    // initialize the state vector x
    x_ = VectorXd::Zero(5);

    // initialize state covariance matrix P
    P_ = MatrixXd::Identity(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1; // TODO tune

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3; // TODO tune

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
    for(int a = 1; a < n_sig; a++){ // 2n + 1
        weights_(a) = 0.5 / (lambda_ + n_aug_);
    }

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, n_sig);

    // number of radar input values
    n_radar = 3;
    // number of laser input values
    n_laser = 2;

    // just doing this defensively, since in c++ an uninitialized variable just takes on whatever was there previously
    is_initialized_ = false;

    // change this as needed
    // unscented_ = true;
    unscented_ = false;

    // for use with linear kalman filter for laser update
    H_ = MatrixXd::Zero(2, 5);
    H_(0, 0) = 1;
    H_(1, 1) = 1;

    // redudant to recreate every time
    // measurement noise covariance for laser
    R_laser_ = MatrixXd::Zero(n_laser, n_laser);
    R_laser_(0, 0) = pow(std_laspx_, 2);
    R_laser_(1, 1) = pow(std_laspy_, 2);

    // measurement noise covariance for radar
    R_radar_ = MatrixXd::Zero(n_radar, n_radar);
    R_radar_(0, 0) = pow(std_radr_, 2);
    R_radar_(1, 1) = pow(std_radphi_, 2);
    R_radar_(2, 2) = pow(std_radr_, 2);

    // Identity matrix for laser update
    I = MatrixXd::Identity(n_x_, n_x_);

    // NIS file names
    radar_file_ = "Radar_NIS_Dataset1.csv";
    // radar_file_ = "Radar_NIS_Dataset1_LASUKF.csv";
    laser_ukf_file_ = "Lidar_UKF_NIS_Dataset1.csv";
    laser_kf_file_ = "Lidar_KF__NIS_Dataset1.csv";

    // radar_file_ = "Radar_NIS_Dataset2.csv";
    // radar_file_ = "Radar_NIS_Dataset2_LASUKF.csv";
    // laser_ukf_file_ = "Lidar_UKF_NIS_Dataset2.csv";
    // laser_kf_file_ = "Lidar_KF__NIS_Dataset2.csv";
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // cout << "input: " << meas_package.raw_measurements_ << endl;
    if(!is_initialized_){
        // this is the first measurement
        // cout << "this is the first measurement" << endl;
        if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
            // proceed with Laser data
            // set values in x_ state vector
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);

        }
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
            // proceed with Radar data (rho, phi, rho_dot)
            // made constants to prevent modification
            const double rho = meas_package.raw_measurements_(0);
            const double phi = meas_package.raw_measurements_(1);
            const double dot = meas_package.raw_measurements_(2);

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
        // cout << "this is the 2nd or subsequent measurement" << endl;

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
            // cout << "updating for Laser" << endl;
            // use linear Kalman filter update like in EKF
            // should yield the same results with less computation.
            if(unscented_) UpdateLidar(meas_package);
            else UpdateLidarLKF(meas_package);
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
            // updating for Radar
            // cout << "updating for Radar" << endl;


            UpdateRadar(meas_package);
        }
    }

    // setting the previous timestamp
    time_us_ = meas_package.timestamp_;

    // write NIS to file
    WriteToFile(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // cout << "UKF::Prediction" << endl;
    // cout << "x_: \n" << x_ << "\n" << endl;
    // cout << "P_: \n" << P_ << "\n" << endl;
    // generating augmented sigma points
    // matrix for generated sigma points
    MatrixXd Xsig_gen_ = MatrixXd::Zero(n_aug_, n_sig);
    // generating augmented sigma points
    // returns a 7 dimentional vector
    GenerateAugmentedSigma(Xsig_gen_);

    // predict sigma points by applying process model
    PredictSigmaPoints(Xsig_gen_, delta_t);

    // predict the mean and covariance
    CalculateMeanCovariance();

    // cout << "x_: \n" << x_ << "\n" << endl;
    // cout << "P_: \n" << P_ << "\n" << endl;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // cout << "UKF::UpdateLidar" << endl;
    // measurement generated sigma points

    // Below is the code for running the unscented transformation for the Lidar measurements
    // this is rele not necessary and inefficient
    MatrixXd Zsig = MatrixXd::Zero(n_laser, n_sig);
    // predicted measurement mean
    VectorXd zpred = VectorXd::Zero(n_laser);
    // measurement covariance matrix
    MatrixXd S = MatrixXd::Zero(n_laser, n_laser);
    // cross correllation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_laser);
    // predict laser measurements
    PredictLaserMeasurements(Zsig, zpred, S, Tc);
    //update state
    UpdateStateLaser(meas_package, Zsig, zpred, S, Tc);
    cout << "NIS_laser_: " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * Linear Kalman Filter
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidarLKF(MeasurementPackage meas_package) {
    // cout << "UKF::UpdateLidarLKF" << endl;
    // measurement generated sigma points

    // Below is the code for running the Linear Kalman Filter for the Lidar measurements
    // predicted measurement mean
    VectorXd zpred = VectorXd::Zero(n_laser);
    // predict laser measurements
    PredictLaserMeasurements(meas_package, zpred);
    //update state
    UpdateStateLaser(zpred);
    cout << "NIS_laser_: " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // cout << "UKF::UpdateRadar" << endl;
    // measurement generated sigma points
    MatrixXd Zsig = MatrixXd::Zero(n_radar, n_sig);

    // predicted measurement mean
    VectorXd zpred = VectorXd::Zero(n_radar);

    // measurement covariance matrix
    MatrixXd S = MatrixXd::Zero(n_radar, n_radar);

    // cross correllation matrix
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_radar);

    // predict radar measurements
    PredictRadarMeasurements(Zsig, zpred, S, Tc);

    //update state
    UpdateStateRadar(meas_package, Zsig, zpred, S, Tc);

    cout << "NIS_radar_: " << NIS_radar_ << endl;
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
    // cout << "UKF::GenerateAugmentedSigma" << endl;
    //generate augmented state vector
    VectorXd x_aug = VectorXd::Zero(n_aug_);

    // cout << "x_aug: \n" << x_aug << "\n" << endl;
    // cout << "x_aug: \n" << x_aug.size() << "\n" << endl;

    x_aug.head(n_x_) = x_;

    // cout << "x_aug: \n" << x_aug << "\n" << endl;

    // generate augmented state covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = pow(std_a_, 2);
    P_aug(6, 6) = pow(std_yawdd_, 2);
    // cout << "P_aug: \n" << P_aug << "\n" << endl;

    // retrieves A matrix of A x A^T
    MatrixXd A = P_aug.llt().matrixL();

    /*
    cout << "A: \n" << A << "\n" << endl;
    cout << "COLS: \n" << A.cols() << "\n" << endl;
    cout << "ROWS: \n" << A.rows() << "\n" << endl;
    */

    // set the columns of sigma points matrix
    Xsig_gen_.col(0) = x_aug;
    // get the remaining sigma points
    // moved outside loop to reduce the number of redundant computations
    double shared = sqrt(lambda_ + n_aug_);
    for(int a = 0; a < n_aug_; a++){
        // common value
        VectorXd var = shared * A.col(a);
        // cout << "var: \n" << var << "\n" << endl;

        Xsig_gen_.col(a + 1) = x_aug + var;
        // cout << "Xsig_gen_.col(" << (a + 1) << "): \n" << Xsig_gen_.col(a + 1) << "\t" << endl;

        Xsig_gen_.col(a + 1 + n_aug_) = x_aug - var;
        // cout << "Xsig_gen_.col(" << (a + 1 + n_aug_) << "): \n" << Xsig_gen_.col(a + 1 + n_x_) << "\t" << endl;
    }

    // cout << "Xsig_gen_: \n" << Xsig_gen_ << "\n" << endl;
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
    // cout << "UKF::PredictSigmaPoints" << endl;
    // could re/initialize or empty Xsig_pred_ here
    // but it is better to just overwrite the values

    //predict sigma values
    for(int a = 0; a < n_sig; a++){
        // extract all 7 features from generated sigma points
        // defining as constants to prevent accidental changes in the current scope
        const double p_x = Xsig_gen_(0, a);
        const double p_y = Xsig_gen_(1, a);
        const double v = Xsig_gen_(2, a);
        const double yaw = Xsig_gen_(3, a);
        const double yaw_d = Xsig_gen_(4, a);
        // process noise values
        const double nu_a = Xsig_gen_(5, a);
        const double nu_yaw_dd = Xsig_gen_(6, a);

        // predicted state(x and y), speed, yaw angle, yaw rate

        // adding noise
        double shared = 0.5 * pow(delta_t, 2);

        double px_p = p_x + nu_a * shared * cos(yaw);
        double py_p = p_y + nu_a * shared * sin(yaw);
        double v_p = v + nu_a * delta_t;
        double yaw_p = yaw + yaw_d * delta_t + shared * nu_yaw_dd;
        double yaw_d_p = yaw_d + nu_yaw_dd * delta_t;

        //avoid division by zero
        if(fabs(yaw_d) > 0.001){ // TODO consider making this a global for more control
            // turning
            //shared value
            double inner = yaw + yaw_d * delta_t;

            px_p += v / yaw_d * (sin(inner) - sin(yaw));
            py_p += v / yaw_d * (cos(yaw) - cos(inner));
        }
        else{
            // going straight
            px_p += v * delta_t * cos(yaw);
            py_p += v * delta_t * sin(yaw);
        }

        //save predicted sigma points
        Xsig_pred_(0, a) = px_p;
        Xsig_pred_(1, a) = py_p;
        Xsig_pred_(2, a) = v_p;
        Xsig_pred_(3, a) = yaw_p;
        Xsig_pred_(4, a) = yaw_d_p;
    }

    // cout << "Xsig_pred_: \n" << Xsig_pred_ << "\n" << endl;
}

/**
 * Method for calculating the predicted state mean and covariance
 */
void UKF::CalculateMeanCovariance(){
    // cout << "UKF::CalculateMeanCovariance" << endl;
    // resetting the values in the state vector and the covariance vector because we no longer need those values
    x_.fill(0);
    P_.fill(0);

    // cout << "x_:\n" << x_ << endl;

    // predicted state mean

    /*
    for(int a = 0; a < n_sig; a++){
        x_ = x_ + weights_(a) * Xsig_pred_.col(a);
    }
     */

    // simpler but more resource intensive it seems
    // (5X1) + (5X15) * (15X1) = (5X1)
    x_ = x_ + Xsig_pred_ * weights_;

    // cout << "x_:\n" << x_ << endl;
    // cout << "Xsig_pred_:\n" << Xsig_pred_ << endl;

    // predicted state covariance
    for(int a = 0; a < n_sig; a++){
        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;

        // angle normalization to keep the filter stable and accurate
        // while(x_diff(3) > M_PI) x_diff(3) -= 2.* M_PI;
        // while(x_diff(3) < -M_PI) x_diff(3) += 2.* M_PI;
        NormalizeAngle(x_diff(3));

        // cout << "normalized x_diff:\n\t" << x_diff << endl;

        P_ = P_ + weights_(a) * x_diff * x_diff.transpose();
    }

    // cout << "P_:\n" << P_ << endl;
}

// linear kalman filter

/**
 * Method for transforming predicted state into measurement state for Laser
 * Linear Kalman Filter
 *
 * @param meas_package
 * @param zpred
 */
void UKF::PredictLaserMeasurements(MeasurementPackage meas_package, VectorXd &zpred){
    // cout << "UKF::PredictLaserMeasurements: " << endl;
    // (2 x 1) = (2 x 1) - (2 x 5) * (5 x 1)
    zpred = meas_package.raw_measurements_ - H_ * x_;

    // cout << "zpred: \n" << zpred << endl;
}

/**
 * Calculates Kalman Gain and NIS for Laser
 * Linear Kalman Filter
 *
 * @param zpred
 */
void UKF::UpdateStateLaser(VectorXd &zpred){
    // cout << "UKF::UpdateStateLaser: " << endl;

    // Applying the Linear Kalman Filter here
    // H transpose
    MatrixXd HT = H_.transpose();


    // cout << "R_laser_:\n" << R_laser_ << endl;
    // cout << "H_:\n" << H_ << endl;
    // cout << "P_:\n" << P_ << endl;
    // cout << "HT:\n" << HT << endl;


    // project system uncertainty into measurement space
    // measurement noise covariance + measurement projection function
    MatrixXd S = R_laser_ + (H_ * P_ * HT);

    // cout << "S:\n" << S << endl;

    // S Inverse
    MatrixXd SI = S.inverse();

    // cout << "SI:\n" << SI << endl;

    // Kalman Filter gain
    MatrixXd K = P_ * HT * SI;

    // cout << "K:\n" << K << endl;

    // x_ estimate
    x_ = x_ + K * zpred;

    // cout << "x_:\n" << x_ << endl;

    // cout << "I:\n" << I << endl;

    // uncertainty matrix
    P_ = (I - K * H_) * P_;

    // cout << "P_:\n" << P_ << endl;

    // calculate NIS
    // returned so that the corresponding tool can record its NIS score
    NIS_laser_ =  zpred.transpose() * S.inverse() * zpred;
}

/**
 * Method for transforming predicted state into measurement state for Laser and cross correlation matrix
 *
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 */
void UKF::PredictLaserMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // cout << "UKF::PredictLaserMeasurements" << endl;
    // transform sigma points into measurement space
    for(int a = 0; a < n_sig; a++){
        // extract values
        Zsig(0, a) = Xsig_pred_(0, a);
        Zsig(1, a) = Xsig_pred_(1, a);
    }

    // mean predicted measurement
    /*
    for(int a = 0; a < n_sig; a++){
        zpred += weights_(a) * Zsig.col(a);
    }
     */
    //(2X1) + (2X15) * (15X1) = (2X1)
    zpred = zpred + Zsig * weights_;

    // measurement covariance matrix
    for(int a = 0; a < n_sig; a++){
        // residual
        VectorXd z_diff = Zsig.col(a) - zpred;

        S = S + weights_(a) * z_diff * z_diff.transpose();

        // calculate cross corellation

        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;
        //angle normalization
        // while(x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        // while(x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
        NormalizeAngle(x_diff(3));

        Tc = Tc + weights_(a) * x_diff * z_diff.transpose();
    }

    // adding noise
    S = S + R_laser_;

    // cout << "S.cols(): " << S.cols() << endl;
    // cout << "S.rows(): " << S.rows() << endl;
}

/**
 * Calculates Kalman Gain and NIS for Laser
 *
 * @param meas_package
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 * @return NIS value as double
 */
void UKF::UpdateStateLaser(MeasurementPackage meas_package, MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // cout << "UKF::UpdateState" << endl;

    // Kalman Gain
    MatrixXd K = Tc * S.inverse();

    VectorXd z_diff = meas_package.raw_measurements_ - zpred;

    // update state mean and covariance
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // calculate NIS
    // returned so that the corresponding tool can record its NIS score
    NIS_laser_ =  z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Method for transforming predicted state into measurement state for Radar and cross correlation matrix
 *
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 */
void UKF::PredictRadarMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // cout << "UKF::PredictRadarMeasurements" << endl;
    // transform sigma points into measurement space
    // cout << "transform sigma points into measurement space" << endl;
    for(int a = 0; a < n_sig; a++){
        // extract values
        double p_x = Xsig_pred_(0, a);
        double p_y = Xsig_pred_(1, a);
        double v = Xsig_pred_(2, a);
        double yaw = Xsig_pred_(3, a);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        //check for zeros
        /*
        if (fabs(p_x) < 0.001) {
            p_x = 0.001;
        }
        if (fabs(p_y) < 0.001) {
            p_y = 0.001;
        }
        */

        double shared = sqrt(pow(p_x, 2) + pow(p_y, 2));

        // measurement model
        Zsig(0, a) = shared; // rho
        Zsig(1, a) = atan2(p_y, p_x); // phi
        Zsig(2, a) = (p_x * v1 + p_y * v2) / shared; // rho dot

        // cout << "\tZsig.col(" << a << "): " << Zsig.col(a) << "\n" << endl;
    }

    // cout << "Zsig: " << Zsig << endl;

    // mean predicted measurement
    // cout << "mean predicted measurement" << endl;
    // cout << "zpred: " << zpred << endl;

    /*
    for(int a = 0; a < n_sig; a++){
        zpred = zpred + weights_(a) * Zsig.col(a);
    }
     */
    //(3X15) * (15X1) = (3X1)
    zpred = zpred + Zsig * weights_;

    // cout << "zpred: " << zpred << endl;

    // measurement covariance matrix
    // cout << "measurement covariance matrix" << endl;
    for(int a = 0; a < n_sig; a++){
        // residual
        VectorXd z_diff = Zsig.col(a) - zpred;

        // cout << "z_diff: " << z_diff << endl;

        // angle normalization
        // while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        // while(z_diff(1) < - M_PI) z_diff(1) += 2. * M_PI;
        NormalizeAngle(z_diff(1));

        // cout << "z_diff: " << z_diff << endl;

        S = S + weights_(a) * z_diff * z_diff.transpose();

        // calculate cross corellation

        // state difference
        VectorXd x_diff = Xsig_pred_.col(a) - x_;
        //angle normalization
        // while(x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        // while(x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
        NormalizeAngle(x_diff(3));

        Tc = Tc + weights_(a) * x_diff * z_diff.transpose();
    }

    // add measurement noise

    S = S + R_radar_;
}

/**
 * Calculates Kalman Gain and NIS for Radar
 *
 * @param meas_package
 * @param Zsig
 * @param zpred
 * @param S
 * @param Tc
 */
void UKF::UpdateStateRadar(MeasurementPackage meas_package, MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc){
    // cout << "UKF::UpdateState" << endl;

    // Kalman Gain
    MatrixXd K = Tc * S.inverse();

    VectorXd z_diff = meas_package.raw_measurements_ - zpred;

    //angle normalization
    // while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    // while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    NormalizeAngle(z_diff(1));

    // update state mean and covariance
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // calculate NIS
    // returned so that the corresponding tool can record its NIS score
    NIS_radar_ =  z_diff.transpose() * S.inverse() * z_diff;
}

/**
     * Method for normalizing angles
     *
     * @param phi
     */
void UKF::NormalizeAngle(double& phi){
    phi = atan2(sin(phi), cos(phi));
}

/**
 * Method for writing NIS values to file
 *
 * @param meas_package
 */
void UKF::WriteToFile(MeasurementPackage meas_package){
 if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
     // write to radar file
     // open file
     file.open(radar_file_);

     // write to file
     file << NIS_radar_ << " ";
 }
 else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
     if(unscented_){
         // write to unscented laser file
         file.open(laser_ukf_file_);
     }
     else{
         // write to linear kalman filter file
         file.open(laser_kf_file_);
     }

     // write to file
     file << NIS_laser_ << " ";
 }

// close file
file.close();
}