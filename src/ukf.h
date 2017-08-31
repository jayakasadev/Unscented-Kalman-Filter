#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;

    ///* state covariance matrix
    MatrixXd P_;

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    ///* time when the state is true, in us
    long long time_us_;

    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;

    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;

    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;

    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;

    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_ ;

    ///* Weights of sigma points
    VectorXd weights_;

    ///* State dimension
    int n_x_;

    ///* Augmented state dimension
    int n_aug_;

    ///* Sigma point spreading parameter
    double lambda_;

    // number of sigma points
    short n_sig;

    // number of radar input values
    short n_radar;

    //number of laser input values
    short n_laser;

    // NIS for radar
    double NIS_radar_;

    // NIS for laser
    double NIS_laser_;

    // turn unscented transformation on or off and use linear kalman filter for laser
    bool unscented_;

    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);

    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);

private:
    /**
     * Method for generating the Augmented Sigma Points
     *
     * @param Xsig_gen_
     */
    void GenerateAugmentedSigma(MatrixXd &Xsig_gen_);

    /**
     * Method for predicting sigma points by applying process model to generated sigma points
     * We are predicting where the sigma points will be at time k+1 given the generated sigma
     * points at time k
     *
     * Takes a 7x15 matrix in and outputs a 5x15 matrix in Xsig_pred_
     *
     * @param Xsig_gen_
     * @param d
     */
    void PredictSigmaPoints(MatrixXd &Xsig_gen_, const double &delta_t);

    /**
     * Method for calculating the predicted state mean and covariance
     */
    void CalculateMeanCovariance();

    // linear kalman filter for laser

    /**
     * Method for transforming predicted state into measurement state for Laser
     */
    void PredictLaserMeasurements();

    /**
     * Calculates Kalman Gain and NIS for Laser
     *
     * @param meas_package
     */
    void UpdateStateLaser(MeasurementPackage meas_package);

    // unscented transformation for lidar

    /**
     * Method for transforming predicted state into measurement state for Laser and cross correlation matrix
     *
     * @param Zsig
     * @param zpred
     * @param S
     * @param Tc
     */
    void PredictLaserMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc);

    /**
     * Calculates Kalman Gain and NIS for Laser
     *
     * @param meas_package
     * @param Zsig
     * @param zpred
     * @param S
     * @param Tc
     */
    void UpdateStateLaser(MeasurementPackage meas_package, MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc);

    // unscented transformation for radar

    /**
     * Method for transforming predicted state into measurement state for Radar and cross correlation matrix
     *
     * @param Zsig
     * @param zpred
     * @param S
     * @param Tc
     */
    void PredictRadarMeasurements(MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc);

    /**
     * Calculates Kalman Gain and NIS for Radar
     *
     * @param meas_package
     * @param Zsig
     * @param zpred
     * @param S
     * @param Tc
     */
    void UpdateStateRadar(MeasurementPackage meas_package, MatrixXd &Zsig, VectorXd &zpred, MatrixXd &S, MatrixXd &Tc);

    /**
     * Method for normalizing angles
     *
     * @param phi
     */
    void NormalizeAngle(double& phi);

};

#endif /* UKF_H */
