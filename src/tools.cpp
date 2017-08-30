#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

/**
 * Method for calculating the RMSE for given inputs
 *
 * @param estimations
 * @param ground_truth
 * @return
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
    // output vector
    VectorXd rmse(4);
    rmse.fill(0);

    // if there are no estimations or the number of estimations do not match the number of
    // ground truths, return empty rmse
    if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
        cout << "\n\tError: Invalid estimation or ground_truth data" << endl;
        return rmse;
    }


    // accumulate squared residuals
    for(unsigned int  i = 0; i < estimations.size(); i++){
        // getting the difference between the calculated and expected
        VectorXd residual = estimations[i] - ground_truth[i];

        residual = residual.array().square();
        rmse += residual;
    }

    // averaging over the number of estimates and returning the square root
    rmse = (rmse / estimations.size()).array().sqrt();

    return rmse;
}