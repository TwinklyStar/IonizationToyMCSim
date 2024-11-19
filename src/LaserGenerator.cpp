//
// Created by Meng Lv on 2024/10/9.
//

#include "LaserGenerator.h"
#include "OBEsolver.h"

// Meyers' Singleton implementation
LaserGenerator& LaserGenerator::GetInstance() {
    static LaserGenerator instance;
    return instance;
}

TVector3 LaserGenerator::GetFieldE(TVector3 r, Double_t t) {
    TVector3 laser_r = BeamToLaserCoord(r);
    Double_t x = laser_r.X();
    Double_t y = laser_r.Y();
    Double_t z = laser_r.Z();
    const Double_t eta = 376.7303134;

    Double_t prefactor = sqrt(2.0 * sqrt(2.0) * eta * energy / (pow(TMath::Pi(), 3.0 / 2.0) * sigma_x * sigma_y  * tau *
            pow(10, -9)));      // convert ns to s

    // Compute the spatial Gaussian components
    double exp_space = exp(-(x * x) / (sigma_x * sigma_x) - (y * y) / (sigma_y * sigma_y));

    // Compute the temporal Gaussian component
    double exp_time = exp(-((t-peak_time) * (t-peak_time)) / (4.0 * tau * tau));

    // Combine everything to compute the electric field
    double E = prefactor * exp_space * exp_time;
//    double E = prefactor * exp_space;

    return {0, E, 0};   // x-polarized in muonium coordinate
}

Double_t LaserGenerator::GetPeakIntensity(TVector3 r) {
    TVector3 laser_r = BeamToLaserCoord(r);
    Double_t x = laser_r.X();
    Double_t y = laser_r.Y();
    Double_t z = laser_r.Z();
    const Double_t eta = 376.7303134;
//
//    Double_t prefactor = sqrt(2.0 * sqrt(2.0) * eta * energy / (pow(TMath::Pi(), 3.0 / 2.0) * sigma_x * sigma_y  * tau
//            * pow(10, -9)));      // convert ns to s
//    // Compute the spatial Gaussian components
//    double exp_space = exp(-(x * x) / (sigma_x * sigma_x) - (y * y) / (sigma_y * sigma_y));
//
//    double E = prefactor * exp_space;
//
//    return pow(E, 2)/(2*eta)*100;
    // Define the prefactor
    double prefactor = (2.0 * energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * sigma_x * sigma_y * tau * pow(10, -9));

    // Compute the spatial Gaussian components
    double exp_space = exp(-2.0 * (x * x) / (sigma_x * sigma_x) - 2.0 * (y * y) / (sigma_y * sigma_y));

    // Combine everything to compute the intensity
    double I = prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2

    return I; // Return the intensity
}

Double_t LaserGenerator::GetPeakIntensity355(TVector3 r) {
    TVector3 laser_r = BeamToLaserCoord(r);
    Double_t x = laser_r.X();
    Double_t y = laser_r.Y();
    Double_t z = laser_r.Z();
    // Define the prefactor
    double prefactor = (2.0 * energy_355) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * sigma_x * sigma_y * tau * pow(10, -9));

    // Compute the spatial Gaussian components
    double exp_space = exp(-2.0 * (x * x) / (sigma_x * sigma_x) - 2.0 * (y * y) / (sigma_y * sigma_y));

    // Combine everything to compute the intensity
    double I = prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2

    return I; // Return the intensity
}

Double_t LaserGenerator::GetIntensity(TVector3 r, Double_t t) {
    TVector3 laser_r = BeamToLaserCoord(r);
    Double_t x = laser_r.X();
    Double_t y = laser_r.Y();
    Double_t z = laser_r.Z();
    // Define the prefactor
    double prefactor = (2.0 * energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * sigma_x * sigma_y * tau * pow(10, -9));

    // Compute the spatial Gaussian components
    double exp_space = exp(-2.0 * (x * x) / (sigma_x * sigma_x) - 2.0 * (y * y) / (sigma_y * sigma_y));

    // Compute the temporal Gaussian component
    double exp_time = exp(-((t-peak_time) * (t-peak_time)) / (2.0 * tau * tau));

    // Combine everything to compute the intensity
    double I = prefactor * exp_space * exp_time * 100;

    return I; // Return the intensity
}

Double_t LaserGenerator::GetIntensity355(TVector3 r, Double_t t) {
    TVector3 laser_r = BeamToLaserCoord355(r);
    Double_t x = laser_r.X();
    Double_t y = laser_r.Y();
    Double_t z = laser_r.Z();
    // Define the prefactor
    double prefactor = (2.0 * energy_355) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * sigma_x * sigma_y * tau * pow(10, -9));

    // Compute the spatial Gaussian components
    double exp_space = exp(-2.0 * (x * x) / (sigma_x * sigma_x) - 2.0 * (y * y) / (sigma_y * sigma_y));

    // Compute the temporal Gaussian component
    double exp_time = exp(-((t-peak_time) * (t-peak_time)) / (2.0 * tau * tau));

    // Combine everything to compute the intensity
    double I = prefactor * exp_space * exp_time * 100;

    return I; // Return the intensity
}

TVector3 LaserGenerator::BeamToLaserCoord(TVector3 r) {
    Double_t x = r.Y() - laser_offset.Y();
    Double_t y = r.Z() - laser_offset.Z();
    Double_t z = r.X() - laser_offset.X();

    Eigen::MatrixXd before_rot(3,1);
    before_rot << x,
                  y,
                  z;

    Eigen::MatrixXd Y(3, 3), P(3, 3), R(3, 3);
    Y << 1, 0, 0,
         0, TMath::Cos(yaw), TMath::Sin(yaw),
         0, -TMath::Sin(yaw), TMath::Cos(yaw);
    P << TMath::Cos(pitch) ,0, -TMath::Sin(pitch),
         0, 1, 0,
         TMath::Sin(pitch), 0, TMath::Cos(pitch);
    R << TMath::Cos(roll), TMath::Sin(roll), 0,
         -TMath::Sin(roll), TMath::Cos(roll), 0,
         0, 0, 1;

    Eigen::MatrixXd after_rot = R*P*Y*before_rot;

    return {after_rot(0,0), after_rot(1,0), after_rot(2,0)};
}

TVector3 LaserGenerator::BeamToLaserCoord355(TVector3 r) {
    Double_t x = r.Y() - laser_offset_355.Y();
    Double_t y = r.Z() - laser_offset_355.Z();
    Double_t z = r.X() - laser_offset_355.X();

    Eigen::MatrixXd before_rot(3,1);
    before_rot << x,
            y,
            z;

    Eigen::MatrixXd Y(3, 3), P(3, 3), R(3, 3);
    Y << 1, 0, 0,
            0, TMath::Cos(yaw_355), TMath::Sin(yaw_355),
            0, -TMath::Sin(yaw_355), TMath::Cos(yaw_355);
    P << TMath::Cos(pitch_355) ,0, -TMath::Sin(pitch_355),
            0, 1, 0,
            TMath::Sin(pitch_355), 0, TMath::Cos(pitch_355);
    R << TMath::Cos(roll_355), TMath::Sin(roll_355), 0,
            -TMath::Sin(roll_355), TMath::Cos(roll_355), 0,
            0, 0, 1;

    Eigen::MatrixXd after_rot = R*P*Y*before_rot;

    return {after_rot(0,0), after_rot(1,0), after_rot(2,0)};
}

TVector3 LaserGenerator::GetWaveVector() {
    Eigen::MatrixXd before_rot(3,1);
    before_rot << 0,
                  0,
                  1;

    Eigen::MatrixXd Y_rev(3, 3), P_rev(3, 3), R_rev(3, 3);
    Y_rev << 1, 0, 0,
             0, TMath::Cos(yaw), -TMath::Sin(yaw),
             0, TMath::Sin(yaw), TMath::Cos(yaw);
    P_rev << TMath::Cos(pitch) ,0, TMath::Sin(pitch),
             0, 1, 0,
             -TMath::Sin(pitch), 0, TMath::Cos(pitch);
    R_rev << TMath::Cos(roll), -TMath::Sin(roll), 0,
             TMath::Sin(roll), TMath::Cos(roll), 0,
             0, 0, 1;

    Eigen::MatrixXd after_rot = Y_rev*P_rev*R_rev*before_rot;

    return {laser_k*after_rot(2,0), laser_k*after_rot(0,0), laser_k*after_rot(1,0)};
    // z'->x. x'->y, y'->z

}
