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

LaserGenerator::LaserGenerator() {
    yaw=0; pitch=0; roll=0;
    yaw_355=0; pitch_355=0; roll_355=0;
    peak_time=0; obe_ptr=nullptr;

    rot_mat_122 = Eigen::Matrix3d::Identity();
    rot_mat_355 = Eigen::Matrix3d::Identity();
    rot_mat_rev_122 = Eigen::Matrix3d::Identity();

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

    Eigen::MatrixXd after_rot = rot_mat_122*before_rot;

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

    Eigen::MatrixXd after_rot = rot_mat_355*before_rot;

    return {after_rot(0,0), after_rot(1,0), after_rot(2,0)};
}

TVector3 LaserGenerator::GetWaveVector() {
    return laser_k*laser_dirc;
    // z'->x. x'->y, y'->z

}

void LaserGenerator::UpdateRotMat() {
    Eigen::Matrix3d Y, P, R;
    Y << 1, 0, 0,
            0, TMath::Cos(yaw), TMath::Sin(yaw),
            0, -TMath::Sin(yaw), TMath::Cos(yaw);
    P << TMath::Cos(pitch) ,0, -TMath::Sin(pitch),
            0, 1, 0,
            TMath::Sin(pitch), 0, TMath::Cos(pitch);
    R << TMath::Cos(roll), TMath::Sin(roll), 0,
            -TMath::Sin(roll), TMath::Cos(roll), 0,
            0, 0, 1;
    rot_mat_122 = R*P*Y;

    Eigen::Matrix3d Y_355, P_355, R_355;
    Y_355 << 1, 0, 0,
            0, TMath::Cos(yaw_355), TMath::Sin(yaw_355),
            0, -TMath::Sin(yaw_355), TMath::Cos(yaw_355);
    P_355 << TMath::Cos(pitch_355) ,0, -TMath::Sin(pitch_355),
            0, 1, 0,
            TMath::Sin(pitch_355), 0, TMath::Cos(pitch_355);
    R_355 << TMath::Cos(roll_355), TMath::Sin(roll_355), 0,
            -TMath::Sin(roll_355), TMath::Cos(roll_355), 0,
            0, 0, 1;
    rot_mat_355 = R_355*P_355*Y_355;

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
    rot_mat_rev_122 = Y_rev*P_rev*R_rev;

    Eigen::MatrixXd laser_dirc_beforerot(3,1);
    laser_dirc_beforerot << 0,
            0,
            1;
    Eigen::MatrixXd after_rot = rot_mat_rev_122*laser_dirc_beforerot;

    laser_dirc = {after_rot(2,0), after_rot(0,0), after_rot(1,0)};

    std::cout << "--- Rotation matrix updated.\n122nm:\n" << rot_mat_122 << "\n355nm:\n" << rot_mat_355 << std::endl;
    std::cout << "The 122nm wave vector direction:\n" << laser_dirc.X() << ", " << laser_dirc.Y() << ", " << laser_dirc.Z() << std::endl;
}
