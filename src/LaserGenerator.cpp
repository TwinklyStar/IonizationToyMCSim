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
//    yaw=0; pitch=0; roll=0;
//    yaw_355=0; pitch_355=0; roll_355=0;
//    peak_time=0; obe_ptr=nullptr;
//
//    rot_mat_122 = Eigen::Matrix3d::Identity();
//    rot_mat_355 = Eigen::Matrix3d::Identity();
//    rot_mat_rev_122 = Eigen::Matrix3d::Identity();

}

TVector3 LaserGenerator::GetFieldE(TVector3 r, Double_t t) {
    double Ex=0, Ey=0, Ez=0;
    for (auto lsr : vec_laser122){
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();
        Double_t z = laser_r.Z();

        const Double_t eta = 376.7303134;

        Double_t prefactor = sqrt(2.0 * sqrt(2.0) * eta * lsr.energy / (pow(TMath::Pi(), 3.0 / 2.0) * lsr.sigma_x * lsr.sigma_y  * lsr.tau *
                                                                    pow(10, -9)));      // convert ns to s

        // Compute the spatial Gaussian components
        double exp_space = exp(-(x * x) / (lsr.sigma_x * lsr.sigma_x) - (y * y) / (lsr.sigma_y * lsr.sigma_y));

        // Compute the temporal Gaussian component
        double exp_time = exp(-((t-lsr.peak_time) * (t-lsr.peak_time)) / (4.0 * lsr.tau * lsr.tau));

        // Combine everything to compute the electric field
        double E = prefactor * exp_space * exp_time;

        Eigen::MatrixXd before_rot(3,1);
        before_rot << E,
                      0,
                      0;

        Eigen::MatrixXd after_rot = lsr.rot_mat_rev*before_rot;

        Ex += after_rot(2, 0);
        Ey += after_rot(0, 0);
        Ez += after_rot(1, 0);

    }
//    double E = prefactor * exp_space;

    return {Ex, Ey, Ez};   // x-polarized in muonium coordinate
}

Double_t LaserGenerator::GetPeakIntensity(TVector3 r) {
    Double_t sum_I=0;
    for (auto lsr : vec_laser122) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();
        Double_t z = laser_r.Z();

        const Double_t eta = 376.7303134;

        // Define the prefactor
        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * pow(10, -9));

        // Compute the spatial Gaussian components
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));

        // Combine everything to compute the intensity
        double I = prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2

        sum_I += I;
    }
    return sum_I; // Return the intensity
}

Double_t LaserGenerator::GetPeakIntensity355(TVector3 r) {
    Double_t sum_I=0;
    for (auto lsr : vec_laser355) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();
        Double_t z = laser_r.Z();

        // Define the prefactor
        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * pow(10, -9));

        // Compute the spatial Gaussian components
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));

        // Combine everything to compute the intensity
        double I = prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2

        sum_I += I;
    }
    return sum_I; // Return the intensity
}

Double_t LaserGenerator::GetIntensity(TVector3 r, Double_t t) {
    Double_t sum_I=0;
    for (auto lsr : vec_laser122) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();
        Double_t z = laser_r.Z();

        // Define the prefactor
        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * pow(10, -9));

        // Compute the spatial Gaussian components
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));

        // Compute the temporal Gaussian component
        double exp_time = exp(-((t - lsr.peak_time) * (t - lsr.peak_time)) / (2.0 * lsr.tau * lsr.tau));

        // Combine everything to compute the intensity
        double I = prefactor * exp_space * exp_time * 100;

        sum_I += I;
    }
    return sum_I; // Return the intensity
}

Double_t LaserGenerator::GetIntensity355(TVector3 r, Double_t t) {
    Double_t sum_I=0;
    for (auto lsr : vec_laser355) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();
        Double_t z = laser_r.Z();

        // Define the prefactor
        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * pow(10, -9));

        // Compute the spatial Gaussian components
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));

        // Compute the temporal Gaussian component
        double exp_time = exp(-((t - lsr.peak_time) * (t - lsr.peak_time)) / (2.0 * lsr.tau * lsr.tau));

        // Combine everything to compute the intensity
        double I = prefactor * exp_space * exp_time * 100;

        sum_I += I;
    }

    return sum_I; // Return the intensity
}

TVector3 LaserGenerator::BeamToLaserCoord(TVector3 r, Laser lsr) {
    Double_t x = r.Y() - lsr.laser_offset.Y();
    Double_t y = r.Z() - lsr.laser_offset.Z();
    Double_t z = r.X() - lsr.laser_offset.X();

    Eigen::MatrixXd before_rot(3,1);
    before_rot << x,
                  y,
                  z;

    Eigen::MatrixXd after_rot = lsr.rot_mat*before_rot;

    return {after_rot(0,0), after_rot(1,0), after_rot(2,0)};
}

//TVector3 LaserGenerator::BeamToLaserCoord355(TVector3 r) {
//    Double_t x = r.Y() - laser_offset_355.Y();
//    Double_t y = r.Z() - laser_offset_355.Z();
//    Double_t z = r.X() - laser_offset_355.X();
//
//    Eigen::MatrixXd before_rot(3,1);
//    before_rot << x,
//            y,
//            z;
//
//    Eigen::MatrixXd after_rot = rot_mat_355*before_rot;
//
//    return {after_rot(0,0), after_rot(1,0), after_rot(2,0)};
//}
//
//TVector3 LaserGenerator::GetWaveVector() {
//    return laser_k*laser_dirc;
//    // z'->x. x'->y, y'->z
//
//}

void LaserGenerator::UpdateRotMat(Laser &lsr) {
    Eigen::Matrix3d Y, P, R;
    Y << 1, 0, 0,
            0, TMath::Cos(lsr.yaw), TMath::Sin(lsr.yaw),
            0, -TMath::Sin(lsr.yaw), TMath::Cos(lsr.yaw);
    P << TMath::Cos(lsr.pitch) ,0, -TMath::Sin(lsr.pitch),
            0, 1, 0,
            TMath::Sin(lsr.pitch), 0, TMath::Cos(lsr.pitch);
    R << TMath::Cos(lsr.roll), TMath::Sin(lsr.roll), 0,
            -TMath::Sin(lsr.roll), TMath::Cos(lsr.roll), 0,
            0, 0, 1;
    lsr.rot_mat = R*P*Y;

    Eigen::Matrix3d Y_rev, P_rev, R_rev;
    Y_rev << 1, 0, 0,
             0, TMath::Cos(lsr.yaw), -TMath::Sin(lsr.yaw),
             0, TMath::Sin(lsr.yaw), TMath::Cos(lsr.yaw);
    P_rev << TMath::Cos(lsr.pitch) ,0, TMath::Sin(lsr.pitch),
             0, 1, 0,
             -TMath::Sin(lsr.pitch), 0, TMath::Cos(lsr.pitch);
    R_rev << TMath::Cos(lsr.roll), -TMath::Sin(lsr.roll), 0,
             TMath::Sin(lsr.roll), TMath::Cos(lsr.roll), 0,
             0, 0, 1;
    lsr.rot_mat_rev = Y_rev*P_rev*R_rev;

    Eigen::MatrixXd laser_dirc_beforerot(3,1);
    laser_dirc_beforerot << 0,
            0,
            1;
    Eigen::MatrixXd after_rot = lsr.rot_mat_rev*laser_dirc_beforerot;

    lsr.laser_dirc = {after_rot(2,0), after_rot(0,0), after_rot(1,0)};

//    std::cout << "--- Rotation matrix updated.\n122nm:\n" << rot_mat_122 << "\n355nm:\n" << rot_mat_355 << std::endl;
//    std::cout << "The 122nm wave vector direction:\n" << laser_dirc.X() << ", " << laser_dirc.Y() << ", " << laser_dirc.Z() << std::endl;
}

void LaserGenerator::AddLaser122(Double_t energy, Double_t pulse_FWHM, Double_t peak_time, Double_t linewidth,
                                 Double_t sigma_x, Double_t sigma_y, Double_t offset_x, Double_t offset_y,
                                 Double_t offset_z, Double_t yaw, Double_t pitch, Double_t roll) {
    Laser lsr_tmp;
    lsr_tmp.energy = energy;
    lsr_tmp.linewidth = linewidth;
    lsr_tmp.peak_time = peak_time;
    lsr_tmp.sigma_x = sigma_x;
    lsr_tmp.sigma_y = sigma_y;
    lsr_tmp.tau = 0.4247 * pulse_FWHM;
    lsr_tmp.laser_offset = {offset_x, offset_y, offset_z};
    lsr_tmp.yaw = yaw * TMath::Pi() / 180;
    lsr_tmp.pitch = pitch * TMath::Pi() / 180;
    lsr_tmp.roll = roll * TMath::Pi() / 180;
    lsr_tmp.wavelength = 122;
    lsr_tmp.laser_k = 2*TMath::Pi()/lsr_tmp.wavelength*1e9;
    lsr_tmp.cen_freq = 299792458/lsr_tmp.wavelength;

    UpdateRotMat(lsr_tmp);

    // Output the details of the laser parameters
    std::ostringstream oss;
    oss << "\n-- Add 122nm laser:\n"
        << "  Energy: " << lsr_tmp.energy << " J\n"
        << "  Pulse FWHM: " << pulse_FWHM << " ns\n"
        << "  Peak Time: " << lsr_tmp.peak_time << " ns\n"
        << "  Linewidth: " << lsr_tmp.linewidth << " GHz\n"
        << "  Sigma X: " << lsr_tmp.sigma_x << " mm\n"
        << "  Sigma Y: " << lsr_tmp.sigma_y << " mm\n"
        << "  Tau: " << lsr_tmp.tau << " ns\n"
        << "  Offset (X, Y, Z): (" << lsr_tmp.laser_offset.X() << ", "
        << lsr_tmp.laser_offset.Y() << ", "
        << lsr_tmp.laser_offset.Z() << ") mm\n"
        << "  Yaw: " << yaw << " deg\n"
        << "  Pitch: " << pitch << " deg\n"
        << "  Roll: " << roll << " deg\n"
//        << "  Wavevector (k): " << lsr_tmp.laser_k << " m^-1\n"
        << "  Wavevector direction (unit vector): (" << lsr_tmp.laser_dirc.X() << ", "
        << lsr_tmp.laser_dirc.Y() << ", "
        << lsr_tmp.laser_dirc.Z() << ")\n"
        << "  Rotation matrix: " << lsr_tmp.rot_mat << std::endl;

    cout << oss.str();

    vec_laser122.push_back(lsr_tmp);
}

void LaserGenerator::AddLaser355(Double_t energy, Double_t pulse_FWHM, Double_t peak_time, Double_t linewidth,
                                 Double_t sigma_x, Double_t sigma_y, Double_t offset_x, Double_t offset_y,
                                 Double_t offset_z, Double_t yaw, Double_t pitch, Double_t roll) {
    Laser lsr_tmp;
    lsr_tmp.energy = energy;
    lsr_tmp.linewidth = linewidth;
    lsr_tmp.peak_time = peak_time;
    lsr_tmp.sigma_x = sigma_x;
    lsr_tmp.sigma_y = sigma_y;
    lsr_tmp.tau = 0.4247 * pulse_FWHM;
    lsr_tmp.laser_offset = {offset_x, offset_y, offset_z};
    lsr_tmp.yaw = yaw * TMath::Pi() / 180;
    lsr_tmp.pitch = pitch * TMath::Pi() / 180;
    lsr_tmp.roll = roll * TMath::Pi() / 180;
    lsr_tmp.wavelength = 122;
    lsr_tmp.laser_k = 2*TMath::Pi()/lsr_tmp.wavelength*1e9;
    lsr_tmp.cen_freq = 299792458/lsr_tmp.wavelength;

    UpdateRotMat(lsr_tmp);

    // Output the details of the laser parameters
    std::ostringstream oss;
    oss << "\n-- Add 355nm laser:\n"
        << "  Energy: " << lsr_tmp.energy << " J\n"
        << "  Pulse FWHM: " << pulse_FWHM << " ns\n"
        << "  Peak Time: " << lsr_tmp.peak_time << " ns\n"
        << "  Linewidth: " << lsr_tmp.linewidth << " GHz\n"
        << "  Sigma X: " << lsr_tmp.sigma_x << " mm\n"
        << "  Sigma Y: " << lsr_tmp.sigma_y << " mm\n"
        << "  Tau: " << lsr_tmp.tau << " ns\n"
        << "  Offset (X, Y, Z): (" << lsr_tmp.laser_offset.X() << ", "
        << lsr_tmp.laser_offset.Y() << ", "
        << lsr_tmp.laser_offset.Z() << ") mm\n"
        << "  Yaw: " << yaw << " deg\n"
        << "  Pitch: " << pitch << " deg\n"
        << "  Roll: " << roll << " deg\n"
//        << "  Wavevector (k): " << lsr_tmp.laser_k << " m^-1\n"
        << "  Wavevector direction (unit vector): (" << lsr_tmp.laser_dirc.X() << ", "
        << lsr_tmp.laser_dirc.Y() << ", "
        << lsr_tmp.laser_dirc.Z() << ")\n"
        << "  Rotation matrix: " << lsr_tmp.rot_mat << std::endl;

    cout << oss.str();

    vec_laser355.push_back(lsr_tmp);
}
