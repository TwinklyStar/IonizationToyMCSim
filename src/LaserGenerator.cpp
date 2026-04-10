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
    for (size_t i = 0; i < vec_laser122.size(); ++i) {
        const Laser& lsr = vec_laser122[i];
        // Temporal Gaussian only — spatial part precomputed in cached_Espatial_122
        double exp_time = exp(-((t - lsr.peak_time) * (t - lsr.peak_time)) / (4.0 * lsr.tau * lsr.tau));
        double E = cached_Espatial_122[i] * exp_time;

        Eigen::Vector3d before_rot(E, 0, 0);
        Eigen::Vector3d after_rot = lsr.rot_mat_rev * before_rot;

        Ex += after_rot(2);
        Ey += after_rot(0);
        Ez += after_rot(1);
    }
    return {Ex, Ey, Ez};   // x-polarized in muonium coordinate
}

Double_t LaserGenerator::GetPeakIntensity(TVector3 r) {
    Double_t sum_I=0;
    for (const auto& lsr : vec_laser122) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();

        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * 1e-9);
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));
        sum_I += prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2
    }
    return sum_I;
}

Double_t LaserGenerator::GetPeakIntensity355(TVector3 r) {
    Double_t sum_I=0;
    for (const auto& lsr : vec_laser355) {
        TVector3 laser_r = BeamToLaserCoord(r, lsr);
        Double_t x = laser_r.X();
        Double_t y = laser_r.Y();

        double prefactor =
                (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * 1e-9);
        double exp_space = exp(-2.0 * (x * x) / (lsr.sigma_x * lsr.sigma_x) - 2.0 * (y * y) / (lsr.sigma_y * lsr.sigma_y));
        sum_I += prefactor * exp_space * 100; // convert W/mm^2 to W/cm^2
    }
    return sum_I;
}

Double_t LaserGenerator::GetIntensity(TVector3 r, Double_t t) {
    Double_t sum_I=0;
    for (size_t i = 0; i < vec_laser122.size(); ++i) {
        const Laser& lsr = vec_laser122[i];
        // Temporal Gaussian only — spatial part precomputed in cached_Ispatial_122
        double exp_time = exp(-((t - lsr.peak_time) * (t - lsr.peak_time)) / (2.0 * lsr.tau * lsr.tau));
        sum_I += cached_Ispatial_122[i] * exp_time * 100;
    }
    return sum_I;
}

Double_t LaserGenerator::GetIntensity355(TVector3 r, Double_t t) {
    Double_t sum_I=0;
    for (size_t i = 0; i < vec_laser355.size(); ++i) {
        const Laser& lsr = vec_laser355[i];
        // Temporal Gaussian only — spatial part precomputed in cached_Ispatial_355
        double exp_time = exp(-((t - lsr.peak_time) * (t - lsr.peak_time)) / (2.0 * lsr.tau * lsr.tau));
        sum_I += cached_Ispatial_355[i] * exp_time * 100;
    }
    return sum_I;
}

TVector3 LaserGenerator::BeamToLaserCoord(TVector3 r, const Laser& lsr) {
    Eigen::Vector3d v(r.Y() - lsr.laser_offset.Y(),
                      r.Z() - lsr.laser_offset.Z(),
                      r.X() - lsr.laser_offset.X());
    Eigen::Vector3d vr = lsr.rot_mat * v;
    return {vr(0), vr(1), vr(2)};
}

void LaserGenerator::PrecomputeAtPosition(TVector3 r) {
    const Double_t eta = 376.7303134;
    const Double_t pi32 = TMath::Pi() * sqrt(TMath::Pi());  // π^(3/2)

    cached_Espatial_122.resize(vec_laser122.size());
    cached_Ispatial_122.resize(vec_laser122.size());
    for (size_t i = 0; i < vec_laser122.size(); ++i) {
        const Laser& lsr = vec_laser122[i];
        TVector3 lr = BeamToLaserCoord(r, lsr);
        double x = lr.X(), y = lr.Y();
        double sx2 = lsr.sigma_x * lsr.sigma_x;
        double sy2 = lsr.sigma_y * lsr.sigma_y;
        // E-field: exp(-(x²/σx² + y²/σy²))
        double exp_s_E = exp(-(x*x)/sx2 - (y*y)/sy2);
        double pre_E = sqrt(2.0 * sqrt(2.0) * eta * lsr.energy / (pi32 * lsr.sigma_x * lsr.sigma_y * lsr.tau * 1e-9));
        cached_Espatial_122[i] = pre_E * exp_s_E;
        // Intensity: exp(-2*(x²/σx² + y²/σy²)) — computed directly to match original formula
        double exp_s_I = exp(-2.0*(x*x)/sx2 - 2.0*(y*y)/sy2);
        double pre_I = (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * 1e-9);
        cached_Ispatial_122[i] = pre_I * exp_s_I;
    }

    cached_Ispatial_355.resize(vec_laser355.size());
    for (size_t i = 0; i < vec_laser355.size(); ++i) {
        const Laser& lsr = vec_laser355[i];
        TVector3 lr = BeamToLaserCoord(r, lsr);
        double x = lr.X(), y = lr.Y();
        double sx2 = lsr.sigma_x * lsr.sigma_x;
        double sy2 = lsr.sigma_y * lsr.sigma_y;
        double exp_s = exp(-(x*x)/sx2 - (y*y)/sy2);
        double pre_I = (2.0 * lsr.energy) / (sqrt(2.0 * TMath::Pi()) * TMath::Pi() * lsr.sigma_x * lsr.sigma_y * lsr.tau * 1e-9);
        cached_Ispatial_355[i] = pre_I * exp_s * exp_s;
    }
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

    Eigen::Vector3d laser_dirc_beforerot(0, 0, 1);
    Eigen::Vector3d after_rot = lsr.rot_mat_rev * laser_dirc_beforerot;

    lsr.laser_dirc = {after_rot(2), after_rot(0), after_rot(1)};

//    std::cout << "--- Rotation matrix updated.\n122nm:\n" << rot_mat_122 << "\n355nm:\n" << rot_mat_355 << std::endl;
//    std::cout << "The 122nm wave vector direction:\n" << laser_dirc.X() << ", " << laser_dirc.Y() << ", " << laser_dirc.Z() << std::endl;
}

void LaserGenerator::AddLaser122(Double_t energy, Double_t pulse_FWHM, Double_t peak_time, Double_t linewidth,
                                 Double_t sigma_x, Double_t sigma_y, Double_t offset_x, Double_t offset_y,
                                 Double_t offset_z, Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning) {
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
    lsr_tmp.detuning = detuning;

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
        << "  Detuning: " << detuning << " GHz\n"
//        << "  Wavevector (k): " << lsr_tmp.laser_k << " m^-1\n"
        << "  Wavevector direction (unit vector): (" << lsr_tmp.laser_dirc.X() << ", "
        << lsr_tmp.laser_dirc.Y() << ", "
        << lsr_tmp.laser_dirc.Z() << ")\n"
        << "  Rotation matrix:\n" << lsr_tmp.rot_mat << std::endl;

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
        << "  Rotation matrix:\n" << lsr_tmp.rot_mat << std::endl;

    cout << oss.str();

    vec_laser355.push_back(lsr_tmp);
}
