//
// Created by Meng Lv on 2024/10/9.
//

#include "LaserGenerator.h"
#include "OBEsolver.h"
#include "RunManager.h"

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

        if (lsr.profile) {
            sum_I += (lsr.energy / (sqrt(2.0 * TMath::Pi()) * lsr.tau * 1e-9))
                     * SampleProfileDensity(lsr.profile.get(), x, y) * 100;
            continue;
        }
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

        if (lsr.profile) {
            sum_I += (lsr.energy / (sqrt(2.0 * TMath::Pi()) * lsr.tau * 1e-9))
                     * SampleProfileDensity(lsr.profile.get(), x, y) * 100;
            continue;
        }
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
        if (lsr.profile) {
            // Measured transverse density h(x,y) [mm^-2] replaces the analytic Gaussian.
            // pre_I * exp_s_I  ==  [E / (sqrt(2*pi) * tau[s])] * h_analytic(x,y), so the
            // normalized measured h drops straight in and sigma_x/sigma_y no longer enter.
            double hval = SampleProfileDensity(lsr.profile.get(), x, y);
            double Ispatial = (lsr.energy / (sqrt(2.0 * TMath::Pi()) * lsr.tau * 1e-9)) * hval;
            cached_Ispatial_122[i] = Ispatial;
            // Analytic path keeps cached_Espatial_122^2 / cached_Ispatial_122 == 2*eta exactly;
            // reuse that identity so the E-field stays consistent with the intensity.
            cached_Espatial_122[i] = sqrt(2.0 * eta * Ispatial);
            continue;
        }
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
        if (lsr.profile) {
            double hval = SampleProfileDensity(lsr.profile.get(), x, y);
            cached_Ispatial_355[i] = (lsr.energy / (sqrt(2.0 * TMath::Pi()) * lsr.tau * 1e-9)) * hval;
            continue;
        }
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
    vec_laser122_nominal.push_back(lsr_tmp);
    vec_laser122_sigma.emplace_back();
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
    vec_laser355_nominal.push_back(lsr_tmp);
    vec_laser355_sigma.emplace_back();
}

void LaserGenerator::SetLaser122Sigma(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                                      Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                                      Double_t offset_x, Double_t offset_y, Double_t offset_z,
                                      Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning) {
    if (vec_laser122_sigma.empty())
        throw std::runtime_error("LaserGenerator::SetLaser122Sigma: no 122nm laser added yet — call AddLaser122 first");

    LaserSigma &sig = vec_laser122_sigma.back();
    sig.energy = energy;
    sig.linewidth = linewidth;
    sig.peak_time = peak_time;
    sig.sigma_x = sigma_x;
    sig.sigma_y = sigma_y;
    sig.tau = 0.4247 * pulse_FWHM;
    sig.offset_x = offset_x;
    sig.offset_y = offset_y;
    sig.offset_z = offset_z;
    sig.yaw = yaw;
    sig.pitch = pitch;
    sig.roll = roll;
    sig.detuning = detuning;

    std::cout << "-- LaserGenerator: Set sigma for the last 122nm laser (energy=" << energy
              << " J, FWHM=" << pulse_FWHM << " ns, peak_time=" << peak_time << " ns, linewidth=" << linewidth
              << " GHz, sigma_x=" << sigma_x << " mm, sigma_y=" << sigma_y << " mm, offset=(" << offset_x
              << ", " << offset_y << ", " << offset_z << ") mm, yaw=" << yaw << " deg, pitch=" << pitch
              << " deg, roll=" << roll << " deg, detuning=" << detuning << " GHz)" << std::endl;
}

void LaserGenerator::SetLaser355Sigma(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                                      Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                                      Double_t offset_x, Double_t offset_y, Double_t offset_z,
                                      Double_t yaw, Double_t pitch, Double_t roll) {
    if (vec_laser355_sigma.empty())
        throw std::runtime_error("LaserGenerator::SetLaser355Sigma: no 355nm laser added yet — call AddLaser355 first");

    LaserSigma &sig = vec_laser355_sigma.back();
    sig.energy = energy;
    sig.linewidth = linewidth;
    sig.peak_time = peak_time;
    sig.sigma_x = sigma_x;
    sig.sigma_y = sigma_y;
    sig.tau = 0.4247 * pulse_FWHM;
    sig.offset_x = offset_x;
    sig.offset_y = offset_y;
    sig.offset_z = offset_z;
    sig.yaw = yaw;
    sig.pitch = pitch;
    sig.roll = roll;

    std::cout << "-- LaserGenerator: Set sigma for the last 355nm laser (energy=" << energy
              << " J, FWHM=" << pulse_FWHM << " ns, peak_time=" << peak_time << " ns, linewidth=" << linewidth
              << " GHz, sigma_x=" << sigma_x << " mm, sigma_y=" << sigma_y << " mm, offset=(" << offset_x
              << ", " << offset_y << ", " << offset_z << ") mm, yaw=" << yaw << " deg, pitch=" << pitch
              << " deg, roll=" << roll << " deg)" << std::endl;
}

void LaserGenerator::LoadProfileInto(std::vector<Laser> &live, std::vector<Laser> &nominal,
                                    const std::string &path, const char *tag) {
    if (live.empty())
        throw std::runtime_error(std::string("LaserGenerator::SetLaser") + tag +
                                 "Profile: no " + tag + " laser added yet — call AddLaser" + tag + " first");

    TFile f(path.c_str(), "READ");
    if (f.IsZombie())
        throw std::runtime_error(std::string("LaserGenerator::LoadProfileInto: cannot open profile file \"") +
                                 path + "\"");
    auto *src = dynamic_cast<TH2D *>(f.Get("h_profile"));
    if (!src)
        throw std::runtime_error(std::string("LaserGenerator::LoadProfileInto: no TH2D named \"h_profile\" in \"") +
                                 path + "\"");

    std::shared_ptr<TH2D> h(static_cast<TH2D *>(src->Clone()));
    h->SetDirectory(nullptr);
    f.Close();

    live.back().profile = h;
    nominal.back().profile = h;   // shared ownership; jitter keeps the same lookup

    const TAxis *ax = h->GetXaxis();
    const TAxis *ay = h->GetYaxis();
    std::cout << "-- LaserGenerator: Loaded real transverse profile for the last " << tag << " laser\n"
              << "     file:      " << path << "\n"
              << "     bins:      " << h->GetNbinsX() << " x " << h->GetNbinsY() << "\n"
              << "     x range:   [" << ax->GetXmin() << ", " << ax->GetXmax() << "] mm (sigma_x direction)\n"
              << "     y range:   [" << ay->GetXmin() << ", " << ay->GetXmax() << "] mm (sigma_y direction)\n"
              << "     integral:  " << h->Integral("width") << " (should be ~1)" << std::endl;
    if (live.back().sigma_x != 0 || live.back().sigma_y != 0)
        std::cout << "   WARNING LaserGenerator: sigma_x/sigma_y are ignored for a " << tag
                  << " laser with a real profile" << std::endl;
}

Double_t LaserGenerator::SampleProfileDensity(const TH2D *h, Double_t x, Double_t y) const {
    const TAxis *ax = h->GetXaxis();
    const TAxis *ay = h->GetYaxis();
    // Stay strictly inside the outermost bin centres: TH2::Interpolate needs a
    // surrounding 2x2 block of bin centres and otherwise prints an Error and returns 0.
    if (x <= ax->GetBinCenter(1) || x >= ax->GetBinCenter(ax->GetNbins()) ||
        y <= ay->GetBinCenter(1) || y >= ay->GetBinCenter(ay->GetNbins()))
        return 0.0;
    return h->Interpolate(x, y);
}

void LaserGenerator::SetLaser122Profile(const std::string &path) {
    LoadProfileInto(vec_laser122, vec_laser122_nominal, path, "122");
}

void LaserGenerator::SetLaser355Profile(const std::string &path) {
    LoadProfileInto(vec_laser355, vec_laser355_nominal, path, "355");
}

void LaserGenerator::ResampleOneLaser(Laser &live, const Laser &nominal, const LaserSigma &sigma, bool has_detuning) {
    RunManager &RM = RunManager::GetInstance();
    const int kMaxRetries = 1000;

    auto sampleNonNegative = [&](Double_t mean, Double_t sd, const char *name) {
        Double_t v = RM.rdm_gen.Gaus(mean, sd);
        for (int attempt = 0; v < 0 && attempt < kMaxRetries; ++attempt) v = RM.rdm_gen.Gaus(mean, sd);
        if (v < 0)
            throw std::runtime_error(std::string("LaserGenerator::ResampleOneLaser: could not sample a non-negative ") +
                                     name + " (mean=" + std::to_string(mean) + ", sigma=" + std::to_string(sd) +
                                     ") after " + std::to_string(kMaxRetries) + " retries — check the configured sigma");
        return v;
    };
    auto samplePositive = [&](Double_t mean, Double_t sd, const char *name) {
        Double_t v = RM.rdm_gen.Gaus(mean, sd);
        for (int attempt = 0; v <= 0 && attempt < kMaxRetries; ++attempt) v = RM.rdm_gen.Gaus(mean, sd);
        if (v <= 0)
            throw std::runtime_error(std::string("LaserGenerator::ResampleOneLaser: could not sample a positive ") +
                                     name + " (mean=" + std::to_string(mean) + ", sigma=" + std::to_string(sd) +
                                     ") after " + std::to_string(kMaxRetries) + " retries — check the configured sigma");
        return v;
    };

    live.energy = sampleNonNegative(nominal.energy, sigma.energy, "energy");
    live.linewidth = sampleNonNegative(nominal.linewidth, sigma.linewidth, "linewidth");
    live.peak_time = RM.rdm_gen.Gaus(nominal.peak_time, sigma.peak_time);
    // sigma_x/sigma_y are unused when a real profile is loaded — skip them so a
    // 0 placeholder in the macro does not trip samplePositive's retry/throw.
    if (!live.profile) {
        live.sigma_x = samplePositive(nominal.sigma_x, sigma.sigma_x, "sigma_x");
        live.sigma_y = samplePositive(nominal.sigma_y, sigma.sigma_y, "sigma_y");
    }
    live.tau = samplePositive(nominal.tau, sigma.tau, "tau");
    live.laser_offset = {RM.rdm_gen.Gaus(nominal.laser_offset.X(), sigma.offset_x),
                         RM.rdm_gen.Gaus(nominal.laser_offset.Y(), sigma.offset_y),
                         RM.rdm_gen.Gaus(nominal.laser_offset.Z(), sigma.offset_z)};
    live.yaw = RM.rdm_gen.Gaus(nominal.yaw, sigma.yaw * TMath::Pi() / 180);
    live.pitch = RM.rdm_gen.Gaus(nominal.pitch, sigma.pitch * TMath::Pi() / 180);
    live.roll = RM.rdm_gen.Gaus(nominal.roll, sigma.roll * TMath::Pi() / 180);
    if (has_detuning) live.detuning = RM.rdm_gen.Gaus(nominal.detuning, sigma.detuning);

    UpdateRotMat(live);
}

void LaserGenerator::ResampleLaserPars() {
    if (!jitter_on) return;

    for (size_t i = 0; i < vec_laser122.size(); ++i)
        ResampleOneLaser(vec_laser122[i], vec_laser122_nominal[i], vec_laser122_sigma[i], /*has_detuning=*/true);
    for (size_t i = 0; i < vec_laser355.size(); ++i)
        ResampleOneLaser(vec_laser355[i], vec_laser355_nominal[i], vec_laser355_sigma[i], /*has_detuning=*/false);
}
