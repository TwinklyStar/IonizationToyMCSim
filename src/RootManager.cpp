//
// Created by Meng Lv on 2024/10/9.
//

#include "RootManager.h"
#include "RunManager.h"

void RootManager::Initialize() {
    std::cout << "--- Output file name: " << outfile_name << std::endl;
    output_file = TFile::Open(outfile_name, "RECREATE");
    if (!output_file || !output_file->IsOpen())
        throw std::runtime_error("RootManager::Initialize: failed to open output file: \"" + outfile_name + "\"");
    output_tree = new TTree("obe", "A tree storing parameters and time evolution of states");

    output_tree->Branch("EventID", &eventID);
    output_tree->Branch("x", &x);
    output_tree->Branch("y", &y);
    output_tree->Branch("z", &z);
    output_tree->Branch("vx", &vx);
    output_tree->Branch("vy", &vy);
    output_tree->Branch("vz", &vz);
//    output_tree->Branch("PulseEnergy", &pulse_energy);
//    output_tree->Branch("PulseEnergy355", &pulse_energy_355);
//    output_tree->Branch("LineWidth", &linewidth);
//    output_tree->Branch("LaserSigmaT", &laser_sigmat);
//    output_tree->Branch("LaserSigmaX", &laser_sigmax);
//    output_tree->Branch("LaserSigmaY", &laser_sigmay);
//    output_tree->Branch("PeakIntensity", &peak_intensity);
//    output_tree->Branch("PeakIntensity355", &peak_intensity_355);
    output_tree->Branch("DoppFreq", &dopp_freq);
    output_tree->Branch("PeakIntensity122", &peak_intensity);
    output_tree->Branch("PeakIntensity355", &peak_intensity_355);
    output_tree->Branch("Step_n", &step_n);
    if (IftOn)          output_tree->Branch("t", &t);
    if (IfRabiFreqOn)   output_tree->Branch("RabiFreq", &rabi_freq);
    if (IfEFieldOn)     output_tree->Branch("EField", &E_field);
    if (IfIntensity122On)output_tree->Branch("Intensity122", &intensity_122);
    if (IfIntensity355On)output_tree->Branch("Intensity355", &intensity_355);
    if (IfGammaIonOn)   output_tree->Branch("GammaIon", &gamma_ion);
    if (Ifrho_ggOn)     output_tree->Branch("rho_gg", &rho_gg);
    if (Ifrho_eeOn)     output_tree->Branch("rho_ee", &rho_ee);
    if (Ifrho_ge_rOn)   output_tree->Branch("rho_ge_r", &rho_ge_r);
    if (Ifrho_ge_iOn)   output_tree->Branch("rho_ge_i", &rho_ge_i);
    if (Ifrho_ionOn)    output_tree->Branch("rho_ion", &rho_ion);
    if (IfLaserPars122On) {
        output_tree->Branch("Laser122_Energy", &laser122_energy);
        output_tree->Branch("Laser122_Linewidth", &laser122_linewidth);
        output_tree->Branch("Laser122_PeakTime", &laser122_peak_time);
        output_tree->Branch("Laser122_SigmaX", &laser122_sigma_x);
        output_tree->Branch("Laser122_SigmaY", &laser122_sigma_y);
        output_tree->Branch("Laser122_Tau", &laser122_tau);
        output_tree->Branch("Laser122_OffsetX", &laser122_offset_x);
        output_tree->Branch("Laser122_OffsetY", &laser122_offset_y);
        output_tree->Branch("Laser122_OffsetZ", &laser122_offset_z);
        output_tree->Branch("Laser122_Yaw", &laser122_yaw);
        output_tree->Branch("Laser122_Pitch", &laser122_pitch);
        output_tree->Branch("Laser122_Roll", &laser122_roll);
        output_tree->Branch("Laser122_Detuning", &laser122_detuning);
    }
    if (IfLaserPars355On) {
        output_tree->Branch("Laser355_Energy", &laser355_energy);
        output_tree->Branch("Laser355_Linewidth", &laser355_linewidth);
        output_tree->Branch("Laser355_PeakTime", &laser355_peak_time);
        output_tree->Branch("Laser355_SigmaX", &laser355_sigma_x);
        output_tree->Branch("Laser355_SigmaY", &laser355_sigma_y);
        output_tree->Branch("Laser355_Tau", &laser355_tau);
        output_tree->Branch("Laser355_OffsetX", &laser355_offset_x);
        output_tree->Branch("Laser355_OffsetY", &laser355_offset_y);
        output_tree->Branch("Laser355_OffsetZ", &laser355_offset_z);
        output_tree->Branch("Laser355_Yaw", &laser355_yaw);
        output_tree->Branch("Laser355_Pitch", &laser355_pitch);
        output_tree->Branch("Laser355_Roll", &laser355_roll);
    }
    output_tree->Branch("LastRho_gg", &last_rho_gg);
    output_tree->Branch("LastRho_ee", &last_rho_ee);
    output_tree->Branch("LastRho_ion", &last_rho_ion);
    output_tree->Branch("IfIonized", &if_ionized);
    output_tree->Branch("IoniTime", &ioni_time);
}

void RootManager::PushTimePoint(Double_t tt, Double_t tE_field, Double_t tintensity_122, Double_t tintensity_355,
                                Double_t trabi_freq, Double_t trho_gg, Double_t trho_ee,
                                Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion, Double_t tgamma_ion) {
    t.push_back(tt);
    E_field.push_back(tE_field);
    intensity_122.push_back(tintensity_122);
    peak_intensity = TMath::Max(peak_intensity, tintensity_122);
    intensity_355.push_back(tintensity_355);
    peak_intensity_355 = TMath::Max(peak_intensity_355, tintensity_355);
    rabi_freq.push_back(trabi_freq);
    rho_gg.push_back(trho_gg);
    rho_ee.push_back(trho_ee);
    rho_ge_r.push_back(trho_ge_r);
    rho_ge_i.push_back(trho_ge_i);
    rho_ion.push_back(trho_ion);
    gamma_ion.push_back(tgamma_ion);
}

void RootManager::SetLaser122Snapshot(Double_t energy, Double_t linewidth, Double_t peak_time,
                                      Double_t sigma_x, Double_t sigma_y, Double_t tau,
                                      Double_t offset_x, Double_t offset_y, Double_t offset_z,
                                      Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning) {
    laser122_energy = energy;
    laser122_linewidth = linewidth;
    laser122_peak_time = peak_time;
    laser122_sigma_x = sigma_x;
    laser122_sigma_y = sigma_y;
    laser122_tau = tau;
    laser122_offset_x = offset_x;
    laser122_offset_y = offset_y;
    laser122_offset_z = offset_z;
    laser122_yaw = yaw;
    laser122_pitch = pitch;
    laser122_roll = roll;
    laser122_detuning = detuning;
}

void RootManager::SetLaser355Snapshot(Double_t energy, Double_t linewidth, Double_t peak_time,
                                      Double_t sigma_x, Double_t sigma_y, Double_t tau,
                                      Double_t offset_x, Double_t offset_y, Double_t offset_z,
                                      Double_t yaw, Double_t pitch, Double_t roll) {
    laser355_energy = energy;
    laser355_linewidth = linewidth;
    laser355_peak_time = peak_time;
    laser355_sigma_x = sigma_x;
    laser355_sigma_y = sigma_y;
    laser355_tau = tau;
    laser355_offset_x = offset_x;
    laser355_offset_y = offset_y;
    laser355_offset_z = offset_z;
    laser355_yaw = yaw;
    laser355_pitch = pitch;
    laser355_roll = roll;
}

void RootManager::SetLaserPars(Double_t E, Double_t E_355, Double_t sigmat, Double_t sigmax, Double_t sigmay, Double_t intensity,
                              Double_t intensity_355, Double_t linw) {
    pulse_energy = E;
    pulse_energy_355 = E_355;
    laser_sigmat = sigmat;
    laser_sigmax = sigmax;
    laser_sigmay = sigmay;
    peak_intensity = intensity;
    peak_intensity_355 = intensity_355;
    linewidth = linw;
}

void RootManager::SetLaserPars(Double_t intensity, Double_t intensity_355) {
    peak_intensity = intensity;
    peak_intensity_355 = intensity_355;
}

void RootManager::SetLastState() {
    if (t.empty())
        throw std::runtime_error("RootManager::SetLastState: no time steps recorded for this event — ODE produced no output");
    last_rho_gg=rho_gg.back();
    last_rho_ee=rho_ee.back();
    last_rho_ion=rho_ion.back();

    TGraph gtemp(t.size(), rho_ion.data(), t.data());

    RunManager &RM = RunManager::GetInstance();
    Double_t uni_0to1 = RM.rdm_gen.Uniform();

    if(uni_0to1<last_rho_ion){
        if_ionized = 1;
        ioni_time = gtemp.Eval(uni_0to1);
    }
    else {
        if_ionized = 0;
        ioni_time = -1;
    }
}

void RootManager::FillEvent() {
    step_n = t.size();
    output_tree->Fill();
    t.clear();
    rabi_freq.clear();
    E_field.clear();
    intensity_122.clear();
    intensity_355.clear();
    rho_gg.clear();
    rho_ee.clear();
    rho_ge_r.clear();
    rho_ge_i.clear();
    rho_ion.clear();
    gamma_ion.clear();

    peak_intensity = 0;
    peak_intensity_355 = 0;

    // shrink_to_fit() removed: clear() retains capacity for reuse next event,
    // avoiding a malloc/free roundtrip on every event.
}

RootManager& RootManager::GetInstance() {
    static RootManager instance;
    return instance;
}
