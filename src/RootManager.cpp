//
// Created by Meng Lv on 2024/10/9.
//

#include "RootManager.h"

void RootManager::Initialize() {
    std::cout << "--- Output file name: " << outfile_name << std::endl;
    output_file = TFile::Open(outfile_name, "RECREATE");
    output_tree = new TTree("obe", "A tree storing parameters and time evolution of states");

    output_tree->Branch("EventID", &eventID);
    output_tree->Branch("x", &x);
    output_tree->Branch("y", &y);
    output_tree->Branch("z", &z);
    output_tree->Branch("vx", &vx);
    output_tree->Branch("vy", &vy);
    output_tree->Branch("vz", &vz);
    output_tree->Branch("PulseEnergy", &pulse_energy);
    output_tree->Branch("PulseEnergy355", &pulse_energy_355);
    output_tree->Branch("LineWidth", &linewidth);
    output_tree->Branch("LaserSigmaT", &laser_sigmat);
    output_tree->Branch("LaserSigmaX", &laser_sigmax);
    output_tree->Branch("LaserSigmaY", &laser_sigmay);
    output_tree->Branch("PeakIntensity", &peak_intensity);
    output_tree->Branch("PeakIntensity355", &peak_intensity_355);
    output_tree->Branch("DoppFreq", &dopp_freq);
    output_tree->Branch("Step_n", &step_n);
    if (IftOn)          output_tree->Branch("t", &t);
    if (IfRabiFreqOn)   output_tree->Branch("RabiFreq", &rabi_freq);
    if (IfEFieldOn)     output_tree->Branch("EField", &E_field);
    if (IfGammaIonOn)   output_tree->Branch("GammaIon", &gamma_ion);
    if (Ifrho_ggOn)     output_tree->Branch("rho_gg", &rho_gg);
    if (Ifrho_eeOn)     output_tree->Branch("rho_ee", &rho_ee);
    if (Ifrho_ge_rOn)   output_tree->Branch("rho_ge_r", &rho_ge_r);
    if (Ifrho_ge_iOn)   output_tree->Branch("rho_ge_i", &rho_ge_i);
    if (Ifrho_ionOn)    output_tree->Branch("rho_ion", &rho_ion);
    output_tree->Branch("LastRho_gg", &last_rho_gg);
    output_tree->Branch("LastRho_ee", &last_rho_ee);
    output_tree->Branch("LastRho_ion", &last_rho_ion);
    output_tree->Branch("IfIonized", &if_ionized);
}

void RootManager::PushTimePoint(Double_t tt, Double_t tE_field, Double_t trabi_freq, Double_t trho_gg, Double_t trho_ee,
                                Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion, Double_t tgamma_ion) {
    t.push_back(tt);
    E_field.push_back(tE_field);
    rabi_freq.push_back(trabi_freq);
    rho_gg.push_back(trho_gg);
    rho_ee.push_back(trho_ee);
    rho_ge_r.push_back(trho_ge_r);
    rho_ge_i.push_back(trho_ge_i);
    rho_ion.push_back(trho_ion);
    gamma_ion.push_back(tgamma_ion);
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

void RootManager::SetLastState() {
    last_rho_gg=rho_gg.back();
    last_rho_ee=rho_ee.back();
    last_rho_ion=rho_ion.back();

    if(randGen.Uniform()<last_rho_ion)
        if_ionized = 1;
    else
        if_ionized = 0;
}

void RootManager::FillEvent() {
    step_n = t.size();
    output_tree->Fill();
    t.clear();
    rabi_freq.clear();
    E_field.clear();
    rho_gg.clear();
    rho_ee.clear();
    rho_ge_r.clear();
    rho_ge_i.clear();
    rho_ion.clear();
    gamma_ion.clear();

    rabi_freq.shrink_to_fit();
    E_field.shrink_to_fit();
    rho_gg.shrink_to_fit();
    rho_ee.shrink_to_fit();
    rho_ge_r.shrink_to_fit();
    rho_ge_i.shrink_to_fit();
    gamma_ion.shrink_to_fit();
}

RootManager& RootManager::GetInstance() {
    static RootManager instance;
    return instance;
}
