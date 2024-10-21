//
// Created by Meng Lv on 2024/10/9.
//

#include "RootManager.h"

RootManager::RootManager(TString name): outfile_name(name){
    std::cout << "--- Output file name: " << outfile_name << std::endl;
    output_file = TFile::Open(name, "RECREATE");
    output_tree = new TTree("obe", "A tree storing parameters and time evolution of states");

    output_tree->Branch("EventID", &eventID);
    output_tree->Branch("x", &x);
    output_tree->Branch("y", &y);
    output_tree->Branch("z", &z);
    output_tree->Branch("vx", &vx);
    output_tree->Branch("vy", &vy);
    output_tree->Branch("vz", &vz);
    output_tree->Branch("DoppFreq", &dopp_freq);
    output_tree->Branch("PeakE", &peak_E);
    output_tree->Branch("Step_n", &step_n);
    output_tree->Branch("t", &t);
    output_tree->Branch("RabiFreq", &rabi_freq);
    output_tree->Branch("EField", &E_field);
    output_tree->Branch("rho_gg", &rho_gg);
    output_tree->Branch("rho_ee", &rho_ee);
    output_tree->Branch("rho_ge_r", &rho_ge_r);
    output_tree->Branch("rho_ge_i", &rho_ge_i);
    output_tree->Branch("rho_ion", &rho_ion);
    output_tree->Branch("LastRho_gg", &last_rho_gg);
    output_tree->Branch("LastRho_ee", &last_rho_ee);
    output_tree->Branch("LastRho_ion", &last_rho_ion);
};

void RootManager::PushTimePoint(Double_t tt, Double_t tE_field, Double_t trabi_freq, Double_t trho_gg, Double_t trho_ee,
                                Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion) {
    t.push_back(tt);
    E_field.push_back(tE_field);
    rabi_freq.push_back(trabi_freq);
    rho_gg.push_back(trho_gg);
    rho_ee.push_back(trho_ee);
    rho_ge_r.push_back(trho_ge_r);
    rho_ge_i.push_back(trho_ge_i);
    rho_ion.push_back(trho_ion);
}

void RootManager::FillEvent() {
    step_n = t.size();
    peak_E = *std::max_element(E_field.begin(), E_field.end());
    output_tree->Fill();
    t.clear();
    rabi_freq.clear();
    E_field.clear();
    rho_gg.clear();
    rho_ee.clear();
    rho_ge_r.clear();
    rho_ge_i.clear();
    rho_ion.clear();

    rabi_freq.shrink_to_fit();
    E_field.shrink_to_fit();
    rho_gg.shrink_to_fit();
    rho_ee.shrink_to_fit();
    rho_ge_r.shrink_to_fit();
    rho_ge_i.shrink_to_fit();
}

RootManager& RootManager::GetInstance(const TString name) {
    static RootManager instance(name);
    return instance;
}
