//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_ROOTMANAGER_H
#define LASERTOYMC_ROOTMANAGER_H
#include "common.h"

class RunManager;

class RootManager {
public:
    // Meyers' Singleton - Get the instance of the class
    static RootManager& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    RootManager(const RootManager&) = delete;
    RootManager& operator=(const RootManager&) = delete;

    void Initialize();

    void SetOutFileName(std::string name) {outfile_name = name;};
    void Sett(bool flag) {IftOn=flag;};
    void SetRabiFreq(bool flag) {IfRabiFreqOn=flag;};
    void SetEField(bool flag) {IfEFieldOn=flag;};
    void SetIntensity122(bool flag) {IfIntensity122On=flag;};
    void SetIntensity355(bool flag) {IfIntensity355On=flag;};
    void SetGammaIon(bool flag) {IfGammaIonOn=flag;};
    void Setrho_gg(bool flag) {Ifrho_ggOn=flag;};
    void Setrho_ee(bool flag) {Ifrho_eeOn=flag;};
    void Setrho_ge_r(bool flag) {Ifrho_ge_rOn=flag;};
    void Setrho_ge_i(bool flag) {Ifrho_ge_iOn=flag;};
    void Setrho_ion(bool flag) {Ifrho_ionOn=flag;};

    void SetEventID(Int_t id){eventID=id;};
    void SetPosition(TVector3 r){x=r.X(); y=r.Y(); z=r.Z();};
    void SetVelocity(TVector3 v){vx=v.X(); vy=v.Y(); vz=v.Z();};
    void SetDoppFreq(Double_t shift){dopp_freq=shift/TMath::Pi()/2;};
    void SetLaserPars(Double_t E, Double_t E_355, Double_t sigmat, Double_t sigmax, Double_t sigmay, Double_t intensity,
                      Double_t intensity_355, Double_t linw);
    void PushTimePoint(Double_t tt, Double_t tE_field, Double_t tintensity_122, Double_t tintensity_355, Double_t trabi_freq,
                       Double_t trho_gg, Double_t trho_ee,
                       Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion, Double_t tgamma_ion);
    void SetLastState();
    void FillEvent();

    void Finalize(){output_file->Write(); output_file->Close();};


private:
    // Private constructor and destructor
    ~RootManager(){};
    RootManager(){outfile_name = "OBE00.root";};

    TString outfile_name;
    TFile *output_file;
    TTree *output_tree;

    Int_t eventID;
    Double_t x;
    Double_t y;
    Double_t z;
    Double_t vx;
    Double_t vy;
    Double_t vz;
    Double_t dopp_freq;
    Double_t pulse_energy;
    Double_t pulse_energy_355;
    Double_t laser_sigmat;
    Double_t laser_sigmax;
    Double_t laser_sigmay;
    Double_t peak_intensity;
    Double_t peak_intensity_355;
    Double_t linewidth;
    Int_t step_n;
    std::vector<Double_t> t;
    std::vector<Double_t> rabi_freq;
    std::vector<Double_t> gamma_ion;
    std::vector<Double_t> E_field;
    std::vector<Double_t> intensity_122;
    std::vector<Double_t> intensity_355;
    std::vector<Double_t> rho_gg;
    std::vector<Double_t> rho_ee;
    std::vector<Double_t> rho_ge_r;
    std::vector<Double_t> rho_ge_i;
    std::vector<Double_t> rho_ion;
    Double_t last_rho_gg;
    Double_t last_rho_ee;
    Double_t last_rho_ion;
    Int_t if_ionized;
    Double_t ioni_time;

    bool IftOn=true;
    bool IfRabiFreqOn=true;
    bool IfEFieldOn=true;
    bool IfIntensity122On=true;
    bool IfIntensity355On=true;
    bool IfGammaIonOn=true;
    bool Ifrho_ggOn=true;
    bool Ifrho_eeOn=true;
    bool Ifrho_ge_rOn=true;
    bool Ifrho_ge_iOn=true;
    bool Ifrho_ionOn=true;


};
#endif //LASERTOYMC_ROOTMANAGER_H
