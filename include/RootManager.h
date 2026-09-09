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
    void SetLaserPars122(bool flag) {IfLaserPars122On=flag;};
    void SetLaserPars355(bool flag) {IfLaserPars355On=flag;};
    bool IsLaserPars122On() const {return IfLaserPars122On;};
    bool IsLaserPars355On() const {return IfLaserPars355On;};

    void SetEventID(Int_t id){eventID=id;};
    void SetPosition(TVector3 r){x=r.X(); y=r.Y(); z=r.Z();};
    void SetVelocity(TVector3 v){vx=v.X(); vy=v.Y(); vz=v.Z();};
    void SetDoppFreq(Double_t shift){dopp_freq=shift/TMath::Pi()/2;};
    void SetLaserPars(Double_t E, Double_t E_355, Double_t sigmat, Double_t sigmax, Double_t sigmay, Double_t intensity,
                      Double_t intensity_355, Double_t linw);
    void SetLaserPars(Double_t intensity, Double_t intensity_355);
    void PushTimePoint(Double_t tt, Double_t tE_field, Double_t tintensity_122, Double_t tintensity_355, Double_t trabi_freq,
                       Double_t trho_gg, Double_t trho_ee,
                       Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion, Double_t tgamma_ion);
    // Record the single 122nm/355nm laser's realized (nominal or jittered) parameter values
    // for the current event. Call once per event, before FillEvent().
    void SetLaser122Snapshot(Double_t energy, Double_t linewidth, Double_t peak_time,
                             Double_t sigma_x, Double_t sigma_y, Double_t tau,
                             Double_t offset_x, Double_t offset_y, Double_t offset_z,
                             Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning);
    void SetLaser355Snapshot(Double_t energy, Double_t linewidth, Double_t peak_time,
                             Double_t sigma_x, Double_t sigma_y, Double_t tau,
                             Double_t offset_x, Double_t offset_y, Double_t offset_z,
                             Double_t yaw, Double_t pitch, Double_t roll);
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
    Double_t peak_intensity=0;
    Double_t peak_intensity_355=0;
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
    // Per-event realized laser parameter snapshot for the single 122nm/355nm laser, written when
    // IfLaserPars122On/IfLaserPars355On is set. Populated via SetLaser122Snapshot/SetLaser355Snapshot.
    // Only one laser of each wavelength is supported (see the guard in RunManager::SolveOBE).
    Double_t laser122_energy, laser122_linewidth, laser122_peak_time;
    Double_t laser122_sigma_x, laser122_sigma_y, laser122_tau;
    Double_t laser122_offset_x, laser122_offset_y, laser122_offset_z;
    Double_t laser122_yaw, laser122_pitch, laser122_roll, laser122_detuning;
    Double_t laser355_energy, laser355_linewidth, laser355_peak_time;
    Double_t laser355_sigma_x, laser355_sigma_y, laser355_tau;
    Double_t laser355_offset_x, laser355_offset_y, laser355_offset_z;
    Double_t laser355_yaw, laser355_pitch, laser355_roll;
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
    bool IfLaserPars122On=false;
    bool IfLaserPars355On=false;


};
#endif //LASERTOYMC_ROOTMANAGER_H
