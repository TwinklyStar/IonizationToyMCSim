//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_ROOTMANAGER_H
#define LASERTOYMC_ROOTMANAGER_H
#include "common.h"

class RootManager {
public:
    // Meyers' Singleton - Get the instance of the class
    static RootManager& GetInstance(const TString name="");

    // Delete copy constructor and assignment operator to avoid copying
    RootManager(const RootManager&) = delete;
    RootManager& operator=(const RootManager&) = delete;

    void SetEventID(Int_t id){eventID=id;};
    void SetPosition(TVector3 r){x=r.X(); y=r.Y(); z=r.Z();};
    void SetVelocity(TVector3 v){vx=v.X(); vy=v.Y(); vz=v.Z();};
    void SetDoppFreq(Double_t freq){dopp_freq=freq;};
    void PushTimePoint(Double_t tt, Double_t trabi_freq,
                       Double_t trho_gg, Double_t trho_ee,
                       Double_t trho_ge_r, Double_t trho_ge_i, Double_t trho_ion);
    void SetLastState(){last_rho_gg=rho_gg.back(); last_rho_ee=rho_ee.back(); last_rho_ion=rho_ion.back();};
    void FillEvent();

    void Finalize(){output_file->Write(); output_file->Close();};


private:
    // Private constructor and destructor
    ~RootManager(){};
    RootManager(TString name);

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
    Int_t step_n;
    std::vector<Double_t> t;
    std::vector<Double_t> rabi_freq;
    std::vector<Double_t> E_field;
    std::vector<Double_t> rho_gg;
    std::vector<Double_t> rho_ee;
    std::vector<Double_t> rho_ge_r;
    std::vector<Double_t> rho_ge_i;
    std::vector<Double_t> rho_ion;
    Double_t last_rho_gg;
    Double_t last_rho_ee;
    Double_t last_rho_ion;


};
#endif //LASERTOYMC_ROOTMANAGER_H
