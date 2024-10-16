//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_LASERGENERATOR_H
#define LASERTOYMC_LASERGENERATOR_H
#include "common.h"

class LaserGenerator {
public:
    // Meyers' Singleton - Get the instance of the class
    static LaserGenerator& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    LaserGenerator(const LaserGenerator&) = delete;
    LaserGenerator& operator=(const LaserGenerator&) = delete;

    void SetEnergy(Double_t E){energy=E;};  // in J
    void SetLinewidth(Double_t l){linewidth=l;};    // in GHz
    void SetSigmaX(Double_t x){sigma_x=x;}; // in mm
    void SetSigmaY(Double_t y){sigma_y=y;}; // in mm
    void SetPulseTimeWidth(Double_t p){tau=p;}; // in ns
    void SetLaserPosition(Double_t z){laser_z=z;};  // in mm
    void SetCentralFreq(Double_t freq){cen_freq=freq; laser_k=2*TMath::Pi()*cen_freq*1e9/299792458;}; // in GHz. k in m^-1
    void SetWaveLength(Double_t wvl){laser_k=2*TMath::Pi()/wvl*1e9; cen_freq=299792458/wvl;}  // in nm
    void SetLaserDirection(TVector3 k){k_dirc=k.Unit();};
    void SetPeakTime(Double_t t){peak_time=t;};

    Double_t GetEnergy(){return energy;};
    Double_t GetLinewidth(){return linewidth;};
    Double_t GetSigmaX(){return sigma_x;};
    Double_t GetSigmaY(){return sigma_y;};
    Double_t GetPulseTimeWidth(){return tau;};
    Double_t GetLaserPosition(){return laser_z;};
    TVector3 GetWaveVector(){return laser_k*k_dirc;};   // in m^-1

    TVector3 GetFieldE(TVector3 r, Double_t t); // in V/mm

private:
    // Private constructor and destructor
    LaserGenerator(){peak_time=0;};
    ~LaserGenerator(){};


    Double_t energy;
    Double_t linewidth;
    Double_t sigma_x;
    Double_t sigma_y;
    Double_t tau;
    Double_t peak_time;

    Double_t laser_z;

    Double_t cen_freq;
    Double_t laser_k;
    TVector3 k_dirc;



};

#endif //LASERTOYMC_LASERGENERATOR_H
