//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_LASERGENERATOR_H
#define LASERTOYMC_LASERGENERATOR_H
#include "common.h"

class OBEsolver;    // Two classes cannot include each other

class LaserGenerator {
public:
    // Meyers' Singleton - Get the instance of the class
    static LaserGenerator& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    LaserGenerator(const LaserGenerator&) = delete;
    LaserGenerator& operator=(const LaserGenerator&) = delete;

    void SetEnergy(Double_t E){energy=E;};  // in J
    void SetEnergy355(Double_t E){energy_355=E;};   // in J
    void SetLinewidth(Double_t l){linewidth=l;};    // in GHz
    void SetSigmaX(Double_t x){sigma_x=x;}; // in mm
    void SetSigmaY(Double_t y){sigma_y=y;}; // in mm
    void SetPulseTimeWidth(Double_t p){tau=p;}; // in ns
    void SetPulseFWHM(Double_t p){tau=0.4247*p;};   // in ns
    void SetLaserOffset(TVector3 r){laser_offset=r;};  // in mm
    void SetLaserOffset355(TVector3 r){laser_offset_355=r;};  // in mm
    void SetYawAngle(Double_t y){yaw = y * TMath::Pi()/180;};     // in deg
    void SetPitchAngle(Double_t p){pitch = p * TMath::Pi()/180;}; // in deg
    void SetRollAngle(Double_t r){roll = r * TMath::Pi()/180;};   // in deg
    void SetYawAngle355(Double_t y){yaw_355 = y * TMath::Pi()/180;};     // in deg
    void SetPitchAngle355(Double_t p){pitch_355 = p * TMath::Pi()/180;}; // in deg
    void SetRollAngle355(Double_t r){roll_355 = r * TMath::Pi()/180;};   // in deg
    void SetCentralFreq(Double_t freq){cen_freq=freq; laser_k=2*TMath::Pi()*cen_freq*1e9/299792458;}; // in GHz. k in m^-1
    void SetWaveLength(Double_t wvl){laser_k=2*TMath::Pi()/wvl*1e9; cen_freq=299792458/wvl;}  // in nm
    void SetPeakTime(Double_t t){peak_time=t;};

    void UpdateRotMat();

    void SetOBESolverPtr(OBEsolver *ptr){obe_ptr=ptr;};

    Double_t GetEnergy(){return energy;};
    Double_t GetEnergy355(){return energy_355;};
    Double_t GetLinewidth(){return linewidth;};
    Double_t GetSigmaX(){return sigma_x;};
    Double_t GetSigmaY(){return sigma_y;};
    Double_t GetPulseTimeWidth(){return tau;};
    Double_t GetPeakIntensity(TVector3 r);      // in W/cm^2
    Double_t GetPeakIntensity355(TVector3 r);   // in W/cm^2
    TVector3 GetWaveVector();   // in m^-1

    TVector3 GetFieldE(TVector3 r, Double_t t); // in V/mm
    Double_t GetIntensity(TVector3 r, Double_t t);    // in W/cm^2
    Double_t GetIntensity355(TVector3 r, Double_t t);    // in W/cm^2

private:
    // Private constructor and destructor
    LaserGenerator();
    ~LaserGenerator(){};

    TVector3 BeamToLaserCoord(TVector3 r); // Transform from target coordinate to laser coordinate
    TVector3 BeamToLaserCoord355(TVector3 r); // Transform from target coordinate to laser coordinate

    OBEsolver *obe_ptr;

    Double_t energy;
    Double_t energy_355;
    Double_t linewidth;
    Double_t sigma_x;
    Double_t sigma_y;
    Double_t tau;
    Double_t peak_time;
    Double_t yaw;
    Double_t pitch;
    Double_t roll;
    Double_t yaw_355;
    Double_t pitch_355;
    Double_t roll_355;

    TVector3 laser_offset;
    TVector3 laser_offset_355;

    Double_t cen_freq;
    Double_t laser_k;
    TVector3 laser_dirc;

    Eigen::Matrix3d rot_mat_122;
    Eigen::Matrix3d rot_mat_355;
    Eigen::Matrix3d rot_mat_rev_122;

};

#endif //LASERTOYMC_LASERGENERATOR_H
