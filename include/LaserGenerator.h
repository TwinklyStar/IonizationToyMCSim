//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_LASERGENERATOR_H
#define LASERTOYMC_LASERGENERATOR_H
#include "common.h"
#include <memory>
#include "TH2D.h"

class OBEsolver;    // Two classes cannot include each other

class LaserGenerator {
public:
    struct Laser {
        Double_t energy;        // in J
        Double_t linewidth;     // in GHz
        Double_t peak_time;     // in ns
        Double_t sigma_x;       // in mm
        Double_t sigma_y;       // in mm
        Double_t tau;           // in ns
        TVector3 laser_offset;  // in mm
        Double_t yaw;           // in deg
        Double_t pitch;         // in deg
        Double_t roll;          // in deg
        Double_t wavelength;    // in nm

        Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d rot_mat_rev = Eigen::Matrix3d::Identity();

        Double_t cen_freq;
        Double_t detuning;
        Double_t laser_k;       // in m^-1
        TVector3 laser_dirc;

        // Measured transverse-intensity density h(x,y) [mm^-2] on the laser frame
        // (x = sigma_x direction, y = sigma_y direction), centred on its own
        // centroid and normalized to unit integral. nullptr => analytic Gaussian.
        // When set, sigma_x/sigma_y are ignored for this laser.
        std::shared_ptr<TH2D> profile = nullptr;
    };

    // Standard deviations for per-event Gaussian randomization of a Laser's parameters.
    // A field left at 0 (the default) keeps that parameter fixed at its nominal value.
    struct LaserSigma {
        Double_t energy = 0;        // in J
        Double_t linewidth = 0;     // in GHz
        Double_t peak_time = 0;     // in ns
        Double_t sigma_x = 0;       // in mm
        Double_t sigma_y = 0;       // in mm
        Double_t tau = 0;           // in ns (= 0.4247 * pulse_FWHM sigma)
        Double_t offset_x = 0;      // in mm
        Double_t offset_y = 0;      // in mm
        Double_t offset_z = 0;      // in mm
        Double_t yaw = 0;           // in deg
        Double_t pitch = 0;         // in deg
        Double_t roll = 0;          // in deg
        Double_t detuning = 0;      // in GHz, 122nm only
    };

public:
    // Meyers' Singleton - Get the instance of the class
    static LaserGenerator& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    LaserGenerator(const LaserGenerator&) = delete;
    LaserGenerator& operator=(const LaserGenerator&) = delete;

    void AddLaser122(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                     Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                     Double_t offset_x, Double_t offset_y, Double_t offset_z,
                     Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning);
    void AddLaser355(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                     Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                     Double_t offset_x, Double_t offset_y, Double_t offset_z,
                     Double_t yaw, Double_t pitch, Double_t roll);
    void SetEnergy(Double_t E){for (auto &itr : vec_laser122) itr.energy=E;};  // in J
    void SetEnergy355(Double_t E){for (auto &itr : vec_laser355) itr.energy=E;};   // in J

    // Set the Gaussian sigma of each parameter for the most-recently-added 122/355nm laser.
    // Same parameter order/units as AddLaser122/AddLaser355.
    void SetLaser122Sigma(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                          Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                          Double_t offset_x, Double_t offset_y, Double_t offset_z,
                          Double_t yaw, Double_t pitch, Double_t roll, Double_t detuning);
    void SetLaser355Sigma(Double_t energy, Double_t pulse_FWHM, Double_t peak_time,
                          Double_t linewidth, Double_t sigma_x, Double_t sigma_y,
                          Double_t offset_x, Double_t offset_y, Double_t offset_z,
                          Double_t yaw, Double_t pitch, Double_t roll);

    // Load a measured transverse-intensity profile (TH2D "h_profile", unit integral,
    // laser-frame axes in mm) for the most-recently-added 122/355 nm laser. Once set,
    // that laser uses a bilinear lookup into the histogram instead of the analytic
    // Gaussian, and its sigma_x/sigma_y are ignored.
    void SetLaser122Profile(const std::string& path);
    void SetLaser355Profile(const std::string& path);

    void SetLaserJitter(bool flag){jitter_on=flag;};

    // Resample every laser's parameters from its nominal value + sigma via a Gaussian draw.
    // No-op unless jitter_on is set. Must be called once per event, before SetMuPosition/PrecomputeAtPosition.
    void ResampleLaserPars();
//    void SetLinewidth(Double_t l){linewidth=l;};    // in GHz
//    void SetSigmaX(Double_t x){sigma_x=x;}; // in mm
//    void SetSigmaY(Double_t y){sigma_y=y;}; // in mm
//    void SetPulseTimeWidth(Double_t p){tau=p;}; // in ns
//    void SetPulseFWHM(Double_t p){tau=0.4247*p;};   // in ns
//    void SetLaserOffset(TVector3 r){laser_offset=r;};  // in mm
//    void SetLaserOffset355(TVector3 r){laser_offset_355=r;};  // in mm
//    void SetYawAngle(Double_t y){yaw = y * TMath::Pi()/180;};     // in deg
//    void SetPitchAngle(Double_t p){pitch = p * TMath::Pi()/180;}; // in deg
//    void SetRollAngle(Double_t r){roll = r * TMath::Pi()/180;};   // in deg
//    void SetYawAngle355(Double_t y){yaw_355 = y * TMath::Pi()/180;};     // in deg
//    void SetPitchAngle355(Double_t p){pitch_355 = p * TMath::Pi()/180;}; // in deg
//    void SetRollAngle355(Double_t r){roll_355 = r * TMath::Pi()/180;};   // in deg
//    void SetCentralFreq(Double_t freq){cen_freq=freq; laser_k=2*TMath::Pi()*cen_freq*1e9/299792458;}; // in GHz. k in m^-1
//    void SetWaveLength(Double_t wvl){laser_k=2*TMath::Pi()/wvl*1e9; cen_freq=299792458/wvl;}  // in nm
//    void SetPeakTime(Double_t t){peak_time=t;};

    void UpdateRotMat(Laser &lsr);

    // Precompute position-dependent (but time-independent) factors for the current event.
    // Must be called once per event whenever the muonium position changes.
    void PrecomputeAtPosition(TVector3 r);

    void SetOBESolverPtr(OBEsolver *ptr){obe_ptr=ptr;};

//    Double_t GetEnergy(){return energy;};
//    Double_t GetEnergy355(){return energy_355;};
//    Double_t GetLinewidth(){return linewidth;};
    Double_t GetLinewidth(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().linewidth;
    }
    Double_t GetSigmaX(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().sigma_x;
    }
    Double_t GetSigmaY(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().sigma_y;
    }
    Double_t GetPulseTimeWidth(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().tau;
    }
    Double_t GetPeakIntensity(TVector3 r);      // in W/cm^2
    Double_t GetPeakIntensity355(TVector3 r);   // in W/cm^2
    TVector3 GetWaveVector(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().laser_k * vec_laser122.front().laser_dirc;
    }   // in m^-1
    Double_t GetDetuning(){
        if (vec_laser122.empty()) throw std::runtime_error("LaserGenerator: no 122 nm laser configured");
        return vec_laser122.front().detuning;
    }  // in GHz

    TVector3 GetFieldE(TVector3 r, Double_t t); // in V/mm
    Double_t GetIntensity(TVector3 r, Double_t t);    // in W/cm^2
    Double_t GetIntensity355(TVector3 r, Double_t t);    // in W/cm^2

    // Live per-event laser parameter vectors (post-resampling if laser jitter is on).
    const std::vector<Laser>& GetLaser122Vec() const {return vec_laser122;};
    const std::vector<Laser>& GetLaser355Vec() const {return vec_laser355;};

private:
    // Private constructor and destructor
    LaserGenerator();
    ~LaserGenerator(){};

    TVector3 BeamToLaserCoord(TVector3 r, const Laser& lsr); // Transform from target coordinate to laser coordinate

    // Load "h_profile" from a ROOT file and attach it (shared) to the back() entry of
    // both the live and nominal laser vectors. tag is "122nm"/"355nm" for messages.
    void LoadProfileInto(std::vector<Laser>& live, std::vector<Laser>& nominal,
                         const std::string& path, const char* tag);
    // Bilinear lookup of the transverse density [mm^-2] at laser-frame (x,y).
    // Returns 0 outside the histogram's bin-centre range (no light, no ROOT error).
    Double_t SampleProfileDensity(const TH2D* h, Double_t x, Double_t y) const;

    // Per-event cached spatial factors, populated by PrecomputeAtPosition().
    // Each entry is (prefactor * exp_space) for the corresponding laser in the vector.
    // E-field and intensity use different prefactor forms, so separate caches are kept.
    std::vector<Double_t> cached_Espatial_122;   // for GetFieldE (E-field amplitude)
    std::vector<Double_t> cached_Ispatial_122;   // for GetIntensity (122 nm)
    std::vector<Double_t> cached_Ispatial_355;   // for GetIntensity355 (355 nm)
//    TVector3 BeamToLaserCoord355(TVector3 r); // Transform from target coordinate to laser coordinate

    OBEsolver *obe_ptr;

    std::vector<Laser> vec_laser122, vec_laser355;

    // Nominal (mean) values snapshotted at AddLaser122/355 time, and per-laser Gaussian sigmas
    // (default zero, set via SetLaser122Sigma/SetLaser355Sigma). ResampleLaserPars() draws
    // vec_laser122/355 from these each event when jitter_on is set.
    std::vector<Laser> vec_laser122_nominal, vec_laser355_nominal;
    std::vector<LaserSigma> vec_laser122_sigma, vec_laser355_sigma;
    bool jitter_on = false;

    // Redraw one Laser's fields from nominal/sigma via RunManager's shared RNG, enforcing
    // physical validity (sigma_x/sigma_y/tau > 0, energy/linewidth >= 0) with bounded retries.
    void ResampleOneLaser(Laser &live, const Laser &nominal, const LaserSigma &sigma, bool has_detuning);

//    Double_t energy;
//    Double_t energy_355;
//    Double_t linewidth;
//    Double_t sigma_x;
//    Double_t sigma_y;
//    Double_t tau;
//    Double_t peak_time;
//    Double_t yaw;
//    Double_t pitch;
//    Double_t roll;
//    Double_t yaw_355;
//    Double_t pitch_355;
//    Double_t roll_355;
//
//    TVector3 laser_offset;
//    TVector3 laser_offset_355;
//
//    Double_t cen_freq;
//    Double_t laser_k;
//    TVector3 laser_dirc;
//
//    Eigen::Matrix3d rot_mat_122;
//    Eigen::Matrix3d rot_mat_355;
//    Eigen::Matrix3d rot_mat_rev_122;

};

#endif //LASERTOYMC_LASERGENERATOR_H
