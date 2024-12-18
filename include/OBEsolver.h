//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_OBESOLVER_H
#define LASERTOYMC_OBESOLVER_H

#include "common.h"
#include "LaserGenerator.h"
#include "MuGenerator.h"
#include "RootManager.h"

using namespace std;
using namespace boost::numeric::odeint;

typedef vector<complex<Double_t>> state_type;

class OBEsolver {
public:
    OBEsolver(Double_t gm1);
    ~OBEsolver(){};

    void solve();

    void SetGamma1(Double_t gm1){gamma_1=gm1;}
    void SetStartTime(Double_t t){start_time=t;};
    void SetEndTime(Double_t t){end_time=t;};
    void SetDt(Double_t t){dt=t;};
    void SetMuPosition(TVector3 pos){Mu_pos = pos;}
    void SetMuVelocity(TVector3 v){Mu_v = v;}
    void SetInitialState(Double_t rho_gg, Double_t rho_ee, Double_t rho_eg, Double_t rho_ion){
        initial_rho[0] = rho_gg;
        initial_rho[1] = rho_ee;
        initial_rho[2] = rho_eg;
        initial_rho[3] = rho_ion;
    }
    void SetAbsErr(Double_t err){abs_err=err;};
    void SetRelErr(Double_t err){rel_err=err;};
    void SetDopplerShift(Double_t shift){detuning=shift; if_set_dopp= true;};


    complex<Double_t> GetRabiFreq(const Double_t t);  // in GHz
    Double_t GetDopplerShift();

    Double_t GetGamma1(){return gamma_1;}
//    Double_t GetGamma2(){return gamma_2;}
    Double_t GetGammaIon(Double_t t);
    TVector3 GetMuPosition(){return Mu_pos;}
    TVector3 GetMuVelocity(){return Mu_v;}


private:
    Double_t gamma_1;
//    Double_t gamma_2;
    Double_t start_time;
    Double_t end_time;
    Double_t dt;
    Double_t abs_err;
    Double_t rel_err;
    LaserGenerator *laser_ptr;
    MuGenerator *Mu_ptr;
    RootManager *ROOT_ptr;
    TVector3 Mu_pos;
    TVector3 Mu_v;

    bool if_set_dopp=false;
    Double_t detuning=0;

    state_type initial_rho; // rho_ee, gg, eg, ion

    void OBE(const state_type &rho, state_type &drhodt, const Double_t t);
    void Observer(const state_type &rho, Double_t t);

};


#endif //LASERTOYMC_OBESOLVER_H
