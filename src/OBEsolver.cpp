//
// Created by Meng Lv on 2024/10/9.
//

#include "OBEsolver.h"

OBEsolver::OBEsolver(Double_t gm1, Double_t gmion):gamma_1(gm1), gamma_ion(gmion)
{
    laser_ptr = &LaserGenerator::GetInstance();
//    laser_ptr->SetOBESolverPtr(this);
    Mu_ptr = &MuGenerator::GetInstance();
    ROOT_ptr = &RootManager::GetInstance();

    UpdateGamma2();

    initial_rho = state_type (4);
    SetInitialState(1, 0, 0, 0);
    SetMuPosition({0,0,0});
    SetMuVelocity({0, 0, 0});
    SetStartTime(0);
    SetEndTime(10);
    SetDt(0.01);
    SetAbsErr(1e-8);
    SetRelErr(1e-6);
}

void OBEsolver::UpdateGamma2() {
    gamma_2 = gamma_1/2 + gamma_ion/2 + laser_ptr->GetLinewidth()*TMath::Pi();
}

void OBEsolver::OBE(const state_type &rho, state_type &drhodt, const Double_t t){
    // Extracting the elements of the density matrix
    complex<Double_t> rho_gg = rho[0];  // ρ_gg
    complex<Double_t> rho_ee = rho[1];  // ρ_ee
    complex<Double_t> rho_ge = rho[2];  // ρ_ge (coherence term)
    complex<Double_t> rho_eg = conj(rho_ge); // ρ_eg (conjugate of ρ_ge)
    complex<Double_t> rho_ion = rho[3]; // ρ_ion (ionization)

    // Time-dependent detuning and Rabi frequency
    complex<Double_t> Omega_t = GetRabiFreq(t);
    Double_t delta_t = GetDopplerShift();

    // Equations based on the system provided
    drhodt[0] = gamma_1 * rho_ee + (complex<Double_t>(0.0, 0.5)) * (conj(Omega_t) * rho_eg - Omega_t * rho_ge);
    drhodt[1] = -(gamma_1 + gamma_ion) * rho_ee - (complex<Double_t>(0.0, 0.5)) * (conj(Omega_t) * rho_eg - Omega_t * rho_ge);
    drhodt[2] = -(gamma_2 + complex<Double_t>(0.0, 1.0) * delta_t + gamma_ion / 2.0) * rho_ge
                + (complex<Double_t>(0.0, 0.5)) * conj(Omega_t) * (rho_ee - rho_gg);
    drhodt[3] = gamma_ion * rho_ee; // dρ_ion/dt
}

complex<Double_t> OBEsolver::GetRabiFreq(const Double_t t) {
    Double_t e_field_y = laser_ptr->GetFieldE(Mu_pos, t).Y();  // in V/mm, assume only linearly polarized in y direction
    return 6.0189 * e_field_y * pow(10, -2);    // Return frequency in GHz

}

Double_t OBEsolver::GetDopplerShift() {
    return 0;
//    return -Mu_v.Dot(laser_ptr->GetWaveVector())*1e-9;

}

void OBEsolver::solve() {
    typedef runge_kutta_cash_karp54<state_type> stepper_type;
    auto stepper = make_controlled(abs_err, rel_err, stepper_type());

    // Use std::bind to bind the system and observer to the class instance so that the function can be used as input parameters
    auto sys = std::bind(&OBEsolver::OBE, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    auto obs = std::bind(&OBEsolver::Observer, this, std::placeholders::_1, std::placeholders::_2);

    auto rho_temp = initial_rho;    // initial_rho will not change during integration
    // Perform integration
    integrate_adaptive(stepper, sys, rho_temp, start_time, end_time, dt, obs);

}

void OBEsolver::Observer(const state_type &rho, Double_t t) {
//    std::cout << rho[0].real() << " ";
    ROOT_ptr->PushTimePoint(t, laser_ptr->GetFieldE(Mu_pos, t).Mag(), GetRabiFreq(t).real(),
                            rho[0].real(), rho[1].real(),
                            rho[2].real(), rho[2].imag(), rho[3].real());
}
//void OBEsolver::Observer(const state_type &rho, Double_t t) {
//    std::cout << t << ", " << GetRabiFreq(t).real() << ","
//}
