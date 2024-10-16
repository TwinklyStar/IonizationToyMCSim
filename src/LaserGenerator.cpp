//
// Created by Meng Lv on 2024/10/9.
//

#include "LaserGenerator.h"

// Meyers' Singleton implementation
LaserGenerator& LaserGenerator::GetInstance() {
    static LaserGenerator instance;
    return instance;
}

TVector3 LaserGenerator::GetFieldE(TVector3 r, Double_t t) {
    Double_t x = r.Y(); // simple version
    Double_t y = r.Z() - laser_z; // simple version
    Double_t z = r.X(); // simple version
    const Double_t eta = 376.7303134;

    Double_t prefactor = sqrt(2.0 * sqrt(2.0) * eta * energy / (pow(TMath::Pi(), 3.0 / 2.0) * sigma_x * sigma_y  * tau *
            pow(10, -9)));      // convert ns to s

    // Compute the spatial Gaussian components
    double exp_space = exp(-(x * x) / (sigma_x * sigma_x) - (y * y) / (sigma_y * sigma_y));

    // Compute the temporal Gaussian component
    double exp_time = exp(-((t-peak_time) * (t-peak_time)) / (4.0 * tau * tau));

    // Combine everything to compute the electric field
    double E = prefactor * exp_space * exp_time;

    return {0, E, 0};   // x-polarized in muonium coordinate
}
