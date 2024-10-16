//
// Created by Meng Lv on 2024/10/9.
//

#include "MuGenerator.h"

// Meyers' Singleton implementation
MuGenerator& MuGenerator::GetInstance() {
    static MuGenerator instance;
    return instance;
}

TVector3 MuGenerator::SampleVelocity() {
    Double_t v_thermal = sqrt(k_B * temperature / m_Mu);

    return {randGen.Gaus(0, v_thermal),
            randGen.Gaus(0, v_thermal),
            randGen.Gaus(0, v_thermal)};
}

TVector3 MuGenerator::SampleLocation() {
    // At first, we assume uniform distribution in +- 3 sigma of laser region
    // And laser field in z direction is uniform. We simply set z as 0
    LaserGenerator &lsr = LaserGenerator::GetInstance();

    return {0,
            randGen.Uniform(-3*lsr.GetSigmaX(), 3*lsr.GetSigmaX()),
            randGen.Uniform(-3*lsr.GetSigmaY(), 3*lsr.GetSigmaY())};
//    return {0,0,0};
}