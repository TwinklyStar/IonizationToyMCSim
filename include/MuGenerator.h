//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_MUGENERATOR_H
#define LASERTOYMC_MUGENERATOR_H
#include "common.h"
#include "LaserGenerator.h"

class MuGenerator {
public:
    // Meyers' Singleton - Get the instance of the class
    static MuGenerator& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    MuGenerator(const MuGenerator&) = delete;
    MuGenerator& operator=(const MuGenerator&) = delete;

    void SetTemperature(Double_t temp){temperature=temp;}; // in Kelvin

    Double_t GetTemperature(){return temperature;}; // in Kelvin

    void SetRndSeed(UInt_t rndseed){randGen.SetSeed(rndseed);};

    TVector3 SampleVelocity();  // in m/s

    TVector3 SampleLocation();  // in mm


private:
    // Private constructor and destructor
    MuGenerator(): randGen(0){temperature=322; };
    ~MuGenerator(){};

    Double_t temperature;
    TRandom randGen;

    const Double_t k_B = 1.380649e-23; // Boltzmann constant in J/K
    const Double_t m_Mu = 1.8926410106e-28;  // Muonium mass in kg

};

#endif //LASERTOYMC_MUGENERATOR_H
