//
// Created by Meng Lv on 2024/10/9.
//

#ifndef LASERTOYMC_MUGENERATOR_H
#define LASERTOYMC_MUGENERATOR_H
#include "common.h"
#include "LaserGenerator.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

class RunManager;

class MuGenerator {
public:
    // Meyers' Singleton - Get the instance of the class
    static MuGenerator& GetInstance();

    // Delete copy constructor and assignment operator to avoid copying
    MuGenerator(const MuGenerator&) = delete;
    MuGenerator& operator=(const MuGenerator&) = delete;

    void SetTemperature(Double_t temp){temperature=temp;}; // in Kelvin

    void SetIfInputMu(const std::string& in){input_Mu = !(in=="off");}

    void SetEventN(std::string n) {event_n = n;}

    void ReadInputFile(std::string infile);

    int GetInputEventNum(){return input_n;};

    Double_t GetTemperature(){return temperature;}; // in Kelvin

    TVector3 SampleVelocity();  // in m/s

    TVector3 SampleLocation();  // in mm

    TVector3 GetInputLocation(int i);   // in mm
    TVector3 GetInputVelocity(int i);   // in m/s


private:
    // Private constructor and destructor
    MuGenerator();
    ~MuGenerator(){};

    Double_t temperature;
//    TRandom randGen;

    const Double_t k_B = 1.380649e-23; // Boltzmann constant in J/K
    const Double_t m_Mu = 1.8926410106e-28;  // Muonium mass in kg

    bool input_Mu = true;
    std::string event_n;

    int input_n=0;
    std::vector<Double_t> input_x;
    std::vector<Double_t> input_y;
    std::vector<Double_t> input_z;
    std::vector<Double_t> input_vx;
    std::vector<Double_t> input_vy;
    std::vector<Double_t> input_vz;

};

#endif //LASERTOYMC_MUGENERATOR_H
