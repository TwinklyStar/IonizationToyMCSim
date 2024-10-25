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
//    return {0, 0, 0};
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

void MuGenerator::ReadInputFile(std::string infile) {
    // Open the file
    std::ifstream file(infile);
    if (!file.is_open()) {
        std::cerr << "Error opening file " << infile << std::endl;
        return;
    }

    std::string line;
    // Read each line of the file
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double value1, value2, value3, value4, value5, value6;

        // Extract the six values from the line, separated by tabs
        if (ss >> value1 >> value2 >> value3 >> value4 >> value5 >> value6) {
            input_x.push_back(value1);
            input_y.push_back(value2);
            input_z.push_back(value3);
            input_vx.push_back(value4/1e3);
            input_vy.push_back(value5/1e3);
            input_vz.push_back(value6/1e3);
        } else {
            std::cerr << "Error reading line: " << line << std::endl;
            return;
        }
    }

    file.close();
    input_n = input_x.size();
}

TVector3 MuGenerator::GetInputLocation(int i) {
    if (input_n == 0){
        std::cerr << "MuGenerator::GetInputLocation: Input file not registered yet. Return {0,0,0}" << std::endl;
        return {0,0,0};
    }
    if (i >= input_n){
        std::cerr << "MuGenerator::GetInputLocation: index " << i << " exceeds total number " << input_n << ". Return {0,0,0}" << std::endl;
        return {0,0,0};
    }
    return {input_x.at(i), input_y.at(i), input_z.at(i)};
}

TVector3 MuGenerator::GetInputVelocity(int i) {
    if (input_n == 0){
        std::cerr << "MuGenerator::GetInputVelocity: Input file not registered yet. Return {0,0,0}" << std::endl;
        return {0,0,0};
    }
    if (i >= input_n){
        std::cerr << "MuGenerator::GetInputVelocity: index " << i << " exceeds total number " << input_n << ". Return {0,0,0}" << std::endl;
        return {0,0,0};
    }
    return {input_vx.at(i), input_vy.at(i), input_vz.at(i)};
}