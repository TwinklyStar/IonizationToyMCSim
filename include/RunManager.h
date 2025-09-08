//
// Created by Meng Lv on 2024/11/20.
//

#ifndef IONIZATIONTOYMCSIM_RUNMANAGER_H
#define IONIZATIONTOYMCSIM_RUNMANAGER_H
#include "common.h"
#include "LaserGenerator.h"
#include "RootManager.h"
#include "MuGenerator.h"
#include "OBEsolver.h"

class RunManager {
public:
    // Meyer's Singleton Instance
    static RunManager& GetInstance();

    // Function to read commands from a file
    void ReadCommandFile(const std::string& file_path);

    // Member functions to get the parsed values
    void InitializeLaserGenerator();
    void InitializeMuGenerator();
    void InitializeRootManager();
    void InitializeOBEsolver();

    // Main simulation functions
    void SolveOBE();
    void SolveOBETest();
    void parTestBench();

    // Random generator
    TRandom rdm_gen;
    void SetRdmSeed(int s) {rdm_gen.SetSeed(s);};

private:
    // Private constructor for Singleton
    RunManager();
    ~RunManager() {};

    RootManager *ROOT_ptr;
    LaserGenerator *lsr_ptr;
    MuGenerator *Mu_ptr;
    OBEsolver *solver;

    // Variables for AddLaser122 command
    double pulse_energy_122 = 0.0, pulse_FWHM_122 = 0.0, peak_time_122 = 0.0;
    double linewidth_122 = 0.0, sigma_x_122 = 0.0, sigma_y_122 = 0.0;
    double offset_x_122 = 0.0, offset_y_122 = 0.0, offset_z_122 = 0.0;
    double yaw_122 = 0.0, pitch_122 = 0.0, roll_122 = 0.0;
    int rdm_seed = 0;
    double runtime;

    // Variables for AddLaser355 command
    double pulse_energy_355 = 0.0, pulse_FWHM_355 = 0.0, peak_time_355 = 0.0;
    double linewidth_355 = 0.0, sigma_x_355 = 0.0, sigma_y_355 = 0.0;
    double offset_x_355 = 0.0, offset_y_355 = 0.0, offset_z_355 = 0.0;
    double yaw_355 = 0.0, pitch_355 = 0.0, roll_355 = 0.0;

    // Variables for MuInputFile and OutputFile commands
    std::string input_Mu="on", input_file_name, output_file_name, event_n;

    // Progress bar
    void loader(int rate);
};

#endif //IONIZATIONTOYMCSIM_RUNMANAGER_H
