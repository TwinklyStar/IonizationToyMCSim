//
// Created by Meng Lv on 2024/11/20.
//
#include "RunManager.h"

RunManager::RunManager() {
    ROOT_ptr = &RootManager::GetInstance();
    Mu_ptr = &MuGenerator::GetInstance();
    lsr_ptr = &LaserGenerator::GetInstance();
    solver = new OBEsolver(0.627);
}

RunManager& RunManager::GetInstance() {
    static RunManager instance;
    return instance;
}

void RunManager::ReadCommandFile(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + file_path);
    }

    std::string line;
    while (std::getline(file, line)) {
        // Trim leading and trailing whitespaces
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // Ignore comment lines starting with #
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string command, sub_command, last_word;
        iss >> command;

        if (command == "AddLaser122") {
            iss >> pulse_energy_122 >> pulse_FWHM_122 >> peak_time_122
                >> linewidth_122 >> sigma_x_122 >> sigma_y_122
                >> offset_x_122 >> offset_y_122 >> offset_z_122
                >> yaw_122 >> pitch_122 >> roll_122;

            lsr_ptr->AddLaser122(pulse_energy_122, pulse_FWHM_122, peak_time_122,
                                 linewidth_122, sigma_x_122, sigma_y_122,
                                 offset_x_122, offset_y_122, offset_z_122,
                                 yaw_122, pitch_122, roll_122);

        } else if (command == "RandomSeed") {
            iss >> rdm_seed;
            rdm_gen.SetSeed(rdm_seed);
            std::cout << "-- RunManager: Set random seed as " << rdm_seed << std::endl;
        } else if (command == "AddLaser355") {
            iss >> pulse_energy_355 >> pulse_FWHM_355 >> peak_time_355
                >> linewidth_355 >> sigma_x_355 >> sigma_y_355
                >> offset_x_355 >> offset_y_355 >> offset_z_355
                >> yaw_355 >> pitch_355 >> roll_355;

            lsr_ptr->AddLaser355(pulse_energy_355, pulse_FWHM_355, peak_time_355,
                                 linewidth_355, sigma_x_355, sigma_y_355,
                                 offset_x_355, offset_y_355, offset_z_355,
                                 yaw_355, pitch_355, roll_355);

        } else if (command == "MuInputFile") {
            iss >> input_file_name;
            std::cout << "-- RunManager: Input file: " << input_file_name << std::endl;
        } else if (command == "OutputFile") {
            iss >> output_file_name;
            std::cout << "-- RunManager: Output file: " << output_file_name << std::endl;
            ROOT_ptr->SetOutFileName(output_file_name);
        } else if (command == "SetRunTime"){
            iss >> runtime;
            solver->SetEndTime(runtime);
            std::cout << "-- RunManager: Set simulation time: " << runtime << " ns" << std::endl;
        } else if (command == "RootOutput") {
            iss >> sub_command >> last_word;
            if (sub_command == "t") {
                ROOT_ptr->Sett((last_word == "on"));
            } else if (sub_command == "RabiFreq") {
                ROOT_ptr->SetRabiFreq((last_word == "on"));
            } else if (sub_command == "EField") {
                ROOT_ptr->SetEField((last_word == "on"));
            } else if (sub_command == "Intensity122") {
                ROOT_ptr->SetIntensity122((last_word == "on"));
            } else if (sub_command == "Intensity355") {
                ROOT_ptr->SetIntensity355((last_word == "on"));
            } else if (sub_command == "GammaIon") {
                ROOT_ptr->SetGammaIon((last_word == "on"));
            } else if (sub_command == "rho_gg") {
                ROOT_ptr->Setrho_gg((last_word == "on"));
            } else if (sub_command == "rho_ee") {
                ROOT_ptr->Setrho_ee((last_word == "on"));
            } else if (sub_command == "rho_ge_r") {
                ROOT_ptr->Setrho_ge_r((last_word == "on"));
            } else if (sub_command == "rho_ge_i") {
                ROOT_ptr->Setrho_ge_i((last_word == "on"));
            } else if (sub_command == "rho_ion") {
                ROOT_ptr->Setrho_ion((last_word == "on"));
            } else {
                std::cerr << "WARNING RunManager: Unknown branch: " << sub_command << std::endl;
            }
        } else {
            std::cerr << "WARNING RunManager: Unknown cammand: " << command << std::endl;
        }
    }
    file.close();

    InitializeMuGenerator();
    InitializeRootManager();
//    InitializeOBEsolver();
}

void RunManager::InitializeMuGenerator() {
    Mu_ptr->ReadInputFile(input_file_name);
}

void RunManager::InitializeRootManager() {
    ROOT_ptr->Initialize();
}

void RunManager::InitializeOBEsolver() {
    solver->SetStartTime(0);        // in ns
    solver->SetEndTime(10);         // in ns
    solver->SetDt(0.001);            // in ns
    solver->SetAbsErr(1e-8);
    solver->SetRelErr(1e-6);
    solver->SetInitialState(1, 0, 0, 0);
}

// Simulation functions
void RunManager::SolveOBE() {
    int eventn = Mu_ptr->GetInputEventNum();
    std::cout << "--- Number of events: " << eventn << std::endl;
    for(int i=0; i<eventn; i++){
        if (eventn >= 100 && i % (eventn/100) == 0) loader(i/(eventn/100));

        solver->SetMuPosition(Mu_ptr->GetInputLocation(i));
        solver->SetMuVelocity(Mu_ptr->GetInputVelocity(i));

        solver->solve();

        ROOT_ptr->SetEventID(i);
//        ROOT_ptr->SetLaserPars(lsr_ptr->GetEnergy(), lsr_ptr->GetEnergy355(), lsr_ptr->GetPulseTimeWidth(),
//                               lsr_ptr->GetSigmaX(), lsr_ptr->GetSigmaY(),
//                               lsr_ptr->GetPeakIntensity(solver->GetMuPosition()),
//                               lsr_ptr->GetPeakIntensity355(solver->GetMuPosition()), lsr_ptr->GetLinewidth());
        ROOT_ptr->SetDoppFreq(solver->GetDopplerShift());
        ROOT_ptr->SetPosition(solver->GetMuPosition());
        ROOT_ptr->SetVelocity(solver->GetMuVelocity());
        ROOT_ptr->SetLastState();

        ROOT_ptr->FillEvent();
    }

    ROOT_ptr->Finalize();
}

void RunManager::SolveOBETest() {
    int eventn = Mu_ptr->GetInputEventNum();
    std::cout << "--- Number of events: " << eventn << std::endl;
    for(int i=0; i<eventn; i++){
        if (eventn >= 100 && i % (eventn/100) == 0) loader(i/(eventn/100));

//        Double_t linewidth_arr[100];
//        Double_t dopp_arr[100];
//        for (int j=0; j<100; j++){
//            linewidth_arr[j] = j;
//            dopp_arr[j] = -100 + j*200/100.;
//        }
//        lsr_ptr->SetEnergy(10e-6);
//        solver->SetMuPosition(Mu_ptr->SampleLocation());
//        solver->SetMuVelocity(Mu_ptr->SampleVelocity());
//        solver->SetMuPosition({0,0,0});
//        solver->SetMuVelocity({0,0,0});
        solver->SetMuPosition(Mu_ptr->GetInputLocation(i));
        solver->SetMuVelocity(Mu_ptr->GetInputVelocity(i));

//        lsr_ptr->SetLinewidth(linewidth_arr[i%100]);
//        lsr_ptr->SetLinewidth(80);
//        solver->SetDopplerShift(dopp_arr[i/100]);

        solver->solve();

        ROOT_ptr->SetEventID(i);
//        ROOT_ptr->SetLaserPars(lsr_ptr->GetEnergy(), lsr_ptr->GetEnergy355(), lsr_ptr->GetPulseTimeWidth(),
//                               lsr_ptr->GetSigmaX(), lsr_ptr->GetSigmaY(),
//                               lsr_ptr->GetPeakIntensity(solver->GetMuPosition()),
//                               lsr_ptr->GetPeakIntensity355(solver->GetMuPosition()), lsr_ptr->GetLinewidth());
        ROOT_ptr->SetDoppFreq(solver->GetDopplerShift());
        ROOT_ptr->SetPosition(solver->GetMuPosition());
        ROOT_ptr->SetVelocity(solver->GetMuVelocity());
        ROOT_ptr->SetLastState();

        ROOT_ptr->FillEvent();
    }

    ROOT_ptr->Finalize();

}

void RunManager::parTestBench() {
    int i;  // event loop index
    Double_t t_pos_x, t_pos_y, t_pos_z, t_v_x, t_v_y, t_v_z, t_dopp, t_Int, t_Int_355;
    TFile *ff = TFile::Open("data/test17.root", "RECREATE");
    TTree *t1 = new TTree("pars", "A tree of parameters");
//    t1->Branch("eventn", &i);
    t1->Branch("x", &t_pos_x);
    t1->Branch("y", &t_pos_y);
    t1->Branch("z", &t_pos_z);
//    t1->Branch("vx", &t_v_x);
//    t1->Branch("vy", &t_v_y);
//    t1->Branch("vz", &t_v_z);
//    t1->Branch("doppler", &t_dopp);
    t1->Branch("Intensity", &t_Int); // at t=0
    t1->Branch("Intensity355", &t_Int_355); // at t=0

    for (int x_idx=0; x_idx<100; x_idx++){
        for (int y_idx=0; y_idx<100; y_idx++){
            for (int z_idx=0; z_idx<100; z_idx++){
                t_pos_x = -10 + 0.2*x_idx;
                t_pos_y = -10 + 0.2*y_idx;
                t_pos_z = -10 + 0.2*z_idx;

                t_Int = lsr_ptr->GetIntensity({t_pos_x, t_pos_y, t_pos_z}, 0);
                t_Int_355 = lsr_ptr->GetIntensity355({t_pos_x, t_pos_y, t_pos_z}, 0);

                t1->Fill();
            }
        }
    }


    ff->Write();
}

void RunManager::loader(int rate) {
    char proc[22];
    memset(proc, '\0', sizeof(proc));

    for (int i = 0; i < rate/5; i++)
    {
        proc[i] = '#';
    }

    printf("\r[%-20s] [%d%%]", proc, rate);        //C语言格式控制时默认右对齐，所以要在前面加-变成左对齐
    fflush(stdout);                                 //刷新屏幕打印
}
