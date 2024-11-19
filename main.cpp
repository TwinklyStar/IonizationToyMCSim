#include "common.h"
#include "OBEsolver.h"
#include "LaserGenerator.h"
#include "MuGenerator.h"
#include "RootManager.h"

// Meyer's Singleton Pattern
LaserGenerator *lsr_ptr = &LaserGenerator::GetInstance();
MuGenerator *Mu_ptr = &MuGenerator::GetInstance();
RootManager *ROOT_ptr = &RootManager::GetInstance("data/OBEtest31.root");
OBEsolver *solver = new OBEsolver(0.627);

void parTestBench(int eventn);
void SolveOBE(int eventn);
void loader(int rate);

int main() {

    lsr_ptr->SetLaserOffset({0, 0, 3.3});

    lsr_ptr->SetYawAngle(0);
    lsr_ptr->SetPitchAngle(0);
    lsr_ptr->SetRollAngle(0);
    lsr_ptr->SetLaserOffset355({0, 0, 0});
    lsr_ptr->SetYawAngle355(0);
    lsr_ptr->SetPitchAngle355(0);
    lsr_ptr->SetRollAngle355(0);

    lsr_ptr->SetSigmaX(4);          // in mm
    lsr_ptr->SetSigmaY(1);          // in mm
//    lsr_ptr->SetPulseTimeWidth(1);  // in ns
    lsr_ptr->SetPulseFWHM(2);       // in ns
    lsr_ptr->SetEnergy(13.5e-6);      // in J
    lsr_ptr->SetEnergy355(8e-3);     // in J
//    lsr_ptr->SetEnergy(0);      // in J
    lsr_ptr->SetLinewidth(80);       // in GHz
    lsr_ptr->SetWaveLength(122);  // in nm
    lsr_ptr->SetLaserDirection({1, 0, 0});  // vector direction
    lsr_ptr->SetPeakTime(5);        // in ns

    Mu_ptr->SetRndSeed(999);
    Mu_ptr->SetTemperature(322);
    Mu_ptr->ReadInputFile("../datasets/test1k.dat");

    ROOT_ptr->SetRndSeed(1000); // A different seed from Mu_ptr

//    parTestBench(1);


    solver->SetStartTime(0);        // in ns
    solver->SetEndTime(10);         // in ns
    solver->SetDt(0.001);            // in ns
    solver->SetAbsErr(1e-8);
    solver->SetRelErr(1e-6);
    solver->SetInitialState(1, 0, 0, 0);

//    int eventn = 10000;
    int eventn = Mu_ptr->GetInputEventNum();
    std::cout << "--- Number of events: " << eventn << std::endl;
//    parTestBench(eventn);
    SolveOBE(eventn);
//    SolveOBE(10);

    std::cout << "\nHello, World!" << std::endl;
    return 0;
}

void SolveOBE(int eventn){
    for(int i=0; i<eventn; i++){
        if (eventn >= 100 && i % (eventn/100) == 0) loader(i/(eventn/100));

        solver->SetMuPosition(Mu_ptr->GetInputLocation(i));
        solver->SetMuVelocity(Mu_ptr->GetInputVelocity(i));

        solver->solve();

        ROOT_ptr->SetEventID(i);
        ROOT_ptr->SetLaserPars(lsr_ptr->GetEnergy(), lsr_ptr->GetEnergy355(), lsr_ptr->GetPulseTimeWidth(),
                               lsr_ptr->GetSigmaX(), lsr_ptr->GetSigmaY(),
                               lsr_ptr->GetPeakIntensity(solver->GetMuPosition()),
                               lsr_ptr->GetPeakIntensity355(solver->GetMuPosition()), lsr_ptr->GetLinewidth());
        ROOT_ptr->SetDoppFreq(solver->GetDopplerShift());
        ROOT_ptr->SetPosition(solver->GetMuPosition());
        ROOT_ptr->SetVelocity(solver->GetMuVelocity());
        ROOT_ptr->SetLastState();

        ROOT_ptr->FillEvent();
    }

    ROOT_ptr->Finalize();

}

void SolveOBETest(int eventn){
    for(int i=0; i<eventn; i++){
        if (eventn >= 100 && i % (eventn/100) == 0) loader(i/(eventn/100));

//        Double_t linewidth_arr[100];
//        Double_t dopp_arr[100];
//        for (int j=0; j<100; j++){
//            linewidth_arr[j] = j;
//            dopp_arr[j] = -100 + j*200/100.;
//        }
        lsr_ptr->SetEnergy(10e-6);
//        solver->SetMuPosition(Mu_ptr->SampleLocation());
//        solver->SetMuVelocity(Mu_ptr->SampleVelocity());
//        solver->SetMuPosition({0,0,0});
//        solver->SetMuVelocity({0,0,0});
        solver->SetMuPosition(Mu_ptr->GetInputLocation(i));
        solver->SetMuVelocity(Mu_ptr->GetInputVelocity(i));

//        lsr_ptr->SetLinewidth(linewidth_arr[i%100]);
        lsr_ptr->SetLinewidth(80);
//        solver->SetDopplerShift(dopp_arr[i/100]);

        solver->solve();

        ROOT_ptr->SetEventID(i);
        ROOT_ptr->SetLaserPars(lsr_ptr->GetEnergy(), lsr_ptr->GetEnergy355(), lsr_ptr->GetPulseTimeWidth(),
                               lsr_ptr->GetSigmaX(), lsr_ptr->GetSigmaY(),
                               lsr_ptr->GetPeakIntensity(solver->GetMuPosition()),
                               lsr_ptr->GetPeakIntensity355(solver->GetMuPosition()), lsr_ptr->GetLinewidth());
        ROOT_ptr->SetDoppFreq(solver->GetDopplerShift());
        ROOT_ptr->SetPosition(solver->GetMuPosition());
        ROOT_ptr->SetVelocity(solver->GetMuVelocity());
        ROOT_ptr->SetLastState();

        ROOT_ptr->FillEvent();
    }

    ROOT_ptr->Finalize();

}

void parTestBench(int eventn) {
    int i;  // event loop index
    Double_t t_pos_x, t_pos_y, t_pos_z, t_v_x, t_v_y, t_v_z, t_dopp, t_Int, t_Int_355;
    TFile *ff = TFile::Open("data/test14.root", "RECREATE");
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

// Progress bar
void loader(int rate)
{
    char proc[22];
    memset(proc, '\0', sizeof(proc));

    for (int i = 0; i < rate/5; i++)
    {
        proc[i] = '#';
    }

    printf("\r[%-20s] [%d%%]", proc, rate);        //C语言格式控制时默认右对齐，所以要在前面加-变成左对齐
    fflush(stdout);                                 //刷新屏幕打印
}
