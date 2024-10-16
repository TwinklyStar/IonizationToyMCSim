#include "common.h"
#include "OBEsolver.h"
#include "LaserGenerator.h"
#include "MuGenerator.h"
#include "RootManager.h"

// Meyer's Singleton Pattern
LaserGenerator *lsr_ptr = &LaserGenerator::GetInstance();
MuGenerator *Mu_ptr = &MuGenerator::GetInstance();
RootManager *ROOT_ptr = &RootManager::GetInstance("OBEtest00.root");
OBEsolver *solver = new OBEsolver(0.627, 1.5);

void parTestBench(int eventn);
void SolveOBE(int eventn);

int main() {

    lsr_ptr->SetLaserPosition(0);
    lsr_ptr->SetSigmaX(4);          // in mm
    lsr_ptr->SetSigmaY(1);          // in mm
    lsr_ptr->SetPulseTimeWidth(1);  // in ns
    lsr_ptr->SetEnergy(10e-6);      // in J
    lsr_ptr->SetLinewidth(8);       // in GHz
    lsr_ptr->SetWaveLength(122);  // in nm
    lsr_ptr->SetLaserDirection({1, 0, 0});  // vector direction
    lsr_ptr->SetPeakTime(5);        // in ns

    Mu_ptr->SetRndSeed(999);
    Mu_ptr->SetTemperature(322);

    solver->SetStartTime(0);        // in ns
    solver->SetEndTime(10);         // in ns
    solver->SetDt(0.001);            // in ns

    int eventn = 10;
//    parTestBench(eventn);
    SolveOBE(eventn);



    std::cout << "Hello, World!" << std::endl;
    return 0;
}

void SolveOBE(int eventn){
    solver->SetInitialState(1, 0, 0, 0);
    for(int i=0; i<eventn; i++){
        solver->SetMuPosition(Mu_ptr->SampleLocation());
        solver->SetMuVelocity(Mu_ptr->SampleVelocity());

        solver->solve();

        ROOT_ptr->SetEventID(i);
        ROOT_ptr->SetDoppFreq(solver->GetDopplerShift());
        ROOT_ptr->SetPosition(solver->GetMuPosition());
        ROOT_ptr->SetVelocity(solver->GetMuVelocity());

        ROOT_ptr->FillEvent();
    }

    ROOT_ptr->Finalize();

}

void parTestBench(int eventn) {
    int i;  // event loop index
    Double_t t_pos_x, t_pos_y, t_pos_z, t_v_x, t_v_y, t_v_z, t_dopp, t_rabi, t_E;
    TFile *ff = TFile::Open("test00.root", "RECREATE");
    TTree *t1 = new TTree("pars", "A tree of parameters");
    t1->Branch("eventn", &i);
    t1->Branch("x", &t_pos_x);
    t1->Branch("y", &t_pos_y);
    t1->Branch("z", &t_pos_z);
    t1->Branch("vx", &t_v_x);
    t1->Branch("vy", &t_v_y);
    t1->Branch("vz", &t_v_z);
    t1->Branch("doppler", &t_dopp);
    t1->Branch("Efield", &t_E); // at t=0
    t1->Branch("rabi", &t_rabi);    // at t=0

    for (i = 0; i < eventn; i++) {
        solver->SetMuPosition(Mu_ptr->SampleLocation());
        solver->SetMuVelocity(Mu_ptr->SampleVelocity());
        t_pos_x = solver->GetMuPosition().X();
        t_pos_y = solver->GetMuPosition().Y();
        t_pos_z = solver->GetMuPosition().Z();
        t_v_x = solver->GetMuVelocity().X();
        t_v_y = solver->GetMuVelocity().Y();
        t_v_z = solver->GetMuVelocity().Z();
        t_dopp = solver->GetDopplerShift();
        t_E = lsr_ptr->GetFieldE(solver->GetMuPosition(), 0).Y();
        t_rabi = solver->GetRabiFreq(0).real();

        t1->Fill();
    }

    TTree *t2 = new TTree("pars_t", "A tree of time-dependent parameters");
    Double_t t; // in ns
    t2->Branch("t", &t);
    t2->Branch("Efield", &t_E);
    t2->Branch("rabi", &t_rabi);
    t=0;
    while (t<15){
        t_E = lsr_ptr->GetFieldE({0,0,0}, t).Y();
        t_rabi = solver->GetRabiFreq(t).real();
        t2->Fill();
        t+=0.01;
    }


    ff->Write();
}
