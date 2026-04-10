// example_analysis.C
//
// Example ROOT macro to read the output of IonizationToyMCSim.
// Run from the run/ana/ directory:
//
//   root -l example_analysis.C
//
// Or with a custom file path:
//
//   root -l 'example_analysis.C("../data/my_output.root")'

#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

void example_analysis(const char* filename = "../data/g2edmIoni_test.root") {

    // -------------------------------------------------------------------------
    // Open the ROOT file and retrieve the tree
    // -------------------------------------------------------------------------
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)f->Get("obe");
    if (!tree) {
        std::cerr << "Error: tree 'obe' not found in file." << std::endl;
        f->Close();
        return;
    }

    std::cout << "Opened: " << filename << std::endl;
    std::cout << "Total events: " << tree->GetEntries() << std::endl;

    // -------------------------------------------------------------------------
    // Declare variables and connect them to tree branches
    // -------------------------------------------------------------------------

    // Per-event scalars (always present)
    Int_t    eventID;
    Double_t x, y, z;            // Initial position [mm]
    Double_t vx, vy, vz;         // Initial velocity [m/s]
    Double_t dopp_freq;           // Peak Doppler frequency [GHz]
    Double_t peak_int_122;        // Peak 122 nm intensity [W/mm^2]
    Double_t peak_int_355;        // Peak 355 nm intensity [W/mm^2]
    Int_t    step_n;              // Number of ODE time steps
    Double_t last_rho_gg;         // Ground state pop. at end of simulation
    Double_t last_rho_ee;         // Excited state pop. at end of simulation
    Double_t last_rho_ion;        // Ionized state pop. at end of simulation
    Int_t    if_ionized;          // Ionization flag (1 = ionized, 0 = not)
    Double_t ioni_time;           // Sampled ionization time [ns] (-1 if not ionized)

    tree->SetBranchAddress("EventID",         &eventID);
    tree->SetBranchAddress("x",               &x);
    tree->SetBranchAddress("y",               &y);
    tree->SetBranchAddress("z",               &z);
    tree->SetBranchAddress("vx",              &vx);
    tree->SetBranchAddress("vy",              &vy);
    tree->SetBranchAddress("vz",              &vz);
    tree->SetBranchAddress("DoppFreq",        &dopp_freq);
    tree->SetBranchAddress("PeakIntensity122",&peak_int_122);
    tree->SetBranchAddress("PeakIntensity355",&peak_int_355);
    tree->SetBranchAddress("Step_n",          &step_n);
    tree->SetBranchAddress("LastRho_gg",      &last_rho_gg);
    tree->SetBranchAddress("LastRho_ee",      &last_rho_ee);
    tree->SetBranchAddress("LastRho_ion",     &last_rho_ion);
    tree->SetBranchAddress("IfIonized",       &if_ionized);
    tree->SetBranchAddress("IoniTime",        &ioni_time);

    // Per-timestep arrays (only present if enabled via RootOutput in the macro)
    // Check whether the branch exists before connecting to avoid errors.
    std::vector<Double_t> *t_arr=0, *rabi_freq=0, *E_field=0;
    std::vector<Double_t> *intensity_122=0, *intensity_355=0, *gamma_ion=0;
    std::vector<Double_t> *rho_gg=0, *rho_ee=0, *rho_ge_r=0, *rho_ge_i=0, *rho_ion=0;

    bool has_t            = tree->GetBranch("t")            != nullptr;
    bool has_rabi         = tree->GetBranch("RabiFreq")     != nullptr;
    bool has_efield       = tree->GetBranch("EField")       != nullptr;
    bool has_int122       = tree->GetBranch("Intensity122")  != nullptr;
    bool has_int355       = tree->GetBranch("Intensity355")  != nullptr;
    bool has_gammaion     = tree->GetBranch("GammaIon")     != nullptr;
    bool has_rho_gg       = tree->GetBranch("rho_gg")       != nullptr;
    bool has_rho_ee       = tree->GetBranch("rho_ee")       != nullptr;
    bool has_rho_ge_r     = tree->GetBranch("rho_ge_r")     != nullptr;
    bool has_rho_ge_i     = tree->GetBranch("rho_ge_i")     != nullptr;
    bool has_rho_ion      = tree->GetBranch("rho_ion")      != nullptr;

    if (has_t)        tree->SetBranchAddress("t",           &t_arr);
    if (has_rabi)     tree->SetBranchAddress("RabiFreq",    &rabi_freq);
    if (has_efield)   tree->SetBranchAddress("EField",      &E_field);
    if (has_int122)   tree->SetBranchAddress("Intensity122", &intensity_122);
    if (has_int355)   tree->SetBranchAddress("Intensity355", &intensity_355);
    if (has_gammaion) tree->SetBranchAddress("GammaIon",    &gamma_ion);
    if (has_rho_gg)   tree->SetBranchAddress("rho_gg",      &rho_gg);
    if (has_rho_ee)   tree->SetBranchAddress("rho_ee",      &rho_ee);
    if (has_rho_ge_r) tree->SetBranchAddress("rho_ge_r",    &rho_ge_r);
    if (has_rho_ge_i) tree->SetBranchAddress("rho_ge_i",    &rho_ge_i);
    if (has_rho_ion)  tree->SetBranchAddress("rho_ion",     &rho_ion);

    // -------------------------------------------------------------------------
    // Event loop
    // -------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    int n_ionized = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // --- Per-event information ---
        std::cout << "Event " << eventID
                  << "  pos=(" << x << ", " << y << ", " << z << ") mm"
                  << "  vel=(" << vx << ", " << vy << ", " << vz << ") m/s"
                  << "  DopplerFreq=" << dopp_freq << " GHz"
                  << "  PeakI_122=" << peak_int_122 << " W/mm^2"
                  << "  PeakI_355=" << peak_int_355 << " W/mm^2"
                  << "  steps=" << step_n
                  << "  LastRho_ion=" << last_rho_ion
                  << "  Ionized=" << if_ionized
                  << "  IoniTime=" << ioni_time << " ns"
                  << std::endl;

        if (if_ionized) n_ionized++;

        // --- Per-timestep information (if arrays were saved) ---
        if (has_t && has_rho_ion && !t_arr->empty()) {
            // Example: print the density matrix evolution at a few time points
            int n_steps = (int)t_arr->size();
            std::cout << "  Time evolution (first/mid/last step):" << std::endl;
            for (int idx : {0, n_steps / 2, n_steps - 1}) {
                std::cout << "    t=" << t_arr->at(idx) << " ns"
                          << "  rho_gg=" << rho_gg->at(idx)
                          << "  rho_ee=" << rho_ee->at(idx)
                          << "  rho_ion=" << rho_ion->at(idx)
                          << std::endl;
            }
        }

        // =====================================================================
        // ADD YOUR ANALYSIS HERE
        // =====================================================================

    }

    // -------------------------------------------------------------------------
    // Summary
    // -------------------------------------------------------------------------
    std::cout << "\n--- Summary ---" << std::endl;
    std::cout << "Total events : " << nEntries << std::endl;
    std::cout << "Ionized      : " << n_ionized << std::endl;
    std::cout << "Ionization probability : " << 100.0 * last_rho_ion << " %" << std::endl;

    f->Close();
}
