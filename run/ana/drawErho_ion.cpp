//
// Created by Meng Lv on 2024/10/17.
//

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include <vector>

void drawErho_ion(){
    TFile *ff = TFile::Open("data/OBEtest00.root", "READ");
    TTree *tt = (TTree*) ff->Get("obe");
    Int_t eventID, step_n;
    Double_t x, y, z, vx, vy, vz, dopp_freq, peak_E, last_rho_gg,
             last_rho_ee, last_rho_ion;

    std::vector<Double_t> *t=0;
    std::vector<Double_t> *E_field=0;
    std::vector<Double_t> *rabi_freq=0;
    std::vector<Double_t> *rho_gg=0;
    std::vector<Double_t> *rho_ee=0;
    std::vector<Double_t> *rho_ge_r=0;
    std::vector<Double_t> *rho_ge_i=0;
    std::vector<Double_t> *rho_ion=0;

    // Set branch address
    {
    tt->SetBranchAddress("EventID", &eventID);
    tt->SetBranchAddress("x", &x);
    tt->SetBranchAddress("y", &y);
    tt->SetBranchAddress("z", &z);
    tt->SetBranchAddress("vx", &vx);
    tt->SetBranchAddress("vy", &vy);
    tt->SetBranchAddress("vz", &vz);
    tt->SetBranchAddress("DoppFreq", &dopp_freq);
    tt->SetBranchAddress("PeakE", &peak_E);
    tt->SetBranchAddress("Step_n", &step_n);
    tt->SetBranchAddress("t", &t);
    tt->SetBranchAddress("RabiFreq", &rabi_freq);
    tt->SetBranchAddress("EField", &E_field);
    tt->SetBranchAddress("rho_gg", &rho_gg);
    tt->SetBranchAddress("rho_ee", &rho_ee);
    tt->SetBranchAddress("rho_ge_r", &rho_ge_r);
    tt->SetBranchAddress("rho_ge_i", &rho_ge_i);
    tt->SetBranchAddress("rho_ion", &rho_ion);
    tt->SetBranchAddress("LastRho_gg", &last_rho_gg);
    tt->SetBranchAddress("LastRho_ee", &last_rho_ee);
    tt->SetBranchAddress("LastRho_ion", &last_rho_ion);}

    tt->GetEntry(0);

    TCanvas *c1 = new TCanvas("c1","hists with different scales",800,600);

    //create/fill draw h1
    gStyle->SetOptStat(kFALSE);
    TGraph *g_ion = new TGraph(t->size(), t->data(), rho_ion->data());
    g_ion->SetTitle(";t [ns]; #rho_{ion}");
    g_ion->GetXaxis()->SetTitleSize(0.05);
    g_ion->GetXaxis()->SetTitleOffset(0.9);
    g_ion->GetYaxis()->SetTitleSize(0.05);
    g_ion->GetYaxis()->SetTitleOffset(0.9);
    g_ion->Draw("apl");
    c1->Update();

    Float_t rightmax = 1.1*(*std::max_element(E_field->begin(), E_field->end()));
    Float_t scale = gPad->GetUymax()/rightmax;
    std::vector<Double_t> E_temp(*E_field);
    // scale the plot
    std::for_each(E_temp.begin(), E_temp.end(), [=](Double_t &el){el *= scale;});
    TGraph *g_E = new TGraph(t->size(), t->data(), E_temp.data());
    g_E->SetMarkerColor(kRed);
    g_E->SetLineColor(kRed);
    g_E->Draw("pl same");

    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetTitle("E [V/mm]");
    axis->SetTitleColor(kRed);
    axis->Draw();


}
