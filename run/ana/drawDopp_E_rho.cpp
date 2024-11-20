//
// Created by Meng Lv on 2024/10/19.
//
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TStyle.h"
#include <vector>

void drawDopp_E_rho() {
    TFile *ff = TFile::Open("data/OBEtest00.root", "READ");
    TTree *tt = (TTree *) ff->Get("obe");

    int n = tt->Draw("DoppFreq:PeakE^2/3.77:LastRho_ion", "");
    TGraph2D *gg = new TGraph2D(n, tt->GetV1(), tt->GetV2(), tt->GetV3());
    gStyle->SetPalette(1);
    gg->SetTitle(";DoppFreq[GHz];Peak Intensity [W/cm^{2}];#rho_{ion}");
    gg->GetXaxis()->SetTitleOffset(1.7);
    gg->GetYaxis()->SetTitleOffset(1.7);
    gg->GetZaxis()->SetTitleOffset(1.2);
    TCanvas c1("c1", "c1", 800, 600);
    gg->SetMarkerStyle(20);
    gg->Draw("pcol");
    c1.SaveAs("plots/Dopp_E_rho_scatter.pdf");
    gg->SetNpx(100);
    gg->SetNpy(100);
    TH2D *hh = gg->GetHistogram();
    TCanvas cc("cc", "cc", 800, 600);
    hh->Draw("colz");
    cc.SaveAs("plots/Dopp_E_rho.pdf");
}