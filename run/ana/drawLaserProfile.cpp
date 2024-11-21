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

void drawLaserProfile() {
    TFile *ff = TFile::Open("../data/test17.root", "READ");
    TTree *tt = (TTree *) ff->Get("pars");

    TCanvas cc("cc", "cc", 600, 600);

    int n = tt->Draw("x:y:Intensity", "z >= 0 && z < 0.1");
    TGraph2D *gg = new TGraph2D(n, tt->GetV1(), tt->GetV2(), tt->GetV3());
    gStyle->SetPalette(1);
    gg->SetTitle(";x [mm];y [mm]; Intensity [W/cm^{2}]");
    gg->GetXaxis()->SetTitleOffset(1.7);
    gg->GetYaxis()->SetTitleOffset(1.7);
    gg->GetZaxis()->SetTitleOffset(1.2);
    gg->SetNpx(100);
    gg->SetNpy(100);
    TH2D *hh = gg->GetHistogram();
    hh->Draw("colz");
    cc.SaveAs("plots/angle/17_xy.pdf");
    delete gg;

    n = tt->Draw("y:z:Intensity", "x >= 0 && x < 0.1");
    gg = new TGraph2D(n, tt->GetV1(), tt->GetV2(), tt->GetV3());
    gStyle->SetPalette(1);
    gg->SetTitle(";y [mm];z [mm]; Intensity [W/cm^{2}]");
    gg->GetXaxis()->SetTitleOffset(1.7);
    gg->GetYaxis()->SetTitleOffset(1.7);
    gg->GetZaxis()->SetTitleOffset(1.2);
    gg->SetNpx(100);
    gg->SetNpy(100);
    hh = gg->GetHistogram();
    hh->Draw("colz");
    cc.SaveAs("plots/angle/17_yz.pdf");
    delete gg;

    n = tt->Draw("z:x:Intensity", "y>=0 && y < 0.1");
    gg = new TGraph2D(n, tt->GetV1(), tt->GetV2(), tt->GetV3());
    gStyle->SetPalette(1);
    gg->SetTitle(";z [mm];x [mm]; Intensity [W/cm^{2}]");
    gg->GetXaxis()->SetTitleOffset(1.7);
    gg->GetYaxis()->SetTitleOffset(1.7);
    gg->GetZaxis()->SetTitleOffset(1.2);
    gg->SetNpx(100);
    gg->SetNpy(100);
    hh = gg->GetHistogram();
    hh->Draw("colz");
    cc.SaveAs("plots/angle/17_zx.pdf");
}
