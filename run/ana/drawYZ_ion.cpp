//
// Created by Meng Lv on 2024/10/23.
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
#include "TExec.h"
#include "TLatex.h"
#include <vector>


TH2F *h2;
TH1D * projh2X;
TH1D * projh2Y;
TPad *right_pad, *top_pad;

void drawYZ_ion()
{

    TFile *ff = TFile::Open("data/OBEtest15.root", "READ");
    TTree *tt = (TTree*) ff->Get("obe");

    int n = tt->Draw("y:z:LastRho_ion", "", "goff");
    TGraph2D *gg = new TGraph2D(n, tt->GetV1(), tt->GetV2(), tt->GetV3());
    gg->SetTitle(";y [mm];z [mm]; #rho_{ion}");
    gg->SetNpx(100);
    gg->SetNpy(100);
    auto h2 = gg->GetHistogram();
    gg->GetZaxis()->SetTitleOffset(0.4);
    gg->GetZaxis()->SetTitleSize(0.04);
    gg->GetXaxis()->SetTitleOffset(0.7);
    gg->GetXaxis()->SetTitleSize(0.05);
    gg->GetYaxis()->SetTitleOffset(0.7);
    gg->GetYaxis()->SetTitleSize(0.05);

    auto c1 = new TCanvas("c1", "c1",1200,900);
    gStyle->SetOptStat(0);

    TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.55,0.6);
    center_pad->Draw();

    right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
    right_pad->Draw();

    top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.55,1.0);
    top_pad->Draw();

    projh2X = h2->ProjectionX();
    projh2Y = h2->ProjectionY();

    center_pad->cd();
    gStyle->SetPalette(1);
    h2->Draw("colz");

    top_pad->cd();
    projh2X->SetFillColor(kBlue+1);
    projh2X->GetXaxis()->SetLabelSize(0);
    projh2X->GetXaxis()->SetTitleSize(0);
    projh2X->GetYaxis()->SetLabelSize(0);
    projh2X->GetYaxis()->SetTitleSize(0);
    projh2X->Draw("bar hist");

    right_pad->cd();
    projh2Y->SetFillColor(kBlue-2);
    projh2Y->GetXaxis()->SetLabelSize(0);
    projh2Y->GetXaxis()->SetTitleSize(0);
    projh2Y->GetYaxis()->SetLabelSize(0);
    projh2Y->GetYaxis()->SetTitleSize(0);
    projh2Y->Draw("hbar hist");

    c1->cd();
    TLatex t;
    t.SetTextFont(42);
    t.SetTextSize(0.02);
//    t.DrawLatex(0.6,0.88,"This example demonstrates how to display");
//    t.DrawLatex(0.6,0.85,"a histogram and its two projections.");
//
//    auto ex = new TExec("zoom","ZoomExec()");
//    h2->GetListOfFunctions()->Add(ex);
}

void ZoomExec()
{
    int xfirst = h2->GetXaxis()->GetFirst();
    int xlast = h2->GetXaxis()->GetLast();
    double xmin = h2->GetXaxis()->GetBinLowEdge(xfirst);
    double xmax = h2->GetXaxis()->GetBinUpEdge(xlast);
    projh2X->GetXaxis()->SetRangeUser(xmin, xmax);
    top_pad->Modified();

    int yfirst = h2->GetYaxis()->GetFirst();
    int ylast = h2->GetYaxis()->GetLast();
    double ymin = h2->GetYaxis()->GetBinLowEdge(yfirst);
    double ymax = h2->GetYaxis()->GetBinUpEdge(ylast);
    projh2Y->GetXaxis()->SetRangeUser(ymin, ymax);
    right_pad->Modified();
}