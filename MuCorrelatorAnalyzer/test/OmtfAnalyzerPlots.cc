/*
 * OmtfAnalyzerPlots.cc
 *
 *  Created on: Nov 9, 2018
 *      Author: kbunkow
 */


#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include <sstream>
#include <string>
#include <vector>

using namespace std;

TCanvas* canvasCompare = new TCanvas("canvasCompare", "canvasCompare", 1200, 800);

void makeEfficiencyPlots(const char* nameLegend, int color, int ptCut, const char* rootFileName);

TLegend* legend = new TLegend(0.3,0.3,0.9,0.6);

bool compareFirst = true;

int OmtfAnalyzerPlots() {
  gStyle->SetOptStat(0);

  canvasCompare->Divide(2, 2);
  int pt_cut = 20;

  ostringstream ostr;

  //drawEff("HighDTQuality gb1_fin0", kBlack, pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb1_fin0_withHighEta/omtfTTAnalysis1_31kEvents.root");
  //drawEff("HighDTQuality gb2",      kBlue,  pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb2_fin1_noHighEta/omtfTTAnalysis1_31kEv_Ana1.root");
  //drawEff("HighDTQuality gb3",      kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb3_fin1_noHighEta/omtfTTAnalysis1_31kEv_highDtQ.root");

  //drawEff("HighDTQuality gb3",      kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb3_fin1_noHighEta/omtfTTAnalysis1_2kEv_hihgDtQ.root");
  //drawEff("HighDTQuality gb3 2",      kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb3_fin1_noHighEta/omtfTTAnalysis1_2kEv_hihgDtQ_2.root");
  //drawEff("HighDTQuality gb3 loose",    kGreen,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_looseMatch_gb3_fin1_highDtQ/results/omtfAnalysisLooseMatch1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 3 stubs",    kGreen,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_3TrackStubs_gb3_fin1_highDtQ/results/omtfAnalysis1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 4 stubs",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_4Stubs_tightMatch_gb3_fin1_highDtQ/results/omtfAnalysis1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 4 stubs, loose match",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_4Stubs_looseMatch_gb3_fin1_highDtQ/results/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch

  //drawEff("HighDTQuality gb3 4 stubs new",    kGreen,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch

  //drawEff("HighDTQuality gb3 4 stubs, loose match, single mu",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_4Stubs_looseMatch_gb3_fin1_highDtQ_singleMuPu200/results/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 4 stubs, loose match, QCD",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_4Stubs_looseMatch_gb3_fin1_highDtQ_QCD/results/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 4 stubs, loose match, SingleNeutrino",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_4Stubs_looseMatch_gb3_fin1_highDtQ_SingleNeutrino/results/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch
  //drawEff("HighDTQuality gb3 4 stubs, loose match, Higgs",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/crab/crab_omtf_MC_analysis_HToZZTo4L_Pu140/results/omtfTTAnalysis1.root");//omtfAnalysisLooseMatch

  //drawEff("singleMuPu200",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1_neutrinoPu200.root"); //muCorrelatorTTAnalysis1.root
  //drawEff("singleMuPu300",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1_neutrinoPu300.root"); //muCorrelatorTTAnalysis1.root

  //drawEff("singleMu",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1.root");
  //drawEff("singleMu",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200/results/muCorrelatorTTAnalysis1.root");//omtfAnalysisLooseMatch

  //drawEff("SingleNeutrino_PU200 noQualCut",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_newPdfBin_noQualCut/results/muCorrelatorTTAnalysis1.root"); //muCorrelatorTTAnalysis1.root
  //drawEff("SingleNeutrino_PU200 withQulaCut",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino/results/muCorrelatorTTAnalysis1.root"); //muCorrelatorTTAnalysis1.root
  //drawEff("SingleNeutrino_PU200 sigma 1.3",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_sigma1p3/results/muCorrelatorTTAnalysis1.root"); //muCorrelatorTTAnalysis1.root
  //drawEff("SingleNeutrino_PU200 sigma 1.3 plus 3BX",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_sigma1p3_plus3BX/results/muCorrelatorTTAnalysis1.root"); //muCorrelatorTTAnalysis1.root

  //drawEff("all DT Quality", kRed,    pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonOverlap/test/results_gb2_fin1_noHighEta/omtfTTAnalysis1_31kEv_allDtQ.root");

  //drawEff("SingleMu_PU200_newPdfBin_withQualCut",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200_newPdfBin_withQualCut/results/muCorrelatorTTAnalysis1.root");
  //drawEff("SingleMu_PU200_newPdfBin_withQualCut",    kBlue,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200_1_newPdfBin_withQualCut_newAnalyzer/results/muCorrelatorTTAnalysis1.root");
  //drawEff("SingleMu_PU200 sigma 1.2",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200_sigma1p2/results/muCorrelatorTTAnalysis1.root");
  //drawEff("SingleMu_PU200 sigma 1.3",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200_sigma1p3/results/muCorrelatorTTAnalysis1.root");
  //drawEff("HSCP",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1_HSCP.root");

  makeEfficiencyPlots("with trigger rules",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1_withTR.root");
  makeEfficiencyPlots("no trigger rules",    kRed,   pt_cut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1_noTR.root");



/*
  c0->cd();
  legend->Draw();
  c1->cd();
  legend->Draw();
*/

  return 0;

}


void makeEfficiencyPlots(const char* nameLegend, int color, int ptCut, const char* rootFileName)
{

  string canvasName = string("canvas1_") + nameLegend + string("_eff");
  TCanvas* canvas1Eff = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvas1Eff->Divide(3, 3);
  canvas1Eff->SetGrid();

  canvasName = string("canvas1_") + nameLegend + string("_rate");
  TCanvas* canvas1Rate = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvas1Rate->Divide(2, 2);
  canvas1Rate->SetGrid();

  canvasName = string("canvas1_") + nameLegend + string("_timing");
  TCanvas* canvas1Timing= new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvas1Timing->Divide(2, 2);

  TFile * file = new TFile(rootFileName);
  file->ls();

  TDirectory* curDir = file;
  file->cd("omtfTTAnalyzer");
  file->ls();

  //if(string(nameLegend).find("gb3") != string::npos)
  {
    curDir = (TDirectory*)file->Get("omtfTTAnalyzer");
  }

  TH1D* gpMuonPt   = (TH1D*)(curDir->Get("gpMuonPt"));
  TH1D* gpMuonEta  = (TH1D*)curDir->Get("gpMuonEta");
  TH1D* gpMuonEta_ptGen20GeV  = (TH1D*)curDir->Get("gpMuonEta_ptGen20GeV");
  TH1D* gpMuonPhi  = (TH1D*)curDir->Get("gpMuonPhi");

  TH1D* ttMuonPt  = (TH1D*)curDir->Get("ttMuonPt");
  TH1D* ttMuonEta = (TH1D*)curDir->Get("ttMuonEta");
  TH1D* ttMuonEta_ptGen20GeV_ptTT18Gev = (TH1D*)curDir->Get("ttMuonEta_ptGen20GeV_ptTT18Gev");
  TH1D* ttMuonPhi = (TH1D*)curDir->Get("ttMuonPhi");

  TH1D* omtfTTTrackPt  = (TH1D*)curDir->Get("omtfTTTrackPt");
  TH1D* omtfTTTrackEta = (TH1D*)curDir->Get("omtfTTTrackEta");
  TH1D* omtfTTTrackPhi = (TH1D*)curDir->Get("omtfTTTrackPhi");

 //nominator
  TH1D* omtfAndTtMuonPt  = (TH1D*)curDir->Get("omtfAndTtMuonPt");
  TH1D* omtfAndTtMuonEta = (TH1D*)curDir->Get("omtfAndTtMuonEta");
  TH1D* omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev = (TH1D*)curDir->Get("omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev");
  TH1D* omtfAndTtMuonPhi = (TH1D*)curDir->Get("omtfAndTtMuonPhi");

  TH1D* lostTtMuonPt  = (TH1D*)curDir->Get("lostTtMuonPt");
  TH1D* lostTtMuonEta = (TH1D*)curDir->Get("lostTtMuonEta");
  TH1D* lostTtMuonPhi = (TH1D*)curDir->Get("lostTtMuonPhi");

  if(gpMuonPt == nullptr) {
    cout<<" gpMuonPt == nullptr"<<" "<<curDir->Get("gpMuonPt")<<endl;
    return;
  }

  const int nIpt = 27;
  double lower[nIpt + 1] = {0., 4., 4.5, 5., 6., 7., 8., 10., 12., 14., 16., 18.,
      20., 22., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120., 140., 500.};
  for (int i = 0; i <= nIpt; i++){
    lower[i] = lower[i] - 0.5;
  }


  const int orgPtBinCnt = 1000;
  const int orgPtBinMax = 500;

  vector<double> ptBinsVec;

  double bin = 0;
  while(bin <= orgPtBinMax) {
    ptBinsVec.push_back(bin);
    //cout<<"bin "<<bin<<endl;
    //bin+= 1.;

    if(bin < 10) {
      bin+= 0.5;
    }
    else if(bin < 20) {
      bin+= 1.;
    }
    else //if(bin <= 30)
    {
      bin+= 1.;
    }
  }

  TH1* gpMuonPtRebined = gpMuonPt->Rebin(ptBinsVec.size()-1, "gpMuonPtRebined", ptBinsVec.data());

  TH1* ttMuonPtRebined = ttMuonPt->Rebin(ptBinsVec.size()-1, "ttMuonPtRebined", ptBinsVec.data());

  TH1* omtfAndTtMuonPtRebined = omtfAndTtMuonPt->Rebin(ptBinsVec.size()-1, "omtfMuonPtRebined", ptBinsVec.data());

  TH1* omtfMuonPtRebined = omtfTTTrackPt->Rebin(ptBinsVec.size()-1, "omtfMuonPtRebined", ptBinsVec.data());

  canvas1Eff->cd(1);
  gpMuonPt->SetLineColor(kBlue);
  gpMuonPt->Draw("hist");

  gpMuonPtRebined->SetLineColor(kRed);
  gpMuonPtRebined->Draw("samehist");

  gpMuonPt->GetXaxis()->SetRangeUser(0, 100);



  canvas1Eff->cd(2);
  ttMuonPtRebined->Divide(gpMuonPtRebined);
  ttMuonPtRebined->SetTitle("ttTrack efficiency");
  ttMuonPtRebined->Draw("hist");

  omtfAndTtMuonPtRebined->Divide(gpMuonPtRebined);
  omtfAndTtMuonPtRebined->SetLineColor(color);
  omtfAndTtMuonPtRebined->SetTitle("ttTrack and OMTF efficiency");
  omtfAndTtMuonPtRebined->Draw("same");
  ttMuonPtRebined->GetXaxis()->SetRangeUser(0, 100);

/////////////////////////////////
  canvas1Eff->cd(3);
  TH1* ttMuonPtRebined2 = ttMuonPt->Rebin(ptBinsVec.size()-1, "ttMuonPtRebined", ptBinsVec.data());
  TH1* omtfAndTtMuonPtRebined2 = omtfAndTtMuonPt->Rebin(ptBinsVec.size()-1, "omtfMuonPtRebined", ptBinsVec.data());
  omtfAndTtMuonPtRebined2->Divide(ttMuonPtRebined2);
  omtfAndTtMuonPtRebined2->SetLineColor(color);
  omtfAndTtMuonPtRebined2->GetXaxis()->SetRangeUser(0, 100);
  omtfAndTtMuonPtRebined2->SetTitle("OMTF tagging efficiency");
  omtfAndTtMuonPtRebined2->Draw("hist");


  double meanEff = 0;
  int binCnt = 0;
  for(int iBin = 1; iBin <= omtfAndTtMuonPtRebined2->GetNbinsX(); iBin++) {
    if(omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) == 100) {
      break;
    }

    double eff = omtfAndTtMuonPtRebined2->GetBinContent(iBin);
    //cout<<iBin<<" "<<omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) <<" "<<eff<<endl;

    if(omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) > 20) {
      meanEff += eff;
      binCnt++;
    }

  }
  cout<<"mean OMTF tagging efficiency "<<meanEff/binCnt<<endl;



  //ptCut = 18+1;
  ptCut = 1+1;

  canvas1Eff->cd(4);


  TH2I* ptGenPtTTMuon= (TH2I*)curDir->Get("ptGenPtTTMuon");

  TH1D* ptGenPtTTMuonNom = ptGenPtTTMuon->ProjectionX("ptGenPtTTMuonNom", ptCut, -1);
  TH1D* ptGenPtTTMuonDenom = ptGenPtTTMuon->ProjectionX("ptGenPtTTMuonDenom", 0, -1);

  ptGenPtTTMuonNom->SetTitle( ("ttTrack efficiency, pT cut = " + to_string(ptCut -1) + " GeV").c_str());
  ptGenPtTTMuonNom->Divide(ptGenPtTTMuonDenom);
  ptGenPtTTMuonNom->GetYaxis()->SetTitle("efficiency");
  ptGenPtTTMuonNom->SetLineColor(kBlue);
  ptGenPtTTMuonNom->Draw("hist");


  TH2I* ptGenPtOMtfMuon= (TH2I*)curDir->Get("ptGenPtOMtfMuon");
  TH1D* ptGenPtOMtfMuonNom = ptGenPtOMtfMuon->ProjectionX("ptGenPtOMtfMuonNom", ptCut, -1);
  TH1D* ptGenPtOMtfMuonDenom = ptGenPtOMtfMuon->ProjectionX("ptGenPtOMtfMuonDenom", 0, -1);

  //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
  ptGenPtOMtfMuonNom->Divide(ptGenPtTTMuonDenom);
  ptGenPtOMtfMuonNom->SetTitle( ("ttTrack and OMTF efficiency, pT cut = " + to_string(ptCut -1) + "GeV").c_str() );
  ptGenPtOMtfMuonNom->SetLineColor(kRed);
  ptGenPtOMtfMuonNom->Draw("samehist");

  canvas1Eff->cd(7);

  if(curDir->Get("ptGenPtTTMuonEv0")) {
    TH2I* ptGenPtTTMuonEv0= (TH2I*)curDir->Get("ptGenPtTTMuonEv0");

    TH1D* ptGenPtTTMuonNomEv0 = ptGenPtTTMuonEv0->ProjectionX("ptGenPtTTMuonNomEv0", ptCut, -1);
    TH1D* ptGenPtTTMuonDenomEv0 = ptGenPtTTMuonEv0->ProjectionX("ptGenPtTTMuonDenomEv0", 0, -1);

    ptGenPtTTMuonNomEv0->SetTitle( ("ttTrack efficiency, Event 0, pT cut = " + to_string(ptCut -1) + " GeV").c_str() );
    ptGenPtTTMuonNomEv0->Divide(ptGenPtTTMuonDenomEv0);
    ptGenPtTTMuonNomEv0->GetYaxis()->SetTitle("efficiency");
    ptGenPtTTMuonNomEv0->SetLineColor(kBlue);
    ptGenPtTTMuonNomEv0->Draw("hist");


    TH2I* ptGenPtOMtfMuonEv0 = (TH2I*)curDir->Get("ptGenPtOMtfMuonEv0");
    TH1D* ptGenPtOMtfMuonNomEv0 = ptGenPtOMtfMuonEv0->ProjectionX("ptGenPtOMtfMuonNomEv0", ptCut, -1);
    TH1D* ptGenPtOMtfMuonDenomEv0 = ptGenPtOMtfMuonEv0->ProjectionX("ptGenPtOMtfMuonDenomEv0", 0, -1);

    //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
    ptGenPtOMtfMuonNomEv0->Divide(ptGenPtTTMuonDenomEv0);
    ptGenPtOMtfMuonNomEv0->SetTitle( ("ttTrack and OMTF efficiency, Event 0, pT cut = " + to_string(ptCut -1 ) + " GeV").c_str() );
    ptGenPtOMtfMuonNomEv0->SetLineColor(kRed);
    ptGenPtOMtfMuonNomEv0->Draw("samehist");
  }

  canvas1Eff->cd(8);
  if(curDir->Get("ptGenPtTTMuonEvPu")) {
    TH2I* ptGenPtTTMuonEvPu= (TH2I*)curDir->Get("ptGenPtTTMuonEvPu");

    TH1D* ptGenPtTTMuonNomEvPu = ptGenPtTTMuonEvPu->ProjectionX("ptGenPtTTMuonNomEvPu", ptCut, -1);
    TH1D* ptGenPtTTMuonDenomEvPu = ptGenPtTTMuonEvPu->ProjectionX("ptGenPtTTMuonDenomEvPu", 0, -1);

    ptGenPtTTMuonNomEvPu->SetTitle( ("ttTrack efficiency, PU Events, pT cut = " + to_string(ptCut -1) + " GeV").c_str() );
    ptGenPtTTMuonNomEvPu->Divide(ptGenPtTTMuonDenomEvPu);
    ptGenPtTTMuonNomEvPu->GetYaxis()->SetTitle("efficiency");
    ptGenPtTTMuonNomEvPu->SetLineColor(kBlue);
    ptGenPtTTMuonNomEvPu->Draw("hist");


    TH2I* ptGenPtOMtfMuonEvPu = (TH2I*)curDir->Get("ptGenPtOMtfMuonEvPu");
    TH1D* ptGenPtOMtfMuonNomEvPu = ptGenPtOMtfMuonEvPu->ProjectionX("ptGenPtOMtfMuonNomEvPu", ptCut, -1);
    TH1D* ptGenPtOMtfMuonDenomEvPu = ptGenPtOMtfMuonEvPu->ProjectionX("ptGenPtOMtfMuonDenomEvPu", 0, -1);

    //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
    ptGenPtOMtfMuonNomEvPu->Divide(ptGenPtTTMuonDenomEvPu);
    ptGenPtOMtfMuonNomEvPu->SetTitle( ("ttTrack and OMTF efficiency, PU Events, pT cut = " + to_string(ptCut - 1) + " GeV").c_str() );
    ptGenPtOMtfMuonNomEvPu->SetLineColor(kRed);
    ptGenPtOMtfMuonNomEvPu->Draw("samehist");
  }


////////////////////////////////
  canvas1Rate->cd(1);
  canvas1Rate->cd(1)->SetLogy();
  canvas1Rate->cd(1)->SetGridx();
  canvas1Rate->cd(1)->SetGridy();

  TH1D* gpPerEvent = (TH1D*)curDir->Get("gpPerEvent");
  int eventsCnt = gpPerEvent->GetEntries();
  cout<<"eventsCnt: "<<eventsCnt<<endl;

  omtfTTTrackPt = (TH1D*)curDir->Get("omtfTTTrackPt");
  omtfTTTrackPt->SetLineColor(kBlack);
  omtfTTTrackPt->Draw("hist");
  omtfTTTrackPt->GetXaxis()->SetRangeUser(0, 100);
  omtfTTTrackPt->GetYaxis()->SetRangeUser(0.1, 10000);

  TH1D* omtfTTTrackMuonPt = (TH1D*)curDir->Get("omtfTTTrackMuonPt");
  omtfTTTrackMuonPt->SetLineColor(kBlue);
  omtfTTTrackMuonPt->Draw("histsame");

  TH1D* fakeTrackOmtfPt = (TH1D*)curDir->Get("fakeTrackOmtfPt");
  fakeTrackOmtfPt->SetLineColor(kMagenta);
  //fakeTrackOmtfPt->SetFillColor(kMagenta);
  fakeTrackOmtfPt->Draw("histsame");

  TH1D* wrongTagOmtfPt = (TH1D*)curDir->Get("wrongTagOmtfPt");
  wrongTagOmtfPt->SetLineColor(kGreen);
  wrongTagOmtfPt->Draw("histsame");

  ///////
  canvas1Rate->cd(2);
  canvas1Rate->cd(2)->SetLogy();
  canvas1Rate->cd(2)->SetGridx();
  canvas1Rate->cd(2)->SetGridy();

  double lhcFillingRatio = 2760./3564.;
  cout<<"lhcFillingRatio "<<lhcFillingRatio<<endl;
  double scale = 1./eventsCnt * 40000000 * lhcFillingRatio;
  TH1* omtfTTTrackPt_rate = omtfTTTrackPt->GetCumulative(false, "_rate");
  omtfTTTrackPt_rate->Scale(scale);

  omtfTTTrackPt_rate->SetLineColor(kBlack);
  omtfTTTrackPt_rate->Draw("hist");
  omtfTTTrackPt_rate->GetXaxis()->SetRangeUser(0, 100);
  omtfTTTrackPt_rate->GetYaxis()->SetRangeUser(100, 50000000);

  TH1* fakeTrackOmtfPt_rate = fakeTrackOmtfPt->GetCumulative(false, "_rate");
  fakeTrackOmtfPt_rate->Scale(scale);
  fakeTrackOmtfPt_rate->SetLineColor(kMagenta);
  //fakeTrackOmtfPt_rate->SetFillColor(kMagenta);
  fakeTrackOmtfPt_rate->Draw("histsame");

  TH1* wrongTagOmtfPt_rate = wrongTagOmtfPt->GetCumulative(false, "_rate");
  wrongTagOmtfPt_rate->Scale(scale);
  wrongTagOmtfPt_rate->SetLineColor(kGreen);
  wrongTagOmtfPt_rate->Draw("histsame");


  TH1* omtfTTTrackMuonPt_rate = omtfTTTrackMuonPt->GetCumulative(false, "_rate");
  omtfTTTrackMuonPt_rate->Scale(scale);
  omtfTTTrackMuonPt_rate->SetLineColor(kBlue);
  omtfTTTrackMuonPt_rate->Draw("histsame");
  //canvas1->cd(6)->SetLogy();


  /////////////////////
/*
  TH1D* omtfTTTrackPtRate = new TH1D("omtfTTTrackPtRate", "omtfTTTrackPtRate", 200, 0, 100);
  TH1D* omtfTTTrackMuonPtRate = new TH1D("omtfTTTrackMuonPtRate", "omtfTTTrackMuonPtRate", 200, 0, 100);
  TH1D* fakeTrackOmtfPtRate = new TH1D("fakeTrackOmtfPtRate", "fakeTrackOmtfPtRate", 200, 0, 100);
  TH1D* wrongTagOmtfPtRate = new TH1D("wrongTagOmtfPtRate", "wrongTagOmtfPtRate", 200, 0, 100);


  double rateOmtfTTTrackPt = 0;
  double rateOmtfTTTrackMuonPtRate = 0;
  double rateFakeTrackOmtfMuonPt = 0;
  double raeteWrongTagOmtfMuonPt = 0;

  for(int iBin = omtfTTTrackPt->GetNbinsX() +1; iBin >= 0; iBin--) {
    double ptCut = omtfTTTrackPt->GetBinLowEdge(iBin);
    rateOmtfTTTrackPt += omtfTTTrackPt->GetBinContent(iBin);
    rateOmtfTTTrackMuonPtRate += omtfTTTrackPt->GetBinContent(iBin);
    rateFakeTrackOmtfMuonPt += omtfTTTrackPt->GetBinContent(iBin);
    raeteWrongTagOmtfMuonPt += omtfTTTrackPt->GetBinContent(iBin);
  }
*/

  //--------------------------------

  canvas1Eff->cd(6);

  TH1D* ttMuonEtaEff = (TH1D*)ttMuonEta->Clone("ttMuonEtaEff");
  ttMuonEtaEff->Divide(gpMuonEta);
  ttMuonEtaEff->Draw("hist");

  TH1D* omtfAndTtMuonEtaEff = (TH1D*)omtfAndTtMuonEta->Clone("omtfAndTtMuonEtaEff");
  omtfAndTtMuonEtaEff->Divide(gpMuonEta);
  omtfAndTtMuonEtaEff->SetLineColor(kRed);
  omtfAndTtMuonEtaEff->Draw("histsame");

//------------------------
  canvas1Eff->cd(9);
  canvas1Eff->cd(9)->SetGridx();
  canvas1Eff->cd(9)->SetGridy();

  TH1D* ttMuonEta_ptCut20GeV_Eff = (TH1D*)ttMuonEta_ptGen20GeV_ptTT18Gev->Clone("ttMuonEta_ptCut20GeV_Eff");
  ttMuonEta_ptCut20GeV_Eff->Divide(gpMuonEta_ptGen20GeV);
  ttMuonEta_ptCut20GeV_Eff->Draw("hist");

  TH1D* omtfAndTtMuonEta_ptCut20GeV_Eff = (TH1D*)omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev->Clone("omtfAndTtMuonEta_ptCut20GeV_Eff");
  omtfAndTtMuonEta_ptCut20GeV_Eff->Divide(gpMuonEta_ptGen20GeV);
  omtfAndTtMuonEta_ptCut20GeV_Eff->SetLineColor(kRed);
  omtfAndTtMuonEta_ptCut20GeV_Eff->Draw("histsame");

  ttMuonEta_ptCut20GeV_Eff->GetYaxis()->SetRangeUser(0.95, 1.01);

//------------------------



  canvasCompare->cd(1);
  canvasCompare->cd(1)->SetGridx();
  canvasCompare->cd(1)->SetGridy();
  if(compareFirst) {
    omtfAndTtMuonPtRebined2->Draw("hist");
  }
  else {
    omtfAndTtMuonPtRebined2->Draw("samehist");
  }

  canvasCompare->cd(2);
  canvasCompare->cd(2)->SetGridx();
  canvasCompare->cd(2)->SetGridy();
  if(compareFirst) {
    ttMuonPtRebined->SetLineColor(kBlack);
    ttMuonPtRebined->Draw("hist");
    ptGenPtOMtfMuonNom->Draw("samehist");
  }
  else {
    ptGenPtOMtfMuonNom->SetLineColor(color);
    ptGenPtOMtfMuonNom->Draw("samehist");
  }


/*  canvasCompare->cd(3);
  canvasCompare->cd(3)->SetLogy();
  canvasCompare->cd(3)->SetGridx();
  canvasCompare->cd(3)->SetGridy();
  fakeOmtfMuonPt->SetLineColor(color);
  if(first) {
    fakeOmtfMuonPt->GetXaxis()->SetRangeUser(0, 100);
    fakeOmtfMuonPt->SetFillColor(color);
    fakeOmtfMuonPt->Draw("hist");
  }
  else {
    fakeOmtfMuonPt->SetMarkerStyle(3);
    fakeOmtfMuonPt->Draw("samehist");
  }*/


  canvas1Timing->cd(1);
  canvas1Timing->cd(1)->SetGridx();
  canvas1Timing->cd(1)->SetGridy();


  TH2I* betaGenBetaL1Mu = (TH2I*)curDir->Get("betaGenBetaL1Mu");

  int bin0 = betaGenBetaL1Mu->GetYaxis()->FindBin(0.);
  TH1D* l1MuVsBetaGen = betaGenBetaL1Mu->ProjectionX("l1MuVsBetaGen", bin0, -1);
  TH1D* allVsBetaGen = betaGenBetaL1Mu->ProjectionX("allVsBetaGen", -1, -1);
  l1MuVsBetaGen->Divide(allVsBetaGen);
  l1MuVsBetaGen->Draw("hist");

  compareFirst = false;
}


