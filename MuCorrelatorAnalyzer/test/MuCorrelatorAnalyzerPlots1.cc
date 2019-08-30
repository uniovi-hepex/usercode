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
#include "TEfficiency.h"
#include "THStack.h"

#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <regex>

using namespace std;



void makeEfficiencyPlots(const char* nameLegend, int color, int ptCut, const char* rootFileName);
void makeCandidatesMatchingPlots(TDirectory* omtfTTAnalyzerDir, const char* nameLegend, double eventsCnt);

TCanvas* canvasCompare = new TCanvas("canvasCompare", "canvasCompare", 1200, 800);

bool compareFirst = true;

double lhcFillingRatio = 2748./3564.;

int MuCorrelatorAnalyzerPlots1() {
  gStyle->SetOptStat(0);


  canvasCompare->Divide(2, 2);
  int ptCut = 20;

  //ptCut = 18+1;
  ptCut = 1+1;
  //ptCut = 10 +1;

  ostringstream ostr;

  //makeEfficiencyPlots("hscp",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1HSCP_100Files.root");
  //makeEfficiencyPlots("hscp",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1HSCP.root");
  //makeEfficiencyPlots("hscp CMSSW_10_5_0_pre1",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1HSCP_100FilesUpdatedCuts.root");

  //makeEfficiencyPlots("hscp_noiRPC",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1HSCP_noiRPC.root");

  //makeEfficiencyPlots("hscp CMSSW_10_6_0_pre4",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1HSCP.root");


  //makeEfficiencyPlots("singleMu",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("singleNeutrinoPU200 without_iRPC",    kGreen,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/without_iRPC/crab_muCorr_MC_analysis_SingleNeutrino_PU200_600Files/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("singleNeutrinoPU300_noChi2DofCut",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/without_iRPC/crab_muCorr_MC_analysis_SingleNeutrinoPU300_noChi2DofCut/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("singleNeutrinoPU300",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrinoPU300/results/muCorrelatorTTAnalysis1.root");


  //makeEfficiencyPlots("singleNeutrinoPU200 CMSSW_10_1_7",  kBlue,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU200_600Files/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("singleNeutrinoPU300",  kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU300/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("singleNeutrinoPU200 CMSSW_10_6_0_pre4",  kBlue,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4//src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU200/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("singleNeutrinoPU300 CMSSW_10_6_0_pre4",  kBlue,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4//src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU200/results/muCorrelatorTTAnalysis1.root");


  //makeEfficiencyPlots("SingleMu_PU200",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("SingleMu_PU200 CMSSW_10_6_0_pre4",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleMu_PU200_v0/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("HToZZTo4L CMSSW_10_6_0_pre4 gb2Stub",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_HToZZTo4L_v0_gb2Stub/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("HToZZTo4L CMSSW_10_6_0_pre4 gb1Stub",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_0_pre4/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_HToZZTo4L_v-3_gb1Stub/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("SingleMu_PU200_without_iRPC",    kRed,   ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_1_7/src/L1Trigger/L1TMuonBayes/test/crab/without_iRPC/crab_muCorr_MC_analysis_SingleMu_PU200/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("HToZZTo4L CMSSW_10_6_0_pre4 gb1Stub",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_GluGluHToZZTo4L_NoPU.root");

  //makeEfficiencyPlots("SingleNeutrino_PU200_t5",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU200_v1_t5/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("SingleNeutrino_PU200_t7",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_SingleNeutrino_PU200_v1_t7/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("TTTo2L2Nu_PU200",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_TTTo2L2Nu_v1_t6/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("DYToMuMuorEleEle",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_DYToMuMuorEleEle_v1_t6/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("MuFlatPt_PU200_t6",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_MuFlatPt_PU200_v1_t6/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("MuFlatPt_PU200_t7",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_MuFlatPt_PU200_v1_t7/results/muCorrelatorTTAnalysis1.root");
  //makeEfficiencyPlots("JPsiToMuMu_Pt0to100_PU200",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_muCorr_MC_analysis_JPsiToMuMu_Pt0to100_PU200_v1_t6/results/muCorrelatorTTAnalysis1.root");

  //makeEfficiencyPlots("GluGluHToZZTo4L_NoPU gb4",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_GluGluHToZZTo4L_NoPU_gb4.root");
  //makeEfficiencyPlots("singleMu_sigma_2p8",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_singleMu_sigma_2p8.root");
  //makeEfficiencyPlots("singleMu_sigma_1p4",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_singleMu_sigma_1p4.root");

  //makeEfficiencyPlots("HSCPppstau_M_200_NoPU",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_HSCPppstau_M_200_NoPU.root");
  //makeEfficiencyPlots("JPsiToMuMu_Pt0to100_NoPU gb4",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_JPsiToMuMu_Pt0to100_NoPU_gb4.root");
  makeEfficiencyPlots("JPsiToMuMu_Pt0to100_NoPU",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_JPsiToMuMu_Pt0to100_NoPU.root");

  //makeEfficiencyPlots("JPsiToMuMu_Pt0to100_NoPU gb3",    kRed,       ptCut,  "/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/muCorrelatorTTAnalysis1_JPsiToMuMu_Pt0to100_NoPU_gb3.root");

/*
  c0->cd();
  legend->Draw();
  c1->cd();
  legend->Draw();
*/

  return 0;

}

void makeEffVsBeta(TDirectory* omtfTTAnalyzerDir, const char* nameLegend) {
  string canvasName = string("canvas_") + nameLegend + string("_beta");
  TCanvas* canvasBetaPlost = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvasBetaPlost->cd();
  canvasBetaPlost->Divide(2, 2);
  canvasBetaPlost->cd(1);
  canvasBetaPlost->cd(1)->SetGridx();
  canvasBetaPlost->cd(1)->SetGridy();
  TLegend* legend = new TLegend(0.3,0.3,0.9,0.6);

  TIter next(omtfTTAnalyzerDir->GetListOfKeys());
  TKey* key;
  bool first = true;
  int color = 2;

  TH1D* l1MuVsBetaGenSingleMu = nullptr;

  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();
      if(dirName.find("EfficiencyAnalyser_") != string::npos) {
        TDirectory* subdir = gDirectory;
        cout<<"making beta  plots for "<<dirName<<" color "<<color<<endl;

        TH2I* betaGenBetaL1Mu = (TH2I*)subdir->Get("betaGenBetaL1Mu");

        int bin0 = betaGenBetaL1Mu->GetYaxis()->FindBin(0.);
        TH1D* l1MuVsBetaGen = betaGenBetaL1Mu->ProjectionX( ( string("l1MuVsBetaGen_") + to_string(color) ).c_str(), bin0, -1);
        TH1D* allVsBetaGen = betaGenBetaL1Mu->ProjectionX("allVsBetaGen", -1, -1);

        if(dirName == "EfficiencyAnalyser_SingleMuAlgo20") {
          l1MuVsBetaGenSingleMu = (TH1D*)l1MuVsBetaGen->Clone("l1MuVsBetaGenSingleMu");
        }

        l1MuVsBetaGen->Divide(allVsBetaGen);

        l1MuVsBetaGen->SetLineColor(color);

        canvasBetaPlost->cd(1);

        if(first)
          l1MuVsBetaGen->Draw("hist");
        else
          l1MuVsBetaGen->Draw("histsame");

        cout<<"l1MuVsBetaGen->GetName() "<<l1MuVsBetaGen->GetName()<<endl;

        legend->AddEntry(l1MuVsBetaGen, dirName.erase(0, 19).c_str());

        first = false;
        color++;
        if(color == 5)
          color = 6;
      }
    }
  }

  canvasBetaPlost->cd(2);
  legend->Draw();

  canvasBetaPlost->cd(3)->SetGridx();
  canvasBetaPlost->cd(3)->SetGridy();
  first = true;
  next.Reset();
  color = 2;

  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();
      if(dirName.find("EfficiencyAnalyser_") != string::npos) {
        TDirectory* subdir = gDirectory;
        cout<<"making beta  plots for "<<dirName<<endl;

        TH2I* betaGenBetaL1Mu = (TH2I*)subdir->Get("betaGenBetaL1Mu");

        int bin0 = betaGenBetaL1Mu->GetYaxis()->FindBin(0.);
        TH1D* l1MuVsBetaGen = betaGenBetaL1Mu->ProjectionX("l1MuVsBetaGen2", bin0, -1);
        TH1D* allVsBetaGen = betaGenBetaL1Mu->ProjectionX("allVsBetaGen2", -1, -1);

        if(dirName != "EfficiencyAnalyser_SingleMuAlgo20") {
          TH1D* l1MuVsBetaGenSum = (TH1D*)l1MuVsBetaGen->Clone("l1MuVsBetaGenSum");
          TH1D* allVsBetaGenSum = (TH1D*)allVsBetaGen->Clone("allVsBetaGenSum");

          l1MuVsBetaGenSum->Add(l1MuVsBetaGenSingleMu);

          l1MuVsBetaGenSum->Divide(allVsBetaGenSum);

          l1MuVsBetaGenSum->SetLineColor(color);

          canvasBetaPlost->cd(3);

          l1MuVsBetaGenSum->SetTitle("total efficiency vs. beta gen");
          if(first)
            l1MuVsBetaGenSum->Draw("hist");
          else
            l1MuVsBetaGenSum->Draw("histsame");

          first = false;

          /*if(color == 7) {
            canvasCompare->cd(3);
            if(compareFirst)
              l1MuVsBetaGenSum->Draw("hist");
            else
              l1MuVsBetaGenSum->Draw("histsame");
          }*/

        }

        color++;
        if(color == 5)
          color = 6;
      }
    }
  }

}

void makeEfficiencyPlots(const char* nameLegend, int color, int ptCut, const char* rootFileName)
{

  string canvasName = string("canvas_") + nameLegend + string("_controlPlost");
  TCanvas* canvasControlPlost = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvasControlPlost->Divide(2, 2);
  canvasControlPlost->SetGrid();

/*  canvasName = string("canvas1_") + nameLegend + string("_timing");
  TCanvas* canvas1Timing= new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
  canvas1Timing->Divide(2, 2);*/

  TFile* file = new TFile(rootFileName);
  file->ls();

  TDirectory* omtfTTAnalyzerDir = file;
  file->cd("omtfTTAnalyzer");
  file->ls();

  omtfTTAnalyzerDir = (TDirectory*)file->Get("omtfTTAnalyzer");

  TH1D* gpPerEvent = (TH1D*)omtfTTAnalyzerDir->Get("gpPerEvent");
  double eventsCnt = gpPerEvent->GetEntries();

  makeEffVsBeta(omtfTTAnalyzerDir, nameLegend);

  makeCandidatesMatchingPlots(omtfTTAnalyzerDir, nameLegend, eventsCnt);


  TH1D* gpMuonPt   = (TH1D*)(omtfTTAnalyzerDir->Get("gpMuonPt"));
  TH1D* gpMuonEta  = (TH1D*)omtfTTAnalyzerDir->Get("gpMuonEta");
  TH1D* gpMuonEta_ptGen20GeV  = (TH1D*)omtfTTAnalyzerDir->Get("gpMuonEta_ptGen20GeV");
  TH1D* gpMuonPhi  = (TH1D*)omtfTTAnalyzerDir->Get("gpMuonPhi");

  TH1D* ttMuonPt  = (TH1D*)omtfTTAnalyzerDir->Get("ttMuonPt");
  TH1D* ttMuonEta = (TH1D*)omtfTTAnalyzerDir->Get("ttMuonEta");
  TH1D* ttMuonEta_ptGen20GeV_ptTT18Gev = (TH1D*)omtfTTAnalyzerDir->Get("ttMuonEta_ptGen20GeV_ptTT18Gev");
  TH1D* ttMuonPhi = (TH1D*)omtfTTAnalyzerDir->Get("ttMuonPhi");

  TH1D* muCandPt  = (TH1D*)omtfTTAnalyzerDir->Get("muCandPt");
  TH1D* muCandEta = (TH1D*)omtfTTAnalyzerDir->Get("muCandEta");
  TH1D* muCandPhi = (TH1D*)omtfTTAnalyzerDir->Get("muCandPhi");



  //8888888888888888888888888888
  canvasControlPlost->cd(1);
  canvasControlPlost->cd(1)->SetLogy();
  gpMuonPt->SetLineColor(kBlue);
  gpMuonPt->Draw("hist");

  ttMuonPt->SetLineColor(kRed);
  ttMuonPt->Draw("histsame");

  gpMuonPt->GetXaxis()->SetRangeUser(0, 100);

  canvasControlPlost->cd(2);
  canvasControlPlost->cd(2)->SetLogy();

  cout<<"eventsCnt: "<<eventsCnt<<endl;
  gpPerEvent->Draw("hist");

  //////////////////////////////////////////////////
  canvasControlPlost->cd(3);
  TH2F* ttTrackEta_Pt = (TH2F*)omtfTTAnalyzerDir->Get("ttTrackEta_Pt");
  TH2F* ttMuonEta_Pt = (TH2F*)omtfTTAnalyzerDir->Get("ttMuonEta_Pt");
  ttMuonEta_Pt->Divide(ttTrackEta_Pt);
  ttMuonEta_Pt->SetTitle("ttMuonEta_Pt / ttTrackEta_Pt");
  ttMuonEta_Pt->Draw("colz");

  //_---------------------------------------------------------------------------
  canvasCompare->cd(1);
  canvasCompare->cd(1)->SetLogy();
  canvasCompare->cd(1)->SetGridx();
  canvasCompare->cd(1)->SetGridy();
  if(omtfTTAnalyzerDir->Get("ttTracksPt")) {
    TH1* ttTracksPt  = (TH1*)omtfTTAnalyzerDir->Get("ttTracksPt");
    TH1* ttTracksPt_rateCumul = ttTracksPt->GetCumulative(false, "_rate");
    ttTracksPt_rateCumul->GetXaxis()->SetTitle("ttTrack pt tresh [GeV]");
    ttTracksPt_rateCumul->GetYaxis()->SetTitle("<ttTracks per event>");
    ttTracksPt_rateCumul->Scale(1./eventsCnt);
    ttTracksPt_rateCumul->SetLineColor(color);
    if(compareFirst) {
      ttTracksPt_rateCumul->Draw("hist");
      ttTracksPt_rateCumul->GetYaxis()->SetRangeUser(0.001, 5000);
    }
    else
      ttTracksPt_rateCumul->Draw("histsame");
  }

  canvasCompare->cd(2);
  canvasCompare->cd(2)->SetLogy();
  canvasCompare->cd(2)->SetGridx();
  canvasCompare->cd(2)->SetGridy();
  TH1* ttTracksFakesPt  = (TH1*)omtfTTAnalyzerDir->Get("ttTracksFakesPt");
  TH1* ttTracksFakesPt_rateCumul = ttTracksFakesPt->GetCumulative(false, "_rate");
  ttTracksFakesPt_rateCumul->GetXaxis()->SetTitle("ttTrack pt tresh [GeV]");
  ttTracksFakesPt_rateCumul->GetYaxis()->SetTitle("<ttTracks per event>");
  ttTracksFakesPt_rateCumul->Scale(1./eventsCnt);
  ttTracksFakesPt_rateCumul->SetLineColor(color);
  if(compareFirst) {
    ttTracksFakesPt_rateCumul->Draw("hist");
    ttTracksFakesPt_rateCumul->GetYaxis()->SetRangeUser(0.001, 5000);
  }
  else
    ttTracksFakesPt_rateCumul->Draw("histsame");
//*************************************************** efficiency ********************************************
  //tracking trigger efficiency, reference for all algos
  TH2I* ptGenPtTTMuon= (TH2I*)omtfTTAnalyzerDir->Get("ptGenPtTTMuon");

  TH1D* ptGenPtTTMuonNom = ptGenPtTTMuon->ProjectionX("ptGenPtTTMuonNom", ptCut, -1);
  TH1D* ptGenPtTTMuonDenom = ptGenPtTTMuon->ProjectionX("ptGenPtTTMuonDenom", 0, -1);

  ptGenPtTTMuonNom->SetTitle( ("ttTrack efficiency, pT cut = " + to_string(ptCut -1) + " GeV").c_str());
  ptGenPtTTMuonNom->Divide(ptGenPtTTMuonDenom);
  ptGenPtTTMuonNom->GetYaxis()->SetTitle("efficiency");
  ptGenPtTTMuonNom->SetLineColor(kBlue);


  TH2I* ptGenPtTTMuonEv0= (TH2I*)omtfTTAnalyzerDir->Get("ptGenPtTTMuonEv0");

  TH1D* ptGenPtTTMuonNomEv0 = ptGenPtTTMuonEv0->ProjectionX("ptGenPtTTMuonNomEv0", ptCut, -1);
  TH1D* ptGenPtTTMuonDenomEv0 = ptGenPtTTMuonEv0->ProjectionX("ptGenPtTTMuonDenomEv0", 0, -1);

  ptGenPtTTMuonNomEv0->SetTitle( ("ttTrack efficiency, Event 0, pT cut = " + to_string(ptCut -1) + " GeV").c_str() );
  ptGenPtTTMuonNomEv0->Divide(ptGenPtTTMuonDenomEv0);
  ptGenPtTTMuonNomEv0->GetYaxis()->SetTitle("efficiency");
  ptGenPtTTMuonNomEv0->SetLineColor(kBlue);


  TH2I* ptGenPtTTMuonEvPu= (TH2I*)omtfTTAnalyzerDir->Get("ptGenPtTTMuonEvPu");

  TH1D* ptGenPtTTMuonNomEvPu = ptGenPtTTMuonEvPu->ProjectionX("ptGenPtTTMuonNomEvPu", ptCut, -1);
  TH1D* ptGenPtTTMuonDenomEvPu = ptGenPtTTMuonEvPu->ProjectionX("ptGenPtTTMuonDenomEvPu", 0, -1);

  ptGenPtTTMuonNomEvPu->SetTitle( ("ttTrack efficiency, PU Events, pT cut = " + to_string(ptCut -1) + " GeV").c_str() );
  ptGenPtTTMuonNomEvPu->Divide(ptGenPtTTMuonDenomEvPu);
  ptGenPtTTMuonNomEvPu->GetYaxis()->SetTitle("efficiency");
  ptGenPtTTMuonNomEvPu->SetLineColor(kBlue);


  TH1D* ttMuonEtaEff = (TH1D*)ttMuonEta->Clone("ttMuonEtaEff");
  ttMuonEtaEff->Divide(gpMuonEta);

  //TDirectory* currentDir = gDirectory;
  TIter next(omtfTTAnalyzerDir->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();

      if(dirName.find("EfficiencyAnalyser_") != string::npos) {
        TDirectory* subdir = gDirectory;
        cout<<"making plots for "<<dirName<<endl;
        canvasName = string("canvas_") + nameLegend + "_" + dirName;
        TCanvas* canvas1Eff = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 1200);
        canvas1Eff->Divide(2, 2);

        canvas1Eff->cd(1);
        canvas1Eff->cd(1)->SetGridx();
        canvas1Eff->cd(1)->SetGridy();

        ptGenPtTTMuonNom->Draw("hist");

        TH2I* ptGenPtMuCandMuons= (TH2I*)subdir->Get("ptGenPtMuCandMuons");
        TH1D* ptGenPtMuCandMuonsNom = ptGenPtMuCandMuons->ProjectionX("ptGenPtMuCandMuonsNom", ptCut, -1);
        TH1D* ptGenPtMuCandMuonsDenom = ptGenPtMuCandMuons->ProjectionX("ptGenPtMuCandMuonsDenom", 0, -1);

        //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
        ptGenPtMuCandMuonsNom->Divide(ptGenPtTTMuonDenom);
        ptGenPtMuCandMuonsNom->SetTitle( ("ttTrack and OMTF efficiency, pT cut = " + to_string(ptCut -1) + "GeV").c_str() );
        ptGenPtMuCandMuonsNom->SetLineColor(kRed);
        ptGenPtMuCandMuonsNom->Draw("samehist");

        canvas1Eff->cd(2);
        canvas1Eff->cd(2)->SetGridx();
        canvas1Eff->cd(2)->SetGridy();

        ptGenPtTTMuonNomEv0->Draw("hist");


        TH2I* ptGenPtMuCandMuonsEv0 = (TH2I*)subdir->Get("ptGenPtMuCandMuonsEv0");
        TH1D* ptGenPtMuCandMuonsEv0Nom = ptGenPtMuCandMuonsEv0->ProjectionX("ptGenPtMuCandMuonsEv0Nom", ptCut, -1);
        TH1D* ptGenPtMuCandMuonsEv0Denom = ptGenPtMuCandMuonsEv0->ProjectionX("ptGenPtMuCandMuonsEv0Denom", 0, -1);

        //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
        ptGenPtMuCandMuonsEv0Nom->Divide(ptGenPtTTMuonDenomEv0);
        ptGenPtMuCandMuonsEv0Nom->SetTitle( ("ttTrack and OMTF efficiency, Event 0, pT cut = " + to_string(ptCut -1 ) + " GeV").c_str() );
        ptGenPtMuCandMuonsEv0Nom->SetLineColor(kRed);
        ptGenPtMuCandMuonsEv0Nom->Draw("samehist");

        //--------------------------------
        canvas1Eff->cd(3);
        canvas1Eff->cd(3)->SetGridx();
        canvas1Eff->cd(3)->SetGridy();

        ptGenPtTTMuonNomEvPu->Draw("hist");


        TH2I* ptGenPtMuCandMuonsPu = (TH2I*)subdir->Get("ptGenPtMuCandMuonsPu");
        TH1D* ptGenPtMuCandMuonsPuNom = ptGenPtMuCandMuonsPu->ProjectionX("ptGenPtMuCandMuonsPuNom", ptCut, -1);
        TH1D* ptGenPtMuCandMuonsPuDenom = ptGenPtMuCandMuonsPu->ProjectionX("ptGenPtMuCandMuonsPuDenom", 0, -1);

        //ptGenPtOMtfMuonNom->Divide(ptGenPtOMtfMuonDenom); //TODO!!!! in principle ptGenPtOMtfMuonDenom and ptGenPtTTMuonDenom should be the same
        ptGenPtMuCandMuonsPuNom->Divide(ptGenPtTTMuonDenomEvPu);
        ptGenPtMuCandMuonsPuNom->SetTitle( ("ttTrack and OMTF efficiency, PU Events, pT cut = " + to_string(ptCut - 1) + " GeV").c_str() );
        ptGenPtMuCandMuonsPuNom->SetLineColor(kRed);
        ptGenPtMuCandMuonsPuNom->Draw("samehist");

        //--------------------------------

        double meanEff = 0;
        int binCnt = 0; //TODO

        TH1D* muCandGenEtaMuons_ptGen20GeV_ptTT18GeV = (TH1D*)subdir->Get("muCandGenEtaMuons_ptGen20GeV_ptTT18GeV");

        double ttMuonEff_ptGen20GeV_ptTT18Gev = ttMuonEta_ptGen20GeV_ptTT18Gev->Integral() / gpMuonEta_ptGen20GeV->Integral();
        double muCadnEffEff_ptGen20GeV_ptTT18Gev = muCandGenEtaMuons_ptGen20GeV_ptTT18GeV->Integral() / gpMuonEta_ptGen20GeV->Integral();
        double taggingEff_ptGen20GeV_ptTT18Gev = muCandGenEtaMuons_ptGen20GeV_ptTT18GeV->Integral() / ttMuonEta_ptGen20GeV_ptTT18Gev->Integral();
        /*  for(int iBin = 1; iBin <= omtfAndTtMuonPtRebined2->GetNbinsX(); iBin++) {
            if(omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) == 100) {
              break;
            }

            double eff = omtfAndTtMuonPtRebined2->GetBinContent(iBin);
            //cout<<iBin<<" "<<omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) <<" "<<eff<<endl;

            if(omtfAndTtMuonPtRebined2->GetBinLowEdge(iBin) > 20) {
              meanEff += eff;
              binCnt++;
            }

          }*/
        cout<<"ttMuonEff_ptGen20GeV_ptTT18Gev    "<<ttMuonEff_ptGen20GeV_ptTT18Gev<<endl;
        cout<<"muCadnEffEff_ptGen20GeV_ptTT18Gev "<<muCadnEffEff_ptGen20GeV_ptTT18Gev<<endl;
        cout<<"taggingEff_ptGen20GeV_ptTT18Gev   "<<taggingEff_ptGen20GeV_ptTT18Gev<<endl<<endl;

        //------------------------

        canvas1Eff->cd(4);
        canvas1Eff->cd(4)->SetGridx();
        canvas1Eff->cd(4)->SetGridy();

        if(0) {
          TH1D* ttMuonEta_ptCut20GeV_Eff = (TH1D*)ttMuonEta_ptGen20GeV_ptTT18Gev->Clone("ttMuonEta_ptCut20GeV_Eff");
          ttMuonEta_ptCut20GeV_Eff->Divide(gpMuonEta_ptGen20GeV);
          ttMuonEta_ptCut20GeV_Eff->Draw("hist");

          TH1D* muCandGenEtaMuons_ptGen20GeV_ptTT18GeV_Eff = (TH1D*)muCandGenEtaMuons_ptGen20GeV_ptTT18GeV->Clone("muCandGenEtaMuons_ptGen20GeV_ptTT18GeV_Eff");
          muCandGenEtaMuons_ptGen20GeV_ptTT18GeV_Eff->Divide(gpMuonEta_ptGen20GeV);
          muCandGenEtaMuons_ptGen20GeV_ptTT18GeV_Eff->SetLineColor(kRed);
          muCandGenEtaMuons_ptGen20GeV_ptTT18GeV_Eff->Draw("histsame");

          ttMuonEta_ptCut20GeV_Eff->GetYaxis()->SetRangeUser(0.95, 1.01);
        }
        else {

          TH1D* muCandGenEtaMuons_ptGenCutm2_ptTTCut = (TH1D*)subdir->Get("muCandGenEtaMuons_ptGenCutm2_ptTTCut");
          TH1D* gpMuonGenEtaMuons_ptGenCutm2_ptTTCut = (TH1D*)subdir->Get("gpMuonGenEtaMuons_ptGenCutm2_ptTTCut");


          TH1D* muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff = (TH1D*)muCandGenEtaMuons_ptGenCutm2_ptTTCut->Clone("muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff");
          muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff->Divide(gpMuonGenEtaMuons_ptGenCutm2_ptTTCut);

          //muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff->GetYaxis()->SetRangeUser(0.95, 1.01);


          if(TEfficiency::CheckConsistency(*muCandGenEtaMuons_ptGenCutm2_ptTTCut, *gpMuonGenEtaMuons_ptGenCutm2_ptTTCut) )
          {
            TEfficiency* ptGenCutm2_ptTTCut_eff = new TEfficiency(*muCandGenEtaMuons_ptGenCutm2_ptTTCut, *gpMuonGenEtaMuons_ptGenCutm2_ptTTCut);
            std::string title = muCandGenEtaMuons_ptGenCutm2_ptTTCut->GetTitle();
            title = std::regex_replace(title, std::regex("\\muCandGenEtaMuons"), "tagging efficiency");
            ptGenCutm2_ptTTCut_eff->SetTitle(title.c_str());
            ptGenCutm2_ptTTCut_eff->SetStatisticOption(TEfficiency::EStatOption::kBUniform );
            ptGenCutm2_ptTTCut_eff->SetPosteriorMode();
            ptGenCutm2_ptTTCut_eff->SetLineColor(15);
            ptGenCutm2_ptTTCut_eff->Draw("APZ");

            muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff->SetLineColor(kRed);
            muCandGenEtaMuons_ptGenCutm2_ptTTCut_Eff->Draw("histsame");
            canvas1Eff->Update();
            ptGenCutm2_ptTTCut_eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.9, 1.01);
            ptGenCutm2_ptTTCut_eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(-2.4, 2.4);
            //ptGenCutm2_ptTTCut_eff->GetGetYaxis()->SetRangeUser(0.9, 1.01);
          }
        }
      //------------------------
        //nominator
/*        TH1D* muCandGenEtaMuons = (TH1D*)subdir->Get("muCandGenEtaMuons");
        TH1D* muCandGenPhiMuons = (TH1D*)subdir->Get("muCandGenPhiMuons");

        TH1D* lostTtMuonPt  = (TH1D*)subdir->Get("lostTtMuonPt");
        TH1D* lostTtMuonEta = (TH1D*)subdir->Get("lostTtMuonEta");
        TH1D* lostTtMuonPhi = (TH1D*)subdir->Get("lostTtMuonPhi");

        canvas1Eff->cd(5);
        canvas1Eff->cd(5)->SetGridx();
        canvas1Eff->cd(5)->SetGridy();


        ttMuonEtaEff->Draw("hist");
        ttMuonEtaEff->GetYaxis()->SetRangeUser(0.8, 1.0);

        TH1D* muCandGenEtaMuonsEff = (TH1D*)muCandGenEtaMuons->Clone("muCandGenEtaMuonsEff");
        muCandGenEtaMuonsEff->Divide(gpMuonEta);
        muCandGenEtaMuonsEff->SetLineColor(kRed);
        muCandGenEtaMuonsEff->Draw("histsame");*/

        //------------------------------------------------------------------------
/*        canvas1Eff->cd(6);
        canvas1Eff->cd(6)->SetGridx();
        canvas1Eff->cd(6)->SetGridy();
        TH2I* betaGenBetaL1Mu = (TH2I*)subdir->Get("betaGenBetaL1Mu");
        //betaGenBetaL1Mu->Draw("colz");
        int bin0 = betaGenBetaL1Mu->GetYaxis()->FindBin(0.);
        TH1D* l1MuVsBetaGen = betaGenBetaL1Mu->ProjectionX("l1MuVsBetaGen", bin0, -1);
        TH1D* allVsBetaGen = betaGenBetaL1Mu->ProjectionX("allVsBetaGen", -1, -1);
        l1MuVsBetaGen->Divide(allVsBetaGen);
        l1MuVsBetaGen->SetTitle("efficiency vs beta gen");
        l1MuVsBetaGen->Draw("hist");*/
      }

      omtfTTAnalyzerDir->cd();
    }
  }


//*************************************************** rate ********************************************

  TH2I* etaGenPtGenBx1 = (TH2I*)omtfTTAnalyzerDir->Get("etaGenPtGenBx1");
  TH1D* ptGenBx1 = etaGenPtGenBx1->ProjectionY("ptGenBx1");

  TH2I* etaGenPtGenBx2 = (TH2I*)omtfTTAnalyzerDir->Get("etaGenPtGenBx2");
  TH1D* ptGenBx2 = etaGenPtGenBx2->ProjectionY("ptGenBx2");

  next.Reset();
  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();

      if(dirName.find("RateAnalyzer_") != string::npos) {
        TDirectory *subdir = gDirectory;
        cout<<"making plots for "<<dirName<<endl;
        canvasName = string("canvas_") + nameLegend  + "_" + dirName;
        TCanvas* canvasRate1 = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 1200);
        canvasRate1->Divide(2, 2);

        canvasRate1->cd(1);

        canvasRate1->cd(1)->SetLogy();
        canvasRate1->cd(1)->SetGridx();
        canvasRate1->cd(1)->SetGridy();

        muCandPt = (TH1D*)subdir->Get("muCandPt");
        muCandPt->SetLineColor(kBlack);
        muCandPt->Draw("hist");
        muCandPt->GetXaxis()->SetRangeUser(0, 100);
        muCandPt->GetYaxis()->SetRangeUser(0.1, 10000);

        TH1D* muCandPtMuons = (TH1D*)subdir->Get("muCandPtMuons");
        muCandPtMuons->SetLineColor(kBlue);
        muCandPtMuons->Draw("histsame");

        TH1D* muCandPtFakes = (TH1D*)subdir->Get("muCandPtFakes");
        muCandPtFakes->SetLineColor(kMagenta);
        //muCandPtFakes->SetFillColor(kMagenta);
        muCandPtFakes->Draw("histsame");

        TH1D* muCandPtWrongTag = (TH1D*)subdir->Get("muCandPtWrongTag");
        muCandPtWrongTag->SetLineColor(kGreen);
        muCandPtWrongTag->Draw("histsame");


        ///////
        canvasRate1->cd(2);
        canvasRate1->cd(2)->SetLogy();
        canvasRate1->cd(2)->SetGridx();
        canvasRate1->cd(2)->SetGridy();

        cout<<"lhcFillingRatio "<<lhcFillingRatio<<endl;
        double scale = 1./eventsCnt * 40000000 * lhcFillingRatio;
        TH1* muCandPt_rateCumul = muCandPt->GetCumulative(false, "_rate");
        muCandPt_rateCumul->Scale(scale);
        muCandPt_rateCumul->GetXaxis()->SetTitle("ttTrack pt tresh [GeV]");
        muCandPt_rateCumul->GetYaxis()->SetTitle("Rate [Hz]");
        muCandPt_rateCumul->SetLineColor(kBlack);
        muCandPt_rateCumul->Draw("hist");
        muCandPt_rateCumul->GetXaxis()->SetRangeUser(0, 100);
        muCandPt_rateCumul->GetYaxis()->SetRangeUser(100, 50000000);

        cout<<" rate at 20GeV "<<muCandPt_rateCumul->GetBinContent(muCandPt_rateCumul->FindBin(20))<<endl; //TODO <<<<<<<<<<<<<

        TH1* muCandPtFakes_rateCumul = muCandPtFakes->GetCumulative(false, "_rate");
        muCandPtFakes_rateCumul->Scale(scale);
        muCandPtFakes_rateCumul->SetLineColor(kMagenta);
        //muCandPtFakes_rateCumul->SetFillColor(kMagenta);
        muCandPtFakes_rateCumul->Draw("histsame");

        TH1* muCandPtWrongTag_rateCumul = muCandPtWrongTag->GetCumulative(false, "_rate");
        muCandPtWrongTag_rateCumul->Scale(scale);
        muCandPtWrongTag_rateCumul->SetLineColor(kGreen);
        muCandPtWrongTag_rateCumul->Draw("histsame");


        TH1* muCandPtMuons_rateCumul = muCandPtMuons->GetCumulative(false, "_rate");
        muCandPtMuons_rateCumul->Scale(scale);
        muCandPtMuons_rateCumul->SetLineColor(kBlue);
        muCandPtMuons_rateCumul->Draw("histsame");
        //canvas1->cd(6)->SetLogy();

        ////////////////////////////////
        canvasRate1->cd(3);
        canvasRate1->cd(3)->SetLogy();
        canvasRate1->cd(3)->SetGridx();
        canvasRate1->cd(3)->SetGridy();

        TH2I* etaGenPtGenLostBx1 = (TH2I*)subdir->Get("etaGenPtGenLostBx1");
        TH1D* ptGenLostBx1 = etaGenPtGenLostBx1->ProjectionY("ptGenLostBx1");
        ptGenLostBx1->SetLineColor(kRed);
        ptGenBx1->Draw("hist");
        //ptGenBx1->GetXaxis()->Set
        ptGenLostBx1->Draw("histsame");

        int ptGenBx1Pt20 = ptGenBx1->Integral(ptGenBx1->FindBin(20), -1);
        int ptGenLostBx1Pt20 = ptGenLostBx1->Integral(ptGenBx1->FindBin(20), -1);
        cout<<" ptGenBx1Pt20 "<<ptGenBx1Pt20<<" ptGenLostBx1Pt20 "<<ptGenLostBx1Pt20<<" eff lost "<<ptGenLostBx1Pt20/(double)ptGenBx1Pt20<<endl;;

        ////////////////////////////////
        canvasRate1->cd(4);
        canvasRate1->cd(4)->SetLogy();
        canvasRate1->cd(4)->SetGridx();
        canvasRate1->cd(4)->SetGridy();

        TH2I* etaGenPtGenLostBx2 = (TH2I*)subdir->Get("etaGenPtGenLostBx2");
        TH1D* ptGenLostBx2 = etaGenPtGenLostBx2->ProjectionY("ptGenLostBx2");
        ptGenLostBx2->SetLineColor(kRed);
        ptGenBx2->Draw("hist");
        //ptGenBx2->GetXaxis()->Set
        ptGenLostBx2->Draw("histsame");

        int ptGenBx2Pt20 = ptGenBx2->Integral(ptGenBx2->FindBin(20), -1);
        int ptGenLostBx2Pt20 = ptGenLostBx2->Integral(ptGenBx2->FindBin(20), -1);
        cout<<" ptGenBx2Pt20 "<<ptGenBx2Pt20<<" ptGenLostBx2Pt20 "<<ptGenLostBx2Pt20<<" eff lost "<<ptGenLostBx2Pt20/(double)ptGenBx2Pt20<<endl;;

      }
    }
  }

  compareFirst = false;
}




void makeCandidatesMatchingPlots(TDirectory* omtfTTAnalyzerDir, const char* nameLegend, double eventsCnt) {
  bool cumulative = false;

  vector< std::pair<std::string, Color_t> > categoryNames;
  //categoryNames.emplace_back("all candidates", kBlack); //0

  categoryNames.emplace_back("muons", kBlue ); //1
  categoryNames.emplace_back("veryLooseMuons", kBlue -2); //2
  //categoryNames.emplace_back("pions", kGreen) ; //3
  categoryNames.emplace_back("pionsDecayedToMu", kGreen); //4
  categoryNames.emplace_back("kaonsDecayedToMu", kCyan + 1); //7

  categoryNames.emplace_back("pionsNotDecayedToMu", kGreen +1); //5
  //categoryNames.emplace_back("kaons", kCyan); //6
  categoryNames.emplace_back("kaonsNotDecayedToMu",kCyan + 2); //8
  categoryNames.emplace_back("otherParts", kGray ); //9
  categoryNames.emplace_back("fakes", kMagenta); //10


  TIter next(omtfTTAnalyzerDir->GetListOfKeys());
  TKey* key;


  TDirectory* subdirAllTTTrack;
  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();
      if(dirName.find("MuCandsMatchingAnalyzer_AllTTTRacks10") != string::npos) {
        subdirAllTTTrack = gDirectory;
      }
    }
  }

  next.Reset();

  double scale = 1./eventsCnt * 40000000 * lhcFillingRatio;

  while ((key = (TKey*)next())) {
    if (key->IsFolder()) {
      omtfTTAnalyzerDir->cd(key->GetName());
      string dirName = key->GetName();
      if(dirName.find("MuCandsMatchingAnalyzer_") != string::npos) {
        TDirectory* subdir = gDirectory;
        cout<<"making purity  plots for "<<dirName<<endl;

        string canvasName = string("canvas__PurityPlots_") + nameLegend  + "_" + dirName;
        TCanvas* canvasPurityPlots = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 800);
        canvasPurityPlots->cd();
        canvasPurityPlots->Divide(2, 2);

        canvasPurityPlots->cd(1);
        canvasPurityPlots->cd(1)->SetGridx();
        canvasPurityPlots->cd(1)->SetGridy();
        canvasPurityPlots->cd(1)->SetLogy();

        TLegend* legend = new TLegend(0.65, 0.55, 0.9, 0.9);

        canvasPurityPlots->cd(1);
        TH1D* allCandPt = (TH1D*)subdir->Get("candPt");
        allCandPt->SetLineColor(kBlack);
        canvasPurityPlots->cd(1);

        TH1* allCandPtCumul = allCandPt->GetCumulative(false, "_cumul");

        if(!cumulative) {
          allCandPt->Draw("C");
          allCandPt->GetXaxis()->SetRangeUser(0, 100);
          allCandPt->GetYaxis()->SetRangeUser(0.5, 100000000);
          allCandPt->GetXaxis()->SetTitle("ttTrack p_{T} [GeV]");
        }
        else {
          allCandPtCumul->Scale(scale);
          allCandPtCumul->Draw("C");
          allCandPtCumul->GetXaxis()->SetRangeUser(0, 100);
          allCandPtCumul->GetYaxis()->SetRangeUser(100, 100000000);
          allCandPtCumul->GetXaxis()->SetTitle("ttTrack p_{T} threshold [GeV]");
          allCandPtCumul->GetYaxis()->SetTitle("rate [Hz]");
        }

        legend->AddEntry(allCandPt, "all candidates");

        ////////////////////////////////////////
        {
          canvasPurityPlots->cd(2);
          canvasPurityPlots->cd(2)->SetLogy();
          canvasPurityPlots->cd(2)->SetGridx();
          canvasPurityPlots->cd(2)->SetGridy();
          TH1D* allCandPt1 = (TH1D*)allCandPt->Clone( (allCandPt->GetName() + std::string("_eff1") ).c_str());

          TH1D* candPtAllTTTracks = (TH1D*)subdirAllTTTrack->Get("candPt");
          candPtAllTTTracks = (TH1D*)candPtAllTTTracks->Clone((candPtAllTTTracks->GetName() + std::string("_denom") ).c_str());
          /*
          allCandPt1->SetTitle("probability to be tagged as muon");
          allCandPt1->Divide(candPtAllTTTracks);
          allCandPt1->SetMarkerStyle(22);
          allCandPt1->SetMarkerSize(0.5);
          allCandPt1->Draw("P hist");
          allCandPt1->GetYaxis()->SetRangeUser(0.0001, 1.1);
          allCandPt1->GetYaxis()->SetTitle("efficiency");

          */

          if(!cumulative) {
            TH1* allCandPtEff = (TH1D*)allCandPt->Clone( (allCandPtCumul->GetName() + std::string("_eff")).c_str() );
            allCandPtEff->Divide(candPtAllTTTracks);
            allCandPtEff->Draw("hist");
            allCandPtEff->GetXaxis()->SetTitle("ttTrack p_{T} [GeV]");
            allCandPtEff->GetYaxis()->SetRangeUser(0.00001, 1.1);
            allCandPtEff->GetYaxis()->SetTitle("efficiency");
          }
          else {
            TH1* candPtAllTTTracksCumul = candPtAllTTTracks->GetCumulative(false, "_cumul");
            TH1* candPtEffCumul = (TH1D*)allCandPtCumul->Clone( (allCandPtCumul->GetName() + std::string("_eff")).c_str() );
            candPtEffCumul->Divide(candPtAllTTTracksCumul);
            candPtEffCumul->Draw("hist");
            candPtEffCumul->GetXaxis()->SetTitle("ttTrack p_{T} threshold [GeV]");
            candPtEffCumul->GetYaxis()->SetRangeUser(0.00001, 1.1);
            candPtEffCumul->GetYaxis()->SetTitle("efficiency");
          }
        }

        THStack * hsPurity = new THStack("hsPurity"," purity - probability that the muon candidate is matched to a given category of tracking particle");
        THStack * hsPurityCumul = new THStack("hsPurity"," hsPurityCumul");

        for(auto& categoryName : categoryNames) {
          std::string histName = "candPt_" + categoryName.first;

          TH1D* candPt = (TH1D*)subdir->Get(histName.c_str());
          if(candPt) {
            candPt->SetLineColor(categoryName.second);
            candPt->SetMarkerColor(categoryName.second);
            /*if(categoryName.first == "veryLooseMuons")
              candPt->SetFillColor(categoryName.second);*/

            legend->AddEntry(candPt, categoryName.first.c_str());

            canvasPurityPlots->cd(1);//events count vs pt for every category of candidates

            if(!cumulative) {
              candPt->Draw("Csame");
            }
            else {
              TH1* candPtCumul = candPt->GetCumulative(false, "_cumul");
              candPtCumul->Scale(scale);
              candPtCumul->Draw("Csame");
            }

            ////////////////////////////////////////
            //efficiency, i.e. probability to tag as a muon the ttTrack of a given matchin category
            canvasPurityPlots->cd(2);

            TH1D* candPtEff = (TH1D*)candPt->Clone((candPt->GetName() + std::string("_eff") ).c_str());

            TH1D* candPtAllTTTracks = (TH1D*)subdirAllTTTrack->Get(histName.c_str());
            candPtAllTTTracks = (TH1D*)candPtAllTTTracks->Clone((candPtAllTTTracks->GetName() + std::string("_denom") ).c_str());
            if(!cumulative) {
              candPtEff->Divide(candPtAllTTTracks);
              candPtEff->SetMarkerStyle(22);
              candPtEff->SetMarkerSize(0.7);
              candPtEff->Draw("P hist same");
            }
            else {
              //efficiency cumulative
              TH1* candPtEffCumul = candPt->GetCumulative(false, "_eff_cumul");
              TH1* candPtAllTTTracksCumul = candPtAllTTTracks->GetCumulative(false, "_cumul");
              candPtEffCumul->Divide(candPtAllTTTracksCumul);
              candPtEffCumul->Draw("histsame");
            }

            /////////////////////////////////////////
            TH1D* candPtPurity = (TH1D*)candPt->Clone( (candPt->GetName() + std::string("_purity") ).c_str() );
            candPtPurity->Divide(allCandPt);
            candPtPurity->SetFillColor(categoryName.second);
            hsPurity->Add(candPtPurity);

            ///////////////////////////////////////////
            TH1* candPtCumul = candPt->GetCumulative(false, "_cumul");
            candPtCumul->Divide(allCandPtCumul);
            candPtCumul->SetFillColor(categoryName.second);
            hsPurityCumul->Add(candPtCumul);


          }
          else {
            cout<<"no hist "<<histName<<endl;
          }
        }

        canvasPurityPlots->cd(1);
        legend->Draw();


        canvasPurityPlots->cd(3);
        hsPurity->Draw("hist");
        hsPurity->GetXaxis()->SetTitle("ttTrack pt [GeV]");
        hsPurity->GetXaxis()->SetRangeUser(0, 100);

        canvasPurityPlots->cd(4);
        hsPurityCumul->Draw("hist");
        hsPurityCumul->GetXaxis()->SetRangeUser(0, 100);
        hsPurityCumul->GetXaxis()->SetTitle("ttTrack pt threshold [GeV]");
      }
    }
  }
}

