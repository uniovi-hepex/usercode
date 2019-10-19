
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


#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <regex>

#include "PlottingTemplate-master/PlotTemplate.C"

using namespace std;

string outPlotsdir = "omf_correlator_plots/";

struct PlotElements {
  TCanvas* canvas = nullptr;
  TEfficiency* omtfEff =  nullptr;
  TEfficiency* correlEff =  nullptr;
  TLegend* legend = nullptr;
};


PlotElements makeEfficiencyPlot(string canvasName, TFile* omtf_tdr_effs, string omtfEffName, string corrFileName, string corrEffName, double xRangeFrom, double xRangeTo, string header) {
  PlotElements plotElements;

  plotElements.omtfEff = (TEfficiency*)omtf_tdr_effs->Get(omtfEffName.c_str());

  plotElements.canvas = CreateCanvas(canvasName.c_str(), false, true);
  plotElements.canvas->cd();

  if(canvasName.find("effVsEta") != string::npos || canvasName.find("effVsAbsEta") != string::npos)
    plotElements.omtfEff->SetTitle(";generated #eta;Efficiency");
  else
    plotElements.omtfEff->SetTitle(";generated p_{T};Efficiency");

  plotElements.omtfEff->Draw("APZ");
  plotElements.canvas->Update();
  plotElements.omtfEff->GetPaintedGraph()->GetXaxis()->SetRangeUser(xRangeFrom, xRangeTo);
  plotElements.omtfEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  plotElements.canvas->Update();


  TFile* corrEffVsPt = new TFile(corrFileName.c_str());

  plotElements.correlEff = (TEfficiency*)corrEffVsPt->Get(corrEffName.c_str());
  plotElements.correlEff->Draw("same PZ");

  DrawCmsSimulationLabel(plotElements.canvas);
  DrawPuLabel(plotElements.canvas, "200 PU");

  TLegend* leg = new TLegend(0.33, 0.14,0.77,0.28);
  leg->SetHeader(header.c_str());
  //leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.03);
  //leg->SetHeader("here is a beautiful header");x
  leg->AddEntry(plotElements.correlEff , "muon correlator", "le");
  leg->AddEntry(plotElements.omtfEff, "OMTF", "lep");
  leg->Draw("same");
  plotElements.canvas->Update();

  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".pdf").c_str());
  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".png").c_str());
  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".root").c_str());

  return plotElements;
}

void MakePlots()
{

  gStyle->SetOptStat(0);

  TFile* omtf_tdr_effs = new TFile("fromCarlos/tdr_effs.root");

  TFile* omtf_tdr_effs_abseta = new TFile("fromCarlos/tdr_effs_abseta.root");


  makeEfficiencyPlot( string("canvas_effVsAbsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000_fullCorr"), omtf_tdr_effs_abseta, "hEta20_q12_cut20PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone.root",
                                                                                            "gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone", 0.8, 1.2,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");



  makeEfficiencyPlot( string("canvas_effVsAbsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15"), omtf_tdr_effs_abseta, "hEta7_15_q12_cut5PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone.root",
                                                                                       "gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone", 0.8, 1.2,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");

 // return;


  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000_fullCorr"), omtf_tdr_effs, "hEta20_q12_cut20PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_clone.root",
                                                                                            "gpMuonGenEtaMuons_withPtCuts_clone", -2.4, 2.4,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");


  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15_fullCorr"), omtf_tdr_effs, "hEta7_15_q12_cut5PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_clone.root",
                                                                                       "gpMuonGenEtaMuons_withPtCuts_clone", -2.4, 2.4,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");



  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000"), omtf_tdr_effs, "hEta20_q12_cut20PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_overlap_clone.root",
                                                                                            "gpMuonGenEtaMuons_withPtCuts_overlap_clone", -1.5, 1.5,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");


  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15"), omtf_tdr_effs, "hEta7_15_q12_cut5PU200",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_overlap_clone.root",
                                                                                       "gpMuonGenEtaMuons_withPtCuts_overlap_clone", -1.5, 1.5,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");



  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut3"), omtf_tdr_effs, "hPt5_SingleMu_PU200_q12_1",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_3.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone", 0, 100,
                                                                                       "single #mu gun, L1 p_{T} #geq 3 GeV");

  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut5"), omtf_tdr_effs, "hPt5_SingleMu_PU200_q12_1",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_5.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone", 0, 100,
                                                                                       "single #mu gun, L1 p_{T} #geq 5 GeV");

  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut10"), omtf_tdr_effs, "hPt10_SingleMu_PU200_q12_1",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_10.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone", 0, 100,
                                                                                       "single #mu gun, L1 p_{T} #geq 10 GeV");

  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut20"), omtf_tdr_effs, "hPt20_SingleMu_PU200_q12_1",
      "plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_20.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone", 0, 100,
                                                                                       "single #mu gun, L1 p_{T} #geq 20 GeV");

  /*
  {
    TEfficiency* hEta7_15_q12_cut5PU200 = (TEfficiency*)tdr_effs->Get("hEta7_15_q12_cut5PU200");

    TCanvas* canvas = CreateCanvas((string("cancas_") +  hEta7_15_q12_cut5PU200->GetName() ).c_str(), false, true);
    canvas->cd();

    hEta7_15_q12_cut5PU200->Draw("APZ");

    TFile* corrEffVsPt = new TFile("plots_MuFlatPt_PU200_t12/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_20.root");

    corrEffVsPt->Get("ptGenPtMuCandMuonsEv0OverlapDenom_clone")->Draw("same PZ");

    DrawCmsSimulationLabel(canvas, "200");
    DrawPuLabel(canvas);

    //SaveCanvas(canvas, plotsDir, name);
  }
   */



  /////////////////////////rate plots//////////////////////////////

/*  {
    TFile* tdr_pt_rate_q12_PU_log200 = new TFile("fromCarlos/tdr_pt_rate_q12_PU_log200.root");
    TCanvas* canvas_tdr_pt_rate_q12_PU_log200 = (TCanvas*)tdr_pt_rate_q12_PU_log200->Get("name");
    canvas_tdr_pt_rate_q12_PU_log200->cd();
    canvas_tdr_pt_rate_q12_PU_log200->ls();
// "Graph"
    canvas_tdr_pt_rate_q12_PU_log200->Draw();

    TFile* corrRate= new TFile("plots_SingleNeutrino_PU200_t13/RateAnalyzer_SingleMuAlgoOverlap20_muCandPt_rate.root");

    TH1* muCandPt_rate = (TH1*)corrRate->Get("muCandPt_rate");
    muCandPt_rate->Scale(1./1000.);
    muCandPt_rate->Draw("same L");
  }*/

  {
    double xRangeFrom = 0;
    double xRangeTo = 50;

    TFile* tdr_rates = new TFile("fromCarlos/tdr_rates.root");
    TGraphAsymmErrors* omtfRate = (TGraphAsymmErrors*)tdr_rates->Get("tdr_pt_rate_PU200__q12");

    TCanvas* canvas = CreateCanvas("rates_omtfRegion", false, true);
    canvas->cd();



    canvas->SetLogy();
    omtfRate->SetTitle(";p_{T} threshold [GeV]; Rate [kHz]");

    omtfRate->SetLineColor(kBlack);
    omtfRate->Draw("APZ");
    canvas->Update();
    omtfRate->GetXaxis()->SetRangeUser(xRangeFrom, xRangeTo);
    omtfRate->GetYaxis()->SetRangeUser(0.1, 5000);
    canvas->Update();


    TFile* corrRate= new TFile("plots_SingleNeutrino_PU200_t13/RateAnalyzer_SingleMuAlgoOverlap20_muCandPt_rate.root");

    TH1* muCandPt_rate = (TH1*)corrRate->Get("muCandPt_rate");
    muCandPt_rate->Scale(1./1000.);

    muCandPt_rate->SetLineColor(kRed);
    muCandPt_rate->Draw("same L");

    TLegend* leg = new TLegend(0.55, 0.65,0.85,0.75);
    //leg->SetBorderSize(0);
    leg->SetFillStyle(1001);
    leg->SetTextSize(0.035);
    //leg->SetHeader("here is a beautiful header");x
    leg->AddEntry(muCandPt_rate , "muon correlator", "le");
    leg->AddEntry(omtfRate, "OMTF", "lep");
    leg->Draw("same");

    DrawCmsSimulationLabel(canvas);
    DrawPuLabel(canvas, "200 PU");
    canvas->Update();

    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".root").c_str());
  }

  {
    TGraphAsymmErrors* muCorrRateVsPu = new TGraphAsymmErrors(4);


    muCorrRateVsPu->SetPoint(0, 140, 1453.33/1000.);  muCorrRateVsPu->SetPointError(0, 0, 0, 303.04/1000., 303.04/1000.);

    muCorrRateVsPu->SetPoint(1, 200, 2245.29/1000.);  muCorrRateVsPu->SetPointError(1, 0, 0, 622.732/1000., 622.732/1000.);

    muCorrRateVsPu->SetPoint(2, 250, 3171.04/1000.);  muCorrRateVsPu->SetPointError(2, 0, 0, 444.035/1000., 444.035/1000.);

    muCorrRateVsPu->SetPoint(4, -100, -100);  muCorrRateVsPu->SetPointError(4, 0.1, 0.1, 0.1, 0.1);

    TCanvas* canvas = CreateCanvas("ratesVsPu_omtfRegion", false, true);

    canvas->cd();


    muCorrRateVsPu->SetLineColor(kRed);
    muCorrRateVsPu->SetMarkerStyle(23);
    muCorrRateVsPu->SetMarkerColor(kRed);
    muCorrRateVsPu->SetTitle("; average pile-up; rete [kHz]");
    muCorrRateVsPu->Draw("APZ");

    canvas->Update();
    muCorrRateVsPu->GetXaxis()->SetRangeUser(0, 300);
    muCorrRateVsPu->GetYaxis()->SetRangeUser(0, 25);

    muCorrRateVsPu->RemovePoint(4);

    DrawCmsSimulationLabel(canvas);
    DrawPuLabel(canvas, "14 TeV");
    canvas->Update();
    canvas->cd();

    canvas->SaveAs( (outPlotsdir + "rateVsPu_omtfRegion" + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + "rateVsPu_omtfRegion" + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + "rateVsPu_omtfRegion" + ".root").c_str());
  }

  {
    string canvasName = "HSCP_efficiecny_vs_beta";
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), false, true);
    canvas->cd();

    TLegend* leg = new TLegend(0.13, 0.22,0.46,0.36);
    leg->SetHeader("#tilde{#tau}, m = 494 and 1599 GeV");
    leg->SetFillStyle(1001);
    leg->SetTextSize(0.03);


    TFile* totalEffFile = new TFile("plots_HSCP/EfficiencyAnalyser_HscpAlgoSoftCuts20_ptGenFrom_20_ptGenTo_100000_allVsBetaGenSum_clone_totalEff_.root");
    TEfficiency* totalEff = (TEfficiency*)totalEffFile->Get("allVsBetaGenSum_clone");
    totalEff->SetTitle(";generated #beta; efficiency");
    totalEff->SetLineColor(kBlack);
    totalEff->SetMarkerColor(kBlack);
    totalEff->SetMarkerStyle(24);
    totalEff->Draw("APZ");
    canvas->Update();
    totalEff->GetPaintedGraph()->GetXaxis()->SetRangeUser(0, 1);
    totalEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
    canvas->Update();
    leg->AddEntry(totalEff , "total efficiency", "lep");

    DrawCmsSimulationLabel(canvas);
    DrawPuLabel(canvas, "no PU");


    TFile* hscpAlgoEffFile = new TFile("plots_HSCP/EfficiencyAnalyser_HscpAlgoSoftCuts20_ptGenFrom_20_ptGenTo_100000_allVsBetaGen_clone_algoEff_.root");
    TEfficiency* hscpAlgoEff = (TEfficiency*)hscpAlgoEffFile->Get("allVsBetaGen_clone");
    hscpAlgoEff->SetLineColor(kBlue);
    hscpAlgoEff->SetMarkerStyle(22);
    hscpAlgoEff->SetMarkerColor(kBlue);
    hscpAlgoEff->Draw("same PZ");
    canvas->Update();
    leg->AddEntry(hscpAlgoEff, "HSCP algo", "lep");


    TFile* sinleMuAlgoEffFile = new TFile("plots_HSCP/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_20_ptGenTo_100000_allVsBetaGen_clone_algoEff_.root");
    TEfficiency* sinleMuAlgoEff = (TEfficiency*)sinleMuAlgoEffFile->Get("allVsBetaGen_clone");
    sinleMuAlgoEff->SetLineColor(kRed);
    sinleMuAlgoEff->SetMarkerStyle(23);
    sinleMuAlgoEff->SetMarkerColor(kRed);
    sinleMuAlgoEff->Draw("same PZ");
    canvas->Update();
    leg->AddEntry(sinleMuAlgoEff, "#mu algo", "lep");


    leg->Draw("same");
    canvas->Update();

    canvas->SaveAs( (outPlotsdir + canvasName + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + canvasName + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + canvasName + ".root").c_str());
  }



  {
    double xRangeFrom = 0;
    double xRangeTo = 50;

    string canvaseNAme = "ratesHSCP_omtfRegion";
    TCanvas* canvas = CreateCanvas(canvaseNAme.c_str(), false, true);
    canvas->cd();


    canvas->SetLogy();

    TFile* corrRate= new TFile("plots_SingleNeutrino_PU200_t13/RateAnalyzer_HscpAlgoSoftCuts20_muCandPt_rate.root");

    TH1* muCandPt_rate = (TH1*)corrRate->Get("muCandPt_rate");
    muCandPt_rate->SetTitle(";p_{T} threshold [GeV]; Rate [kHz]");
    muCandPt_rate->Scale(1./1000.);

    muCandPt_rate->SetLineColor(kRed);

    muCandPt_rate->GetXaxis()->SetRangeUser(0, 50);
    muCandPt_rate->GetYaxis()->SetRangeUser(0.1, 5000);
    muCandPt_rate->Draw("L"); //APZ

    TLegend* leg = new TLegend(0.5, 0.35, 0.87, 0.4);
    //leg->SetBorderSize(0);
    leg->SetFillStyle(1001);
    leg->SetTextSize(0.035);
    leg->SetHeader("HSCP algo rate, |#eta| < 2.4");
    //leg->AddEntry(muCandPt_rate , "muon correlator", "le");
    leg->Draw("same");

    DrawCmsSimulationLabel(canvas);
    DrawPuLabel(canvas, "200 PU");
    canvas->Update();

    canvas->SaveAs( (outPlotsdir + canvaseNAme + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + canvaseNAme + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + canvaseNAme + ".root").c_str());
  }


  return;

}
