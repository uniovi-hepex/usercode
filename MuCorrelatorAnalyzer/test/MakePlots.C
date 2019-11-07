
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

string outPlotsdir = "omtf_correlator_plots_L1TkMuonsTP/";

struct PlotElements {
  TCanvas* canvas = nullptr;
  TEfficiency* omtfEff =  nullptr;
  TEfficiency* correlEff =  nullptr;
  TEfficiency* l1TkMuEff =  nullptr;

  TLegend* legend = nullptr;

  void update() {
    //canvas->Update();
    canvas->Modified();

    canvas->SaveAs( (outPlotsdir + canvas->GetName() + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + canvas->GetName()  + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + canvas->GetName()  + ".root").c_str());
  }
};

int lineWidth = 2;
bool showTamu = true;
bool showBayes = true;

PlotElements makeEfficiencyPlot(string canvasName, TFile* omtf_tdr_effs, string omtfEffName,
    string corrFileName, string corrEffName,
    string bayesCorrDir, string tkMuDir,
    double xRangeFrom, double xRangeTo, string header)
{
  PlotElements plotElements;

  plotElements.omtfEff = (TEfficiency*)omtf_tdr_effs->Get(omtfEffName.c_str());

  plotElements.canvas = CreateCanvas(canvasName.c_str(), false, true);
  plotElements.canvas->cd();

  if(canvasName.find("effVsEta") != string::npos)
    plotElements.omtfEff->SetTitle(";generated #eta;Efficiency");
  else if(canvasName.find("effVsAbsEta") != string::npos)
    plotElements.omtfEff->SetTitle(";generated #left|#eta#right|;Efficiency");
  else
    plotElements.omtfEff->SetTitle(";generated p_{T};Efficiency");

  plotElements.omtfEff->SetLineWidth(lineWidth);
  plotElements.omtfEff->SetMarkerStyle(20);
  plotElements.omtfEff->SetMarkerColor(kBlack);

  plotElements.omtfEff->Draw("APZ");
  plotElements.canvas->Update();
  plotElements.omtfEff->GetPaintedGraph()->GetXaxis()->SetRangeUser(xRangeFrom, xRangeTo);
  plotElements.omtfEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.05);
  plotElements.canvas->Update();




  DrawCmsSimulationLabel(plotElements.canvas);
  DrawPuLabel(plotElements.canvas, "200 PU");

  TLegend* leg = new TLegend(0.33, 0.13,0.77,0.38);
  leg->SetHeader(header.c_str());
  //leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetMargin(0.2);
  //leg->SetHeader("here is a beautiful header");x



  if(showBayes) {
    TFile* corrEffVsPt = new TFile((bayesCorrDir + corrFileName).c_str());

    plotElements.correlEff = (TEfficiency*)corrEffVsPt->Get(corrEffName.c_str());
    plotElements.correlEff->SetLineWidth(lineWidth);
    plotElements.correlEff->SetMarkerStyle(22);
    plotElements.correlEff->SetMarkerColor(kRed);
    plotElements.correlEff->Draw("same PZ");
    leg->AddEntry(plotElements.correlEff , "Track + Stubs", "lep");
  }

  leg->AddEntry(plotElements.omtfEff, "OMTF", "lep");

  if(tkMuDir.size() != 0 && showTamu) {
    TFile* l1TkMuEffVsPt = new TFile((tkMuDir + corrFileName).c_str());

    plotElements.l1TkMuEff = (TEfficiency*)l1TkMuEffVsPt->Get(corrEffName.c_str());
    plotElements.l1TkMuEff->SetLineColor(kGreen +2);
    plotElements.l1TkMuEff->SetLineWidth(lineWidth);
    plotElements.l1TkMuEff->SetMarkerStyle(23);
    plotElements.l1TkMuEff->SetMarkerColor(kGreen +2);
    plotElements.l1TkMuEff->Draw("same PZ");
    leg->AddEntry(plotElements.l1TkMuEff, "Track + OMTF", "lep");
  }

  leg->Draw("same");

/*  plotElements.canvas->Update();

  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".pdf").c_str());
  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".png").c_str());
  plotElements.canvas->SaveAs( (outPlotsdir + canvasName + ".root").c_str());*/

  plotElements.legend = leg;
  plotElements.canvas->Update();

  return plotElements;
}

void MakePlots()
{

  gStyle->SetOptStat(0);

  TFile* omtf_tdr_effs = new TFile("fromCarlos/tdr_effs_pt3_redoAnalyzer_alleta.root");

  TFile* omtf_tdr_effs_abseta = new TFile("fromCarlos/tdr_effs_pt3_redoAnalyzer.root");


  //string l1TkMuonsDir = "plots_MuFlatPt_PU200_v1_t17_bayesOMTFonly_L1TkMuonsTP/";
  string l1TkMuonsDir = "plots_MuFlatPt_PU200_v1_t17_bayesOMTFonly_L1TkMuonsTP/";


  PlotElements plotElements;

  plotElements =
  makeEfficiencyPlot( string("canvas_effVsAbsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000_fullCorr"), omtf_tdr_effs_abseta, "hEta20_q12_cut20PU200",
                       "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone.root",
                                                                                                             "gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone",
                                                                                                             "plots_MuFlatPt_PU200_t12/",
                                                                                                             l1TkMuonsDir,
                                                                                            0.8, 1.2,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");
  plotElements.update();
  plotElements.canvas->Update();

  plotElements =
  makeEfficiencyPlot( string("canvas_effVsAbsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15"), omtf_tdr_effs_abseta, "hEta7_15_q12_cut5PU200",
                       "EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone.root",
                                                                                                        "gpMuonGenEtaMuons_withPtCuts_overlap_abs_clone",
                                                                                                        "plots_MuFlatPt_PU200_t12/",
                                                                                                        l1TkMuonsDir,
                                                                                       0.8, 1.2,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");
  plotElements.update();
 // return;


  plotElements =
  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000_fullCorr"), omtf_tdr_effs, "hEta20_q12_cut20PU200",
                       "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_clone.root",
                                                                                                             "gpMuonGenEtaMuons_withPtCuts_clone",
                                                                                                             "plots_MuFlatPt_PU200_t12/",
                                                                                                             l1TkMuonsDir,
                                                                                            -2.4, 2.4,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");
  plotElements.legend->SetX1NDC(0.33); plotElements.legend->SetY1NDC(0.44);
  plotElements.legend->SetX2NDC(0.77); plotElements.legend->SetY2NDC(0.61);
  plotElements.update();


  plotElements =
  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15_fullCorr"), omtf_tdr_effs, "hEta7_15_q12_cut5PU200",
                       "EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_clone.root",
                                                                                       "gpMuonGenEtaMuons_withPtCuts_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       -2.4, 2.4,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");

  plotElements.legend->SetX1NDC(0.40); plotElements.legend->SetY1NDC(0.44);
  plotElements.legend->SetX2NDC(0.84); plotElements.legend->SetY2NDC(0.61);
  plotElements.update();


  plotElements =
  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut20_ptGenFrom_25_ptGenTo_10000"), omtf_tdr_effs, "hEta20_q12_cut20PU200",
      "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_gpMuonGenEtaMuons_withPtCuts_overlap_clone.root",
                                                                                            "gpMuonGenEtaMuons_withPtCuts_overlap_clone",
                                                                                            "plots_MuFlatPt_PU200_t12/",
                                                                                            l1TkMuonsDir,
                                                                                            -1.5, 1.5,
                                                                                            "p_{T}^{gen} > 25 GeV, L1 p_{T} #geq 20 GeV");
  plotElements.legend->SetX1NDC(0.31); plotElements.legend->SetY1NDC(0.44);
  plotElements.legend->SetX2NDC(0.75); plotElements.legend->SetY2NDC(0.61);
  plotElements.update();

  plotElements =
  makeEfficiencyPlot( string("canvas_effVsEta_L1ptCut5_ptGenFrom_7_ptGenTo_15"), omtf_tdr_effs, "hEta7_15_q12_cut5PU200",
      "EfficiencyAnalyser_SingleMuAlgo5_ptGenFrom_7_ptGenTo_15_gpMuonGenEtaMuons_withPtCuts_overlap_clone.root",
                                                                                       "gpMuonGenEtaMuons_withPtCuts_overlap_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       -1.5, 1.5,
                                                                                       "7 < p_{T}^{gen} < 15 GeV, L1 p_{T} #geq 5 GeV");
  plotElements.legend->SetX1NDC(0.34); plotElements.legend->SetY1NDC(0.44);
  plotElements.legend->SetX2NDC(0.78); plotElements.legend->SetY2NDC(0.61);
  plotElements.update();

  plotElements =
  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut3"), omtf_tdr_effs, "eff_pt3_q12_1GeVPU200",
      "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_3.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       0, 100,
                                                                                       "#splitline{single #mu gun, 0.82 < #left|#eta^{gen}#right| < 1.24}{L1 p_{T} #geq 3 GeV}");
  plotElements.update();
 /*
  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut5"), omtf_tdr_effs, "eff_pt5_q12_1GeVPU200",
      "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_5.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       0, 100,
                                                                                       "#splitline{single #mu gun, 0.82 < #left|#eta^{gen}#right| < 1.24}{L1 p_{T} #geq 5 GeV}");
*/
  plotElements =
  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut10"), omtf_tdr_effs, "eff_pt10_q12_1GeVPU200",
      "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_10.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       0, 100,
                                                                                       "#splitline{single #mu gun, 0.82 < #left|#eta^{gen}#right| < 1.24}{L1 p_{T} #geq 10 GeV}");
  plotElements.update();

  plotElements =
  makeEfficiencyPlot( string("canvas_effVsPt_L1ptCut20"), omtf_tdr_effs, "eff_pt20_q12_1GeVPU200",
      "EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_25_ptGenTo_10000_ptGenPtMuCandMuonsEv0OverlapDenom_clone_ptCut_20.root",
                                                                                       "ptGenPtMuCandMuonsEv0OverlapDenom_clone",
                                                                                       "plots_MuFlatPt_PU200_t12/",
                                                                                       l1TkMuonsDir,
                                                                                       0, 100,
                                                                                       "#splitline{single #mu gun, 0.82 < #left|#eta^{gen}#right| < 1.24}{L1 p_{T} #geq 20 GeV}");
  plotElements.update();
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

  {//rate vs pt
    double xRangeFrom = 0;
    double xRangeTo = 50;

    TFile* tdr_rates = new TFile("fromCarlos/tdr_rates.root");
    TGraphAsymmErrors* omtfRate = (TGraphAsymmErrors*)tdr_rates->Get("tdr_pt_rate_PU200__q12");

    TCanvas* canvas = CreateCanvas("rates_omtfRegion", false, true);
    canvas->cd();


    TLegend* leg = new TLegend(0.51, 0.61,0.81, 0.76);
    leg->SetHeader("0.82 < #left|#eta^{L1}#right| < 1.24");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    //leg->SetHeader("here is a beautiful header");x
    leg->AddEntry(omtfRate, "OMTF", "lep");

    canvas->SetLogy();
    omtfRate->SetTitle(";p_{T} threshold [GeV]; Rate [kHz]");

    omtfRate->SetLineColor(kBlack);
    omtfRate->SetLineWidth(lineWidth);
    omtfRate->SetMarkerStyle(20);
    omtfRate->SetMarkerColor(kBlack);
    omtfRate->Draw("APZ");
    canvas->Update();
    omtfRate->GetXaxis()->SetRangeUser(xRangeFrom, xRangeTo);
    omtfRate->GetYaxis()->SetRangeUser(0.1, 5000);
    canvas->Update();

    if(l1TkMuonsDir.size() != 0 && showTamu) {
      TFile* corrRate= new TFile("plots_SingleNeutrino_PU200_v1_t17_bayesOMTF_L1TkMuonsTP/RateAnalyzer_SingleMuAlgoOverlap20_muCandPt_rate.root");

      TH1* muCandPt_rate = (TH1*)corrRate->Get("muCandPt_rate");
      muCandPt_rate->Scale(1./1000.);

      muCandPt_rate->SetLineColor(kGreen +2);
      muCandPt_rate->SetLineWidth(lineWidth);
      muCandPt_rate->SetMarkerStyle(22);
      muCandPt_rate->SetMarkerColor(kGreen +2);
      muCandPt_rate->Draw("same PZ");

      leg->AddEntry(muCandPt_rate , "Tracks + OMTF", "le");
    }

    if(showBayes) {
      TFile* corrRate= new TFile("plots_SingleNeutrino_PU200_t17/RateAnalyzer_SingleMuAlgoOverlap20_muCandPt_rate.root");

      TH1* muCandPt_rate = (TH1*)corrRate->Get("muCandPt_rate");
      muCandPt_rate->Scale(1./1000.);

      muCandPt_rate->SetLineColor(kRed);
      muCandPt_rate->SetLineWidth(lineWidth);
      muCandPt_rate->SetMarkerStyle(22);
      muCandPt_rate->SetMarkerColor(kRed);
      muCandPt_rate->Draw("same PZ");
      leg->AddEntry(muCandPt_rate , "Tracks + Stubs", "le");
    }

    DrawCmsSimulationLabel(canvas);
    DrawPuLabel(canvas, "200 PU");

    leg->Draw("same");
    canvas->Update();

    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".pdf").c_str());
    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".png").c_str());
    canvas->SaveAs( (outPlotsdir + "rates_omtfRegion" + ".root").c_str());
  }

  {//rates vs pu
    TFile* tdr_rates = new TFile("fromCarlos/tdr_rates_pu.root");
    TGraphAsymmErrors* omtfRateVsPu = (TGraphAsymmErrors*)tdr_rates->Get("mc_rates_pu_pt20__q12");

    TGraphAsymmErrors* muCorrRateVsPu = new TGraphAsymmErrors(3);


    muCorrRateVsPu->SetPoint(0, 140, 1453.33/1000.);  muCorrRateVsPu->SetPointError(0, 0, 0, 303.04/1000., 303.04/1000.);

    muCorrRateVsPu->SetPoint(1, 200, 2114.03/1000.);  muCorrRateVsPu->SetPointError(1, 0, 0, 362.553/1000., 362.553/1000.); //20GeV 2114.03 error 362.553

    muCorrRateVsPu->SetPoint(2, 250, 3171.04/1000.);  muCorrRateVsPu->SetPointError(2, 0, 0, 444.035/1000., 444.035/1000.);

    //muCorrRateVsPu->SetPoint(4, -100, -100);  muCorrRateVsPu->SetPointError(4, 0.1, 0.1, 0.1, 0.1);

    TCanvas* canvas = CreateCanvas("ratesVsPu_omtfRegion", false, true);

    canvas->cd();

    omtfRateVsPu->SetLineColor(kBlack);
    omtfRateVsPu->SetLineWidth(lineWidth);
    omtfRateVsPu->SetMarkerStyle(20);
    omtfRateVsPu->SetMarkerColor(kBlack);
    omtfRateVsPu->Draw("APZ");

    omtfRateVsPu->GetFunction("f")->SetLineColor(kBlack);
    omtfRateVsPu->GetFunction("f")->SetLineWidth(1);
    omtfRateVsPu->GetFunction("f")->SetLineStyle(7);

    canvas->Update();
    omtfRateVsPu->GetXaxis()->SetRangeUser(0, 350);
    omtfRateVsPu->GetYaxis()->SetRangeUser(0, 25);

    muCorrRateVsPu->SetLineColor(kRed);
    muCorrRateVsPu->SetLineWidth(lineWidth);
    muCorrRateVsPu->SetMarkerStyle(23);
    muCorrRateVsPu->SetMarkerColor(kRed);
    muCorrRateVsPu->SetTitle("; average pile-up; Rate [kHz]");

    TF1* fit = new TF1("muCorrRateVsPu_fit", "pol1", 0, 350);
    fit->SetLineStyle(7);
    fit->SetLineWidth(1);
    fit->SetParLimits(0, -0.01, 0.01);
    muCorrRateVsPu->Fit(fit);
    muCorrRateVsPu->Draw("same PZ");
    fit->Draw("same");
    fit->SetLineStyle(7);
    fit->SetLineWidth(1);

    //muCorrRateVsPu->RemovePoint(4);

    TLegend* leg = new TLegend(0.2, 0.6,0.6,0.8);
    leg->SetHeader("#splitline{L1 p_{T} #geq 20 GeV}{0.82< #left|#eta^{gen}#right|<1.24}");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);

    leg->AddEntry(omtfRateVsPu, "OMTF", "lep");
    leg->AddEntry(muCorrRateVsPu , "Tracks + Stubs", "le");
    leg->Draw("same");

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

    TLegend* leg = new TLegend(0.13, 0.45,0.46,0.65);
    //leg->SetHeader("#splitline{#tilde{#tau} m = 494 and 1599 GeV}{#left|#eta^{gen}#right|<2.4}");
    leg->SetHeader("#tilde{#tau} m = 494 and 1599 GeV #left|#eta^{gen}#right|<2.4");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);


    TFile* totalEffFile = new TFile("plots_HSCP/EfficiencyAnalyser_HscpAlgoSoftCuts20_ptGenFrom_20_ptGenTo_100000_allVsBetaGenSum_clone_totalEff_.root");
    TEfficiency* totalEff = (TEfficiency*)totalEffFile->Get("allVsBetaGenSum_clone");
    totalEff->SetTitle(";generated #beta; efficiency");
    totalEff->SetLineColor(kBlack);
    totalEff->SetLineWidth(lineWidth);
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
    hscpAlgoEff->SetLineWidth(lineWidth);
    hscpAlgoEff->SetMarkerStyle(22);
    hscpAlgoEff->SetMarkerColor(kBlue);
    hscpAlgoEff->Draw("same PZ");
    canvas->Update();
    leg->AddEntry(hscpAlgoEff, "HSCP algorithm, L1 p_{T} #geq 20 GeV", "lep");


    TFile* sinleMuAlgoEffFile = new TFile("plots_HSCP/EfficiencyAnalyser_SingleMuAlgo20_ptGenFrom_20_ptGenTo_100000_allVsBetaGen_clone_algoEff_.root");
    TEfficiency* sinleMuAlgoEff = (TEfficiency*)sinleMuAlgoEffFile->Get("allVsBetaGen_clone");
    sinleMuAlgoEff->SetLineColor(kRed);
    sinleMuAlgoEff->SetMarkerStyle(23);
    sinleMuAlgoEff->SetMarkerColor(kRed);
    sinleMuAlgoEff->Draw("same PZ");
    canvas->Update();
    leg->AddEntry(sinleMuAlgoEff, "#mu algorithm, L1 p_{T} #geq 20 GeV", "lep");


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
    muCandPt_rate->SetLineWidth(lineWidth);

    muCandPt_rate->GetXaxis()->SetRangeUser(0, 50);
    muCandPt_rate->GetYaxis()->SetRangeUser(0.1, 5000);
    muCandPt_rate->Draw("L"); //APZ

    TLegend* leg = new TLegend(0.5, 0.35, 0.87, 0.4);
    //leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->SetHeader("#splitline{HSCP algorithm rate}{#left|#eta^{L1}#right| < 2.4}");
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
