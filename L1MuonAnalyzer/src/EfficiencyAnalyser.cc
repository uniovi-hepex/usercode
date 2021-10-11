/*
 * EfficiencyAnalyser.cpp
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/interface/EfficiencyAnalyser.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

namespace L1MuAn {

EfficiencyAnalyser::EfficiencyAnalyser(bool useUpt): useUpt(useUpt) {
  // TODO Auto-generated constructor stub

}

EfficiencyAnalyser::~EfficiencyAnalyser() {
  // TODO Auto-generated destructor stub
}


PtGenVsPtCand::PtGenVsPtCand(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, int nBins, double binsFrom, double binsTo, bool useUpt):
    EfficiencyAnalyser(useUpt), etaFrom(etaFrom), etaTo(etaTo), qualityCut(qualityCut) {
  std::ostringstream histName;

  string ptPrefix = this->useUpt ? "upt" : "pt";

  histName<<name<<"_"<<ptPrefix<<"_ptGenVsPtCand_eta_";
  histName<<etaFrom<<"_"<<etaTo<<"_qualityCut_"<<qualityCut;

  std::ostringstream histTitle;
  histTitle<<name<<" "<<ptPrefix<<" ptGenVsPtCand eta: ";
  histTitle<<etaFrom<<" - "<<etaTo<<" quality >= "<<qualityCut<<";ptGen [GeV]; "<<ptPrefix<<" L1 cand [GeV]";

  ptGenVsPtCand = subDir.make<TH2D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo, 200, 0, 200);

}

void PtGenVsPtCand::fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) {
  if(ptGen >= ptGenVsPtCand->GetXaxis()->GetXmax())
    ptGen = ptGenVsPtCand->GetXaxis()->GetXmax() - 0.01;

  if( (etaFrom <= fabs(etaGen) ) && (fabs(etaGen) <= etaTo) ) {
    if(l1MuonCand.hwQual >= qualityCut) {

      double candPt = useUpt ? l1MuonCand.uptGev : l1MuonCand.ptGev;
      //std::cout<<"useUpt "<<useUpt<<" l1MuonCand.uptGev  "<<l1MuonCand.uptGev<<" l1MuonCand.ptGev "<<l1MuonCand.ptGev<<std::endl;
      if(candPt >= ptGenVsPtCand->GetYaxis()->GetXmax())
        candPt = ptGenVsPtCand->GetYaxis()->GetXmax() - 0.01;

      ptGenVsPtCand->Fill(ptGen, candPt);
    }
    else
      ptGenVsPtCand->Fill(ptGen, 0);
  }
}




EfficiencyVsPhi::EfficiencyVsPhi(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt):
    EfficiencyAnalyser(useUpt), etaFrom(etaFrom), etaTo(etaTo), qualityCut(qualityCut), ptGenCut(ptGenCut), ptL1Cut(ptL1Cut) {
  std::ostringstream histName;

  histName<<name<<"_allCandsPhi_eta_"<<etaFrom<<"_"<<etaTo<<"_qualityCut_"<<qualityCut<<"_ptGenCut_"<<ptGenCut;

  std::ostringstream histTitle;
  histTitle<<name<<" allCands phi, eta: "<<etaFrom<<" - "<<etaTo<<" quality >= "<<qualityCut<<" ptGenCut "<<ptGenCut<<";phi; entries";

  allCands = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, -M_PI, M_PI);

  histName.str("");
  histName<<name<<"_aceptedCandsPhi_eta_"<<etaFrom<<"_"<<etaTo<<"_qualityCut_"<<qualityCut<<"_ptGenCut_"<<ptGenCut<<"_ptL1Cut_"<<ptL1Cut;
  histTitle.str("");
  histTitle<<name<<" allCands Phi, eta: "<<etaFrom<<" - "<<etaTo<<" quality >= "<<qualityCut<<" ptGenCut "<<ptGenCut<<" ptL1Cut "<<ptL1Cut<<";phi; entries";

  aceptedCands = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, -M_PI, M_PI);
}

void EfficiencyVsPhi::fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) {
  if( (etaFrom <= fabs(etaGen) ) && (fabs(etaGen) <= etaTo) && ptGen >= ptGenCut ) {
    allCands->Fill(phiGen);
    double candPt = useUpt ? l1MuonCand.uptGev : l1MuonCand.ptGev;
    if(l1MuonCand.hwQual >= qualityCut && candPt >=  ptL1Cut) {
      aceptedCands->Fill(phiGen);
    }
  }
}




EfficiencyVsEta::EfficiencyVsEta(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt):
  EfficiencyAnalyser(useUpt), qualityCut(qualityCut), ptGenCut(ptGenCut), ptL1Cut(ptL1Cut) {
  std::ostringstream histName;

  histName<<name<<"_allCandsEta_"<<"_qualityCut_"<<qualityCut<<"_ptGenCut_"<<ptGenCut;

  std::ostringstream histTitle;
  histTitle<<name<<" allCands eta : "<<" quality >= "<<qualityCut<<" ptGenCut "<<ptGenCut<<";eta; entries";

  allCands = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, -2.1, 2.1);

  histName.str("");
  histName<<name<<"_aceptedCandsEta_"<<"_qualityCut_"<<qualityCut<<"_ptGenCut_"<<ptGenCut<<"_ptL1Cut_"<<ptL1Cut;
  histTitle.str("");
  histTitle<<name<<" allCands eta: "<<" quality >= "<<qualityCut<<" ptGenCut "<<ptGenCut<<" ptL1Cut "<<ptL1Cut<<";eta; entries";

  aceptedCands = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, -2.1, 2.1);
}

void EfficiencyVsEta::fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) {
  if( ptGen >= ptGenCut ) {
    allCands->Fill(etaGen);
    double candPt = useUpt ? l1MuonCand.uptGev : l1MuonCand.ptGev;
    if(l1MuonCand.hwQual >= qualityCut && candPt >=  ptL1Cut) {
      aceptedCands->Fill(etaGen);
    }
  }
}




EfficiencyPtGenVsDxy::EfficiencyPtGenVsDxy(TFileDirectory& subDir, std::string name, int qualityCut, double ptL1Cut, int nBinsPt, int nBinsDxy, bool ifPtBelowCut, bool useUpt):
    EfficiencyAnalyser(useUpt), qualityCut(qualityCut), ptL1Cut(ptL1Cut), ifPtBelowCut(ifPtBelowCut) {
  std::ostringstream histName;

  string ptPrefix = this->useUpt ? "upt" : "pt";

  string ifPtBelowCutPrefix = this->ifPtBelowCut ? "_ifPtBelowCut" : "";

  histName<<name<<"_allCandsPtGenVsDxy"<<ifPtBelowCutPrefix<<"_qualityCut_"<<qualityCut<<"_"<<ptPrefix<<"L1Cut_"<<ptL1Cut;

  std::ostringstream histTitle;
  histTitle<<name<<" allCandsPtGenVsDxy "<<ifPtBelowCutPrefix;
  histTitle<<" quality >= "<<qualityCut<<" "<<ptPrefix<<" L1Cut "<<ptL1Cut<<" [GeV] "<<";ptGen [GeV]; dxy [cm]";

  allCands = subDir.make<TH2D>(histName.str().c_str(), histTitle.str().c_str(), nBinsPt, 0, 200, 300/20, 0, 300);

  histName.str("");
  histTitle.str("");

  histName<<name<<"_aceptedCandsPtGenVsDxy"<<ifPtBelowCutPrefix<<"_qualityCut_"<<qualityCut<<"_"<<ptPrefix<<"L1Cut_"<<ptL1Cut;
  histTitle<<name<<" aceptedCandsPtGenVsDxy "<<ifPtBelowCutPrefix;
  histTitle<<" quality >= "<<qualityCut<<" "<<ptPrefix<<" L1Cut "<<ptL1Cut<<" [GeV] "<<";ptGen [GeV]; dxy [cm]";

  aceptedCands = subDir.make<TH2D>(histName.str().c_str(), histTitle.str().c_str(), nBinsPt, 0, 200, nBinsDxy, 0, 300);

  histName.str("");
  histTitle.str("");

  histName<<name<<"_effPtGenVsDxy"<<ifPtBelowCutPrefix<<"_qualityCut_"<<qualityCut<<"_"<<ptPrefix<<"L1Cut_"<<ptL1Cut;
  histTitle<<name<<" effPtGenVsDxy "<<ifPtBelowCutPrefix;
  histTitle<<" quality >= "<<qualityCut<<" "<<ptPrefix<<" L1Cut "<<ptL1Cut<<" [GeV] "<<";ptGen [GeV]; dxy [cm]";

  efficiency = subDir.make<TH2D>(histName.str().c_str(), histTitle.str().c_str(), nBinsPt, 0, 200, nBinsDxy, 0, 300);
}

EfficiencyPtGenVsDxy::~EfficiencyPtGenVsDxy() {
  bool res = efficiency->Divide(aceptedCands, allCands);
  LogTrace("l1MuonAnalyzerOmtf") <<"~EfficiencyPtGenVsDxy divide result "<<res;
}

void EfficiencyPtGenVsDxy::fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) {
  if(ptGen >= allCands->GetXaxis()->GetXmax())
    ptGen = allCands->GetXaxis()->GetXmax() - 0.01;

  allCands->Fill(ptGen, dxyGen);

  double candPt = useUpt ? l1MuonCand.uptGev : l1MuonCand.ptGev;
  if(l1MuonCand.hwQual >= qualityCut && candPt >=  ptL1Cut) {
    if(ifPtBelowCut) {
      if(l1MuonCand.ptGev < ptL1Cut)
        aceptedCands->Fill(ptGen, dxyGen);
    }
    else
      aceptedCands->Fill(ptGen, dxyGen);
  }
}

//this function is not called, the division is done in destructor
void EfficiencyPtGenVsDxy::write() {
  //bool res = efficiency->Divide(aceptedCands, allCands);

  allCands->Write();
  aceptedCands->Write();

  //LogTrace("l1MuonAnalyzerOmtf") <<"EfficiencyPtGenVsDxy::write() divide result "<<res;
  efficiency->Write();
}

LikelihoodDistribution::LikelihoodDistribution(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt):
  EfficiencyAnalyser(useUpt), qualityCut(qualityCut), ptGenCut(ptGenCut), ptL1Cut(ptL1Cut) {
  std::ostringstream histName;
  histName<<name<<"_likelihood_"<<"_qualityCut_"<<qualityCut<<"_ptGenCut_"<<ptGenCut<<"_ptL1Cut_"<<ptL1Cut;
  std::ostringstream histTitle;
  histTitle<<name<<" Likelihood Distribution "<<" quality >= "<<qualityCut<<" ptGenCut "<<ptGenCut<<" ptL1Cut "<<ptL1Cut<<";Likelihood; refLayer";

  distribution = subDir.make<TH2I>(histName.str().c_str(), histTitle.str().c_str(), nBins, 0, 1200, 8, -0.5, 7.5); //TODO change the bins if needed
}

void LikelihoodDistribution::fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) {
  if( ptGen >= ptGenCut ) {
    double candPt = useUpt ? l1MuonCand.uptGev : l1MuonCand.ptGev;
    if(l1MuonCand.hwQual >= qualityCut && candPt >=  ptL1Cut) {
      distribution->Fill(l1MuonCand.likelihood, l1MuonCand.refLayer);
      LogTrace("l1tMuBayesEventPrint") <<" LikelihoodDistribution::fill likelihood "<<l1MuonCand.likelihood<<" refLayer "<<l1MuonCand.refLayer <<std::endl;
    }
  }
}

} /* namespace L1MuAn */
