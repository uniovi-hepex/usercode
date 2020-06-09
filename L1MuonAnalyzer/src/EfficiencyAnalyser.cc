/*
 * EfficiencyAnalyser.cpp
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/interface/EfficiencyAnalyser.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

namespace L1MuAn {

EfficiencyAnalyser::EfficiencyAnalyser() {
  // TODO Auto-generated constructor stub

}

EfficiencyAnalyser::~EfficiencyAnalyser() {
  // TODO Auto-generated destructor stub
}


PtGenVsPtCand::PtGenVsPtCand(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, int nBins, double binsFrom, double binsTo): etaFrom(etaFrom), etaTo(etaTo), qualityCut(qualityCut) {
  std::ostringstream histName;

  histName<<name<<"_ptGenVsPtCand_eta_"<<etaFrom<<"_"<<etaTo<<"_qualityCut_"<<qualityCut;

  std::ostringstream histTitle;
  histTitle<<name<<" ptGenVsPtCand eta: "<<etaFrom<<" - "<<etaTo<<" quality >= "<<qualityCut<<";ptGen [GeV]; pt L1 cand [GeV]";

  ptGenVsPtCand = subDir.make<TH2D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo, nBins, binsFrom, binsTo);

}

void PtGenVsPtCand::fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand) {
  if(ptGen >= ptGenVsPtCand->GetXaxis()->GetXmax())
    ptGen = ptGenVsPtCand->GetXaxis()->GetXmax() - 0.01;

  if( (etaFrom <= fabs(etaGen) ) && (fabs(etaGen) <= etaTo) ) {
    if(l1MuonCand.hwQual >= qualityCut) {

      double candPt = l1MuonCand.ptGev;
      if(candPt >= ptGenVsPtCand->GetYaxis()->GetXmax())
        candPt = ptGenVsPtCand->GetYaxis()->GetXmax() - 0.01;

      ptGenVsPtCand->Fill(ptGen, candPt);
    }
    else
      ptGenVsPtCand->Fill(ptGen, 0);
  }
}






EfficiencyVsPhi::EfficiencyVsPhi(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, double ptGenCut, double ptL1Cut, int nBins):
    etaFrom(etaFrom), etaTo(etaTo), qualityCut(qualityCut), ptGenCut(ptGenCut), ptL1Cut(ptL1Cut) {
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

void EfficiencyVsPhi::fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand) {
  if( (etaFrom <= fabs(etaGen) ) && (fabs(etaGen) <= etaTo) && ptGen >= ptGenCut ) {
    allCands->Fill(phiGen);
    if(l1MuonCand.hwQual >= qualityCut && l1MuonCand.ptGev >=  ptL1Cut) {
      aceptedCands->Fill(phiGen);
    }
  }
}




EfficiencyVsEta::EfficiencyVsEta(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins):
  qualityCut(qualityCut), ptGenCut(ptGenCut), ptL1Cut(ptL1Cut) {
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

void EfficiencyVsEta::fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand) {
  if( ptGen >= ptGenCut ) {
    allCands->Fill(etaGen);
    if(l1MuonCand.hwQual >= qualityCut && l1MuonCand.ptGev >=  ptL1Cut) {
      aceptedCands->Fill(etaGen);
    }
  }
}


} /* namespace L1MuAn */
