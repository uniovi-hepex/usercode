/*
 * RateAnalyser.cc
 *
 *  Created on: Apr 1, 2020
 *      Author: kbunkow
 */

#include <usercode/L1MuonAnalyzer/interface/RateAnalyser.h>

namespace L1MuAn {

RateAnalyser::RateAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo): qualityCut(qualityCut) {
  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<name<<"candPt_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand pt, "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candPt = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo);

  histName.str("");
  histTitle.str("");

  histName<<name<<"candEta_PtCut1GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 1 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut1 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<name<<"candEta_PtCut10GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 10 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut10 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<name<<"candEta_PtCut20GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 20 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut20 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

}

RateAnalyser::~RateAnalyser() {
  // TODO Auto-generated destructor stub
}

///center of eta bin
double hwEtaToEta(int hwEta) {
  double etaUnit = 0.010875; //TODO from the interface note, should be defined somewhere globally

  return (hwEta * etaUnit);
}

void RateAnalyser::fill(L1MuonCand& l1MuonCand) {
  auto candPtGev = l1MuonCand.ptGev;
  if(candPtGev >= candPt->GetXaxis()->GetXmax())
    candPtGev = candPt->GetXaxis()->GetXmax() - 0.01;

  if(l1MuonCand.hwQual >= qualityCut) {
    candPt->Fill(candPtGev);

    auto eta = hwEtaToEta(l1MuonCand.hwEta);
    if(candPtGev >= 1.)
      candEtaPtCut1->Fill(eta);

    if(candPtGev >= 10.)
      candEtaPtCut10->Fill(eta);

    if(candPtGev >= 20.)
      candEtaPtCut20->Fill(eta);
  }

}

void RateAnalyser::write() {
  candPt->Write();

  candEtaPtCut1->Write();
  candEtaPtCut10->Write();
  candEtaPtCut20->Write();
}

} /* namespace L1MuAn */
