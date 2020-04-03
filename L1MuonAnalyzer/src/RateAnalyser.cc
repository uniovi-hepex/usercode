/*
 * RateAnalyser.cc
 *
 *  Created on: Apr 1, 2020
 *      Author: kbunkow
 */

#include <usercode/L1MuonAnalyzer/interface/RateAnalyser.h>

namespace L1MuAn {

RateAnalyser::RateAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo) {
  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<name<<"_candPt_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand pt, "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candPt = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo);

  histName.str("");
  histTitle.str("");

  histName<<name<<"_candEta_PtCut1GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 1 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut1 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);


  histName<<name<<"_candEta_PtCut10GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 10 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut10 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);


  histName<<name<<"_candEta_PtCut20GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 20 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut10 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

}

RateAnalyser::~RateAnalyser() {
  // TODO Auto-generated destructor stub
}

void RateAnalyser::fill(L1MuonCand& l1MuonCand) {

}

} /* namespace L1MuAn */
