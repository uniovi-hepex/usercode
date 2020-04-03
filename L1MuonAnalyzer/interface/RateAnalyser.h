/*
 * RateAnalyser.h
 *
 *  Created on: Apr 1, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_RATEANALYSER_H_
#define INTERFACE_RATEANALYSER_H_

#include "usercode/L1MuonAnalyzer/interface/L1MuonCand.h"
#include "TH1.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

namespace L1MuAn {

class RateAnalyser {
public:
  RateAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo);
  virtual ~RateAnalyser();

  virtual void fill(L1MuonCand& l1MuonCand);

  virtual void write() = 0;

private:
  int qualityCut = 0;


  TH1* candPt = nullptr;

  TH1* candEtaPtCut1 = nullptr;
  TH1* candEtaPtCut10 = nullptr;
  TH1* candEtaPtCut22 = nullptr;

};

} /* namespace L1MuAn */

#endif /* INTERFACE_RATEANALYSER_H_ */
