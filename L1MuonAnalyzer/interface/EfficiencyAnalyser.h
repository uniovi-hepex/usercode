/*
 * EfficiencyAnalyser.h
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_EFFICIENCYANALYSER_H_
#define INTERFACE_EFFICIENCYANALYSER_H_

#include "usercode/L1MuonAnalyzer/interface/L1MuonCand.h"
#include "TH2.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


namespace L1MuAn {

class EfficiencyAnalyser {
public:
  EfficiencyAnalyser();
  virtual ~EfficiencyAnalyser();

  virtual void fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand) = 0;

  virtual void write() = 0;
private:

};

class PtGenVsPtCand : public EfficiencyAnalyser{
public:
  PtGenVsPtCand(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, int nBins, double binsFrom, double binsTo);

  virtual void fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand);

  virtual void write() {
    ptGenVsPtCand->Write();
  }

private:
  double etaFrom = 0;
  double etaTo = 0;

  int qualityCut = 0;

  TH2* ptGenVsPtCand = nullptr;
};

class EfficiencyVsPhi : public EfficiencyAnalyser{
public:
  EfficiencyVsPhi(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, double ptGenCut, double ptL1Cut, int nBins);

  virtual void fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand);

  virtual void write() {
    allCands->Write();
    aceptedCands->Write();
  }

private:
  double etaFrom = 0;
  double etaTo = 0;

  int qualityCut = 0;

  double ptGenCut = 0;
  double ptL1Cut = 0;

  TH1* allCands = nullptr;
  TH1* aceptedCands = nullptr;
};


class EfficiencyVsEta : public EfficiencyAnalyser{
public:
  EfficiencyVsEta(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins);

  virtual void fill(double ptGen, double etaGen, double phiGen, L1MuonCand& l1MuonCand);

  virtual void write() {
    allCands->Write();
    aceptedCands->Write();
  }

private:
  int qualityCut = 0;

  double ptGenCut = 0;
  double ptL1Cut = 0;

  TH1* allCands = nullptr;
  TH1* aceptedCands = nullptr;
};


} /* namespace L1MuAn */

#endif /* INTERFACE_EFFICIENCYANALYSER_H_ */
