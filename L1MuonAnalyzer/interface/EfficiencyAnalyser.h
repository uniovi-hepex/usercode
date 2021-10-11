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
  EfficiencyAnalyser(bool useUpt);
  virtual ~EfficiencyAnalyser();

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand) = 0;

  virtual void write() = 0;

protected:
  bool useUpt = false;
};

class PtGenVsPtCand : public EfficiencyAnalyser{
public:
  PtGenVsPtCand(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, int nBins, double binsFrom, double binsTo, bool useUpt = false);

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand);

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
  EfficiencyVsPhi(TFileDirectory& subDir, std::string name, double etaFrom, double etaTo, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt = false);

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand);

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
  EfficiencyVsEta(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt = false);

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand);

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

class EfficiencyPtGenVsDxy : public EfficiencyAnalyser {
public:
  EfficiencyPtGenVsDxy(TFileDirectory& subDir, std::string name, int qualityCut, double ptL1Cut, int nBinsPt, int nBinsDxy, bool ifPtBelowCut, bool useUpt = false);

  virtual ~EfficiencyPtGenVsDxy();

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand);

  virtual void write();
private:
  int qualityCut = 0;

  double ptL1Cut = 0;

  bool ifPtBelowCut  = false;

  TH2* allCands = nullptr;
  TH2* aceptedCands = nullptr;

  TH2* efficiency = nullptr;
};

class LikelihoodDistribution : public EfficiencyAnalyser {
public:
  LikelihoodDistribution(TFileDirectory& subDir, std::string name, int qualityCut, double ptGenCut, double ptL1Cut, int nBins, bool useUpt = false);

  virtual void fill(double ptGen, double etaGen, double phiGen, double dxyGen, L1MuonCand& l1MuonCand);

  virtual void write() {
    distribution->Write();
  }

private:
  int qualityCut = 0;

  double ptGenCut = 0;
  double ptL1Cut = 0;

  TH2* distribution = nullptr;
};


} /* namespace L1MuAn */

#endif /* INTERFACE_EFFICIENCYANALYSER_H_ */
