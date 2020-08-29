/*
 * RateAnalyser.h
 *
 *  Created on: Apr 1, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_RATEANALYSER_H_
#define INTERFACE_RATEANALYSER_H_

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "usercode/L1MuonAnalyzer/interface/L1MuonCand.h"
#include "usercode/L1MuonAnalyzer/interface/EfficiencyAnalyser.h"

#include "TH2.h"



namespace L1MuAn {

class RateAnalyser {
public:
  RateAnalyser(TFileDirectory subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo);
  virtual ~RateAnalyser();

  virtual void fill(L1MuonCand& l1MuonCand);

  virtual void write();

private:
  int qualityCut = 0;


  TH1* candPt = nullptr;

  TH1* candEtaPtCut1 = nullptr;
  TH1* candEtaPtCut10 = nullptr;
  TH1* candEtaPtCut20 = nullptr;

};

class CandsMatchingAnalyser {
public:
  CandsMatchingAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo);

  virtual ~CandsMatchingAnalyser() {};

  //virtual void fill(L1MuonCand& l1MuonCand, const SimTrack* matchedSimTrack, const edm::SimVertexContainer* simVertices);


  virtual void fill(L1MuonCand& l1MuonCand, const TrackingParticle* matchedTrackingParticle);


  virtual void write();

  class CandsMatchingHists {
  public:
    TFileDirectory subDir;
    RateAnalyser rateAn;
    PtGenVsPtCand ptGenVsPtCand;
    TH2* simVertexRhoVsPtGen = nullptr;

    CandsMatchingHists(TFileDirectory& parrentSubDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo);
    virtual ~CandsMatchingHists() {}

    virtual void fill(L1MuonCand& l1MuonCand, const TrackingParticle* matchedTrackingParticle, double simVertexRho);
    virtual void write();
  };

private:
  CandsMatchingHists promptMuons;
  CandsMatchingHists nonPromptMuons;
  CandsMatchingHists muonsFromPions;
  CandsMatchingHists muonsFromKaons;

  RateAnalyser notMatchedRateAn;

  //PtGenVsPtCand etaCand12PtGenVsPtCand;

  TH2* simVertexRhoVsPtGen = nullptr;

  LikelihoodDistribution likelihoodDistribution;
};

} /* namespace L1MuAn */

#endif /* INTERFACE_RATEANALYSER_H_ */
