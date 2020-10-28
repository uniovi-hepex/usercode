/*
 * MuonMatcher.h
 *
 *  Created on: Jun 16, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_MUONMATCHER_H_
#define INTERFACE_MUONMATCHER_H_

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TH2F.h"

class TFileDirectory;

namespace L1MuAn {

double hwGmtPhiToGlobalPhi(int phi);

double foldPhi(double phi);

class MatchingResult {
public:
  enum class ResultType: short {
    propagationFailed = -1,
    notMatched = 0,
    matched = 1,
    duplicate = 2
  };

  MatchingResult() {}

  MatchingResult(const SimTrack& simTrack): simTrack(&simTrack) {
    pdgId = simTrack.type();
    genPt = simTrack.momentum().pt();
    genEta = simTrack.momentum().eta();
    genPhi = simTrack.momentum().phi();
  }

  MatchingResult(const TrackingParticle& trackingParticle): trackingParticle(&trackingParticle) {
    pdgId = trackingParticle.pdgId();
    genPt = trackingParticle.pt();
    genEta = trackingParticle.momentum().eta();
    genPhi = trackingParticle.momentum().phi();
  }

  ResultType result = ResultType::notMatched;
  //bool propagationFailed = false;
  double deltaPhi = 0;
  double deltaEta = 0;

  double propagatedPhi = 0;
  double propagatedEta = 0;

  double matchingLikelihood = 0;

  const l1t::RegionalMuonCand* muonCand = nullptr;

  //to avoid using simTrack or trackingParticle
  double pdgId = 0;
  double genPt = 0;
  double genEta = 0;
  double genPhi = 0;

  const SimTrack* simTrack = nullptr;
  const TrackingParticle* trackingParticle = nullptr;

};

class MuonMatcher {
public:
  MuonMatcher(const edm::ParameterSet& edmCfg);

  void setup(const edm::EventSetup& eventSetup);

  virtual ~MuonMatcher();

  void saveHists();

  FreeTrajectoryState simTrackToFts(const SimTrack& simTrack, const SimVertex& simVertex);

  FreeTrajectoryState simTrackToFts(const TrackingParticle& trackingParticle);

  TrajectoryStateOnSurface atStation2(FreeTrajectoryState ftsStart, float eta) const;

  TrajectoryStateOnSurface propagate(const SimTrack& simTrack, const edm::SimVertexContainer* simVertices);

  TrajectoryStateOnSurface propagate(const TrackingParticle& trackingParticle);

  //tsof should be the result of track propagation
  MatchingResult match(const l1t::RegionalMuonCand* omtfCand, const SimTrack& simTrack, TrajectoryStateOnSurface& tsof);

  MatchingResult match(const l1t::RegionalMuonCand* omtfCand, const TrackingParticle& trackingParticle, TrajectoryStateOnSurface& tsof);

  std::vector<MatchingResult> cleanMatching(std::vector<MatchingResult> matchingResults, std::vector<const l1t::RegionalMuonCand*>& muonCands);

  std::vector<MatchingResult> match(std::vector<const l1t::RegionalMuonCand*>& muonCands, const edm::SimTrackContainer* simTracks, const edm::SimVertexContainer* simVertices,
      std::function<bool(const SimTrack& )> const& simTrackFilter);

  std::vector<MatchingResult> match(std::vector<const l1t::RegionalMuonCand*>& muonCands, const TrackingParticleCollection* trackingParticles,
      std::function<bool(const TrackingParticle& )> const& simTrackFilter);

  //propagations is not used
  void fillHists(std::vector<const l1t::RegionalMuonCand*>& muonCands, const edm::SimTrackContainer* simTracks,  std::function<bool(const SimTrack& )> const& simTrackFilter);


  MatchingResult match(const l1t::RegionalMuonCand* omtfCand, const TrackingParticle& trackingParticle);

  //propagations is not used, but inseted the histograms deltaPhiVertexCand_Mean_pos and deltaPhiVertexCand_StdDev_pos
  std::vector<MatchingResult> matchWithoutPorpagation(std::vector<const l1t::RegionalMuonCand*>& muonCands, const TrackingParticleCollection* trackingParticles,
      std::function<bool(const TrackingParticle& )> const& simTrackFilter);

private:
  //const edm::EventSetup& eventSetup;

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  edm::ESHandle<MagneticField> magField;

  edm::ESHandle<Propagator> propagator;

  TH2F* deltaPhiPropCand = nullptr; //delta phi between propagated track and muon candidate, stores the likelihood
  TH2F* deltaPhiVertexProp = nullptr; //delta phi between phi at vertex and propagated track phi

  TH1D* deltaPhiPropCandMean = nullptr;
  TH1D* deltaPhiPropCandStdDev = nullptr;

  TH1D* ptGen_pos = nullptr;
  TH1D* deltaPhiVertexCand_Mean_pos = nullptr;
  TH1D* deltaPhiVertexCand_StdDev_pos = nullptr;

  TH1D* ptGen_neg = nullptr;
  TH1D* deltaPhiVertexCand_Mean_neg = nullptr;
  TH1D* deltaPhiVertexCand_StdDev_neg = nullptr;

  bool fillMean = false;
  bool matchUsingPropagation = true;
};

}
#endif /* INTERFACE_MUONMATCHER_H_ */
