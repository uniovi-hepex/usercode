/*
 * L1MuonAnalyzerOmtf.h
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#ifndef PLUGINS_L1MUONANALYZEROMTF_H_
#define PLUGINS_L1MUONANALYZEROMTF_H_

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "DataFormats/L1TMuon/interface/BayesMuCorrelatorTrack.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

#include "usercode/L1MuonAnalyzer/interface/MuonMatcher.h"
#include "usercode/L1MuonAnalyzer/interface/EfficiencyAnalyser.h"
#include "usercode/L1MuonAnalyzer/interface/RateAnalyser.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>
#include <stddef.h>
#include "boost/dynamic_bitset.hpp"


namespace L1MuAn {

class L1MuonAnalyzerOmtf: public edm::EDAnalyzer {
public:
  explicit L1MuonAnalyzerOmtf(const edm::ParameterSet& conf);
  virtual ~L1MuonAnalyzerOmtf();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void analyzeEfficiency(const edm::Event&, std::vector<MatchingResult>& matchingResults);

  //virtual void analyzeRate(const edm::Event& event, std::vector<MatchingResult>& matchingResults, const edm::SimVertexContainer* simVertices);
  virtual void analyzeRate(const edm::Event& event, std::vector<MatchingResult>& matchingResults);

  virtual void analyzeRate(const edm::Event&, const edm::EventSetup&);

  virtual void endJob();

  //simplified ghost busting
  //only candidates in the bx=0 are included
  std::vector<const l1t::RegionalMuonCand*> ghostBust(const l1t::RegionalMuonCandBxCollection* mtfCands);

  double getDeltaR(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand);

  bool matched(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand);

private:
  std::string analysisType;

  bool matchUsingPropagation = true;

  bool fillMatcherHists = false;

  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken;

  edm::EDGetTokenT<edm::SimVertexContainer> simVertexesToken;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken;

  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection > omtfToken;


  MuonMatcher muonMatcher;


  std::vector<std::unique_ptr<EfficiencyAnalyser> > omtfEfficiencyAnalysers;

  std::vector<std::unique_ptr<EfficiencyAnalyser> > omtfNNEfficiencyAnalysers;


  std::vector<std::unique_ptr<RateAnalyser> > omtfRateAnalysers;

  std::vector<std::unique_ptr<RateAnalyser> > omtfNNRateAnalysers;


  std::vector<std::unique_ptr<CandsMatchingAnalyser> > omtfCandsMatchingAnalysers;

  std::vector<std::unique_ptr<CandsMatchingAnalyser> > omtfNNCandsMatchingAnalysers;

  double maxDeltaR = 0.3;

  double etaCutFrom = 0.82;
  double etaCutTo = 1.24;

  vector<double> nn_pThresholds;

  double hwPtToPtGeV(int hwPt) {
    return ( (hwPt - 1.)/2.);
  }

  TH1* candPerEvent = nullptr;

  TH1* firedLayersEventCntOmtf = nullptr;
  TH1* firedLayersEventCntNN = nullptr;

  TH1* ptGenHist = nullptr;
 };

}
#endif /* PLUGINS_L1MUONANALYZEROMTF_H_ */
