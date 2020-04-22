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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
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

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
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
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/MuCorrelatorConfig.h"
#include "L1Trigger/L1TMuonBayes/plugins/L1TMuonBayesMuCorrelatorTrackProducer.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"


#include "usercode/L1MuonAnalyzer/interface/EfficiencyAnalyser.h"
#include "usercode/L1MuonAnalyzer/interface/RateAnalyser.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>
#include <stddef.h>

namespace L1MuAn {

class L1MuonAnalyzerOmtf: public edm::EDAnalyzer {
public:
  explicit L1MuonAnalyzerOmtf(const edm::ParameterSet& conf);
  virtual ~L1MuonAnalyzerOmtf();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void analyzeEfficiency(const edm::Event&, const edm::EventSetup&);
  virtual void analyzeRate(const edm::Event&, const edm::EventSetup&);

  virtual void endJob();

  double getDeltaR(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand);

  bool matched(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand);

private:
  std::string analysisType;

  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken;
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection > omtfToken;

  std::vector<std::unique_ptr<EfficiencyAnalyser> > omtfEfficiencyAnalysers;

  std::vector<std::unique_ptr<EfficiencyAnalyser> > omtfNNEfficiencyAnalysers;


  std::vector<std::unique_ptr<RateAnalyser> > omtfRateAnalysers;

  std::vector<std::unique_ptr<RateAnalyser> > omtfNNRateAnalysers;


  double maxDeltaR = 0.3;

  double etaCutFrom = 0.82;
  double etaCutTo = 1.24;

  vector<double> nn_pThresholds;

  double hwPtToPtGeV(int hwPt) {
    return ( (hwPt - 1.)/2.);
  }

  TH1* candPerEvent = nullptr;

  TH1* firedPlanesEventCntOmtf = nullptr;
  TH1* firedPlanesEventCntNN = nullptr;
 };

}
#endif /* PLUGINS_L1MUONANALYZEROMTF_H_ */
