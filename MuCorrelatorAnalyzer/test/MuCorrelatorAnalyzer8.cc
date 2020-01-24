//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

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

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/MuCorrelatorConfig.h"
#include "L1Trigger/L1TMuonBayes/plugins/L1TMuonBayesMuCorrelatorTrackProducer.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>
#include <stddef.h>

#define EDM_ML_LOGDEBUG

using namespace std;

std::ostream & operator<< (std::ostream &out, const l1t::BayesMuCorrelatorTrack&  muCand);

class TriggerAlgo {
public:
  TriggerAlgo(std::string name, double ptCut): name(name), ptCut(ptCut) {};
  virtual ~TriggerAlgo() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack) = 0;

  std::string name;

  double ptCut = 0;

  int L1Tk_nPar = 4; //todo take from config
};


class AllTTTRacks: public TriggerAlgo {
public:
  AllTTTRacks(double ptCut): TriggerAlgo("AllTTTRacks" + std::to_string((int)ptCut), ptCut) {};
  virtual ~AllTTTRacks() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
      return true;
  }
};

class AllTTTRacksBarrel: public TriggerAlgo {
public:
  AllTTTRacksBarrel(double ptCut): TriggerAlgo("AllTTTRacksBarrel" + std::to_string((int)ptCut), ptCut) {};
  virtual ~AllTTTRacksBarrel() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if(abs(muCorrelatorTrack.getEta() ) < 0.82)
      return true;
    return false;
  }
};

class SingleMuAlgo: public TriggerAlgo {
public:
  SingleMuAlgo(double ptCut): TriggerAlgo("SingleMuAlgo" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgo() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if(muCorrelatorTrack.hwQual() >= 12 && muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack)
      return true;
    return false;
  }
};

class SingleMuAlgoBarrel: public TriggerAlgo {
public:
  SingleMuAlgoBarrel(double ptCut): TriggerAlgo("SingleMuAlgoBarrel" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoBarrel() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if(muCorrelatorTrack.hwQual() >= 12 && muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        abs(muCorrelatorTrack.getEta() ) < 0.82)
      return true;
    return false;
  }
};

class SingleMuAlgoOverlap: public TriggerAlgo {
public:
  SingleMuAlgoOverlap(double ptCut): TriggerAlgo("SingleMuAlgoOverlap" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoOverlap() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if(muCorrelatorTrack.hwQual() >= 12 && muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        abs(muCorrelatorTrack.getEta() ) >= 0.82 && abs(muCorrelatorTrack.getEta() ) < 1.24 )
      return true;
    return false;
  }
};

class SingleMuAlgoEndcap: public TriggerAlgo {
public:
  SingleMuAlgoEndcap(double ptCut): TriggerAlgo("SingleMuAlgoEndcap" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoEndcap() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if(muCorrelatorTrack.hwQual() >= 12 && muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        abs(muCorrelatorTrack.getEta() ) >= 1.24 )
      return true;
    return false;
  }
};


class SingleMuAlgoSoftCuts: public TriggerAlgo {
public:
  SingleMuAlgoSoftCuts(double ptCut): TriggerAlgo("SingleMuAlgoSoftCuts" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoSoftCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1000) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1100) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) &&
        ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
            ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
              muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
           muCorrelatorTrack.getTtTrackPtr().isNull() )
    )
      return true;
    return false;
  }
};

class SingleMuAlgoSoftCutsOverlap: public TriggerAlgo {
public:
  SingleMuAlgoSoftCutsOverlap(double ptCut): TriggerAlgo("SingleMuAlgoSoftCutsOverlap" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoSoftCutsOverlap() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1000) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1100) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) &&
        ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
            ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
              muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
           muCorrelatorTrack.getTtTrackPtr().isNull() ) &&
        (abs(muCorrelatorTrack.getEta() ) >= 0.82 && abs(muCorrelatorTrack.getEta() ) < 1.24 )
    )
      return true;
    return false;
  }
};

//the same cuts as in the L1TMuonBayesMuCorrelatorTrackProducer::produce
class SingleMuAlgoSoftCuts1: public TriggerAlgo {
public:
  SingleMuAlgoSoftCuts1(double ptCut): TriggerAlgo("SingleMuAlgoSoftCuts1_ptCut" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoSoftCuts1() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1100) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1400) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) &&
        ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 200 ) || muCorrelatorTrack.getTtTrackPtr().isNull() )
    )
      return true;
    return false;
  }
};

class SingleMuAlgoSoftCutsOverlap1: public TriggerAlgo {
public:
  SingleMuAlgoSoftCutsOverlap1(double ptCut): TriggerAlgo("SingleMuAlgoSoftCutsOverlap1_ptCut" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoSoftCutsOverlap1() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1100) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1400) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) &&
        ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 200 ) || muCorrelatorTrack.getTtTrackPtr().isNull() ) &&
        (abs(muCorrelatorTrack.getEta() ) >= 0.82 && abs(muCorrelatorTrack.getEta() ) < 1.24 )
    )
      return true;
    return false;
  }
};

class SingleMuAlgoPdfSumSoftCuts: public TriggerAlgo {
public:
  SingleMuAlgoPdfSumSoftCuts(double ptCut): TriggerAlgo("SingleMuAlgoPdfSumSoftCuts" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoPdfSumSoftCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1000) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1100) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) //&&
        //( (muCorrelatorTrack.getTtTrackPtr().isNonnull() && muCorrelatorTrack.getTtTrackPtr()->getChi2Red(L1Tk_nPar) < 200 ) || muCorrelatorTrack.getTtTrackPtr().isNull() )
    )
      return true;
    return false;
  }
};

class SingleMuAlgoHardCuts: public TriggerAlgo {
public:
  SingleMuAlgoHardCuts(double ptCut): TriggerAlgo("SingleMuAlgoHardCuts" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgoHardCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack){
    if( muCorrelatorTrack.hwQual() >= 12 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::fastTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1300) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1900) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 4) &&
           ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
               ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
                 muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
              muCorrelatorTrack.getTtTrackPtr().isNull() )
       )
      return true;
    return false;
  }
};

class HscpAlgo: public TriggerAlgo {
public:
  HscpAlgo(double ptCut): TriggerAlgo("HscpAlgo"+ std::to_string((int)ptCut), ptCut) {};
  virtual ~HscpAlgo() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack) {
    if(muCorrelatorTrack.hwQual() >= 12 && muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::slowTrack)
      return true;
    return false;
  }
};

class HscpAlgoHardCuts: public TriggerAlgo {
public:
  HscpAlgoHardCuts(double ptCut): TriggerAlgo("HscpAlgoHardCuts"+ std::to_string((int)ptCut), ptCut) {};
  virtual ~HscpAlgoHardCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack) {
    if( muCorrelatorTrack.hwQual() >= 13 &&
        muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::slowTrack &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1800 && muCorrelatorTrack.getBetaLikelihood() >= 9) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 4 && muCorrelatorTrack.pdfSum() > 2000 && muCorrelatorTrack.getBetaLikelihood() >= 10) ||
         muCorrelatorTrack.getFiredLayerCnt() >= 5  ) &&
         ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
             ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
               muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
            muCorrelatorTrack.getTtTrackPtr().isNull() )
      )
    {
      return true;
    }
    return false;
  }
};

class HscpAlgoSoftCuts: public TriggerAlgo {
public:
  HscpAlgoSoftCuts(double ptCut): TriggerAlgo("HscpAlgoSoftCuts"+ std::to_string((int)ptCut), ptCut) {};
  virtual ~HscpAlgoSoftCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack) {
    if( muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::slowTrack &&
        muCorrelatorTrack.hwQual() >= 13 &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1000 && muCorrelatorTrack.getBetaLikelihood() >= 6) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1100 && muCorrelatorTrack.getBetaLikelihood() >= 7) ||
          (muCorrelatorTrack.getFiredLayerCnt() == 4 && muCorrelatorTrack.pdfSum() > 1200 && muCorrelatorTrack.getBetaLikelihood() >= 9) ||
           muCorrelatorTrack.getFiredLayerCnt() >= 5) &&
       ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
           ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
             muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
          muCorrelatorTrack.getTtTrackPtr().isNull() )
      )
    {
      return true;
    }
    return false;
  }
};

class HscpAlgoPdfSumCuts: public TriggerAlgo {
public:
  HscpAlgoPdfSumCuts(double ptCut): TriggerAlgo("HscpAlgoPdfSumCuts"+ std::to_string((int)ptCut), ptCut) {};
  virtual ~HscpAlgoPdfSumCuts() {};

  virtual bool accept(const l1t::BayesMuCorrelatorTrack& muCorrelatorTrack) {
    if( muCorrelatorTrack.getCandidateType() == l1t::BayesMuCorrelatorTrack::slowTrack &&
        muCorrelatorTrack.hwQual() >= 13 &&
        ( (muCorrelatorTrack.getFiredLayerCnt() == 2 && muCorrelatorTrack.pdfSum() > 1000) || // && muCorrelatorTrack.getBetaLikelihood() >= 6
          (muCorrelatorTrack.getFiredLayerCnt() == 3 && muCorrelatorTrack.pdfSum() > 1100) || // && muCorrelatorTrack.getBetaLikelihood() >= 7
          (muCorrelatorTrack.getFiredLayerCnt() == 4 && muCorrelatorTrack.pdfSum() > 1200) || // && muCorrelatorTrack.getBetaLikelihood() >= 9
           muCorrelatorTrack.getFiredLayerCnt() >= 5) &&
       ( (muCorrelatorTrack.getTtTrackPtr().isNonnull() &&
           ((muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() == 4 && muCorrelatorTrack.getTtTrackPtr()->getChi2(L1Tk_nPar) < 100 ) ||
             muCorrelatorTrack.getTtTrackPtr()->getStubRefs().size() > 4) ) ||
          muCorrelatorTrack.getTtTrackPtr().isNull() )
      )
    {
      return true;
    }
    return false;
  }
};


///////////////////////////////////


std::ostream & operator<< (std::ostream &out, const l1t::BayesMuCorrelatorTrack&  muCand) {
  out
  <<" hwPt "<<muCand.hwPt()
  <<" pt "<<muCand.getPt()<<" GeV "
  <<" eta "<<muCand.getEta()
  <<" phi "<<muCand.getPhi()
  <<" hwSign "<<muCand.hwSign()
  <<" quality "<<muCand.hwQual()
  <<" pdfSum "<< muCand.pdfSum()
  <<" type "<<muCand.getCandidateType()
  <<" beta "<<muCand.getBeta()
  <<" betaLikelihood "<<muCand.getBetaLikelihood()
  <<" "<<muCand.getFiredLayerBits(30);

  return out;
}

//stupid but otherwise does not work with LogTrace
std::string toString(const l1t::BayesMuCorrelatorTrack&  muCand) {
  std::ostringstream ostr;
  ostr<<muCand;
  return ostr.str();
}

auto printTrackigParticleShort(const edm::Ptr< TrackingParticle >& tpPtr) {
  std::ostringstream ostr;
  ostr<< "Tracking particle  pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
              << " eta: " <<setw(7)<< tpPtr->eta()
              << " phi: " <<setw(7)<< tpPtr->phi()
              << " beta: " <<setw(7)<< tpPtr->p4().Beta()
              << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
              << " bx: "<<tpPtr->eventId().bunchCrossing()
              <<dec<<" key "<<tpPtr.key();
  return ostr.str();
};

auto printTrackigParticleShort(const edm::Ptr< TrackingParticle >& tpPtr, const edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& MCTruthTTTrackHandle) {
  std::ostringstream ostr;
  ostr<< "Tracking particle  pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
              << " eta: " <<setw(7)<< tpPtr->eta()
              << " phi: " <<setw(7)<< tpPtr->phi()
              << " beta: " <<setw(7)<< tpPtr->p4().Beta()
              << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
              << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr).size()
              <<" bx: "<<tpPtr->eventId().bunchCrossing()
              <<dec<<" key "<<tpPtr.key();
  return ostr.str();
};

auto printTTTRack(const edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > >& ttTrackPtr, bool genuine, bool loosegenuine, int L1Tk_nPar) {
  std::ostringstream ostr;
  ostr << "matched ttTrack       pt: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).perp()
                   <<" charge: "<<setw(2)<<(ttTrackPtr->getRInv() > 0 ? 1 : -1)
                   << " eta: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).eta()
                   << " phi: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).phi()
                   << " chi2: " <<setw(7)<< ttTrackPtr->getChi2(L1Tk_nPar)
                   << " consistency: " <<setw(7)<< ttTrackPtr->getStubPtConsistency(L1Tk_nPar)
                   << " z0: " <<setw(7)<< ttTrackPtr->getPOCA(L1Tk_nPar).z()
                   << " nstub: " <<setw(7)<< ttTrackPtr->getStubRefs().size()
                   << (genuine  ? " genuine " : "")
                   << (loosegenuine ? " loose genuine " : "");
  return ostr.str();
};



class AnalyserBase {
public:
  AnalyserBase(std::shared_ptr<TriggerAlgo>& triggerAlgo): triggerAlgo(triggerAlgo) {

  }

  virtual ~AnalyserBase() {}

  virtual void reset() {
    ptOfBestL1MuCand = -1;
    bestL1MuCand = nullptr;
    notAcceptedL1MuCand = nullptr;
    numberOfAcceptedCandidates = 0;
  }

  virtual void takeCanidate(const l1t::BayesMuCorrTrackBxCollection::const_iterator& itL1MuCand);

protected:
  std::shared_ptr<TriggerAlgo> triggerAlgo;

  l1t::BayesMuCorrelatorTrack const* bestL1MuCand = nullptr;
  l1t::BayesMuCorrelatorTrack const* notAcceptedL1MuCand = nullptr;

  double ptOfBestL1MuCand = -1;

  int numberOfAcceptedCandidates = 0;
};

void AnalyserBase::takeCanidate(const l1t::BayesMuCorrTrackBxCollection::const_iterator& itL1MuCand) {
  if(triggerAlgo->accept(*itL1MuCand) ) {
    if(ptOfBestL1MuCand < itL1MuCand->getPt() ) {
      if(ptOfBestL1MuCand > 0)
        edm::LogImportant("l1tMuBayesEventPrint") <<"L:"<<__LINE__<<" AnalyserBase::takeCanidate: another track already exists with pt "<<ptOfBestL1MuCand<<" new track has pt "<<itL1MuCand->getPt() <<std::endl  ;
      numberOfAcceptedCandidates++;
      ptOfBestL1MuCand = itL1MuCand->getPt();
      bestL1MuCand = &(*itL1MuCand);
      //LogTrace("l1tMuBayesEventPrint")<<" itBestL1MuCand set"<<endl;
    }
  }
  else {
    notAcceptedL1MuCand = &(*itL1MuCand);
  }
}

class EfficiencyAnalyser: public AnalyserBase {
public:
    //takes the ownership of the triggerAlgo
  //EfficiencyAnalyser() {}

  EfficiencyAnalyser(std::shared_ptr<TriggerAlgo>& triggerAlgo, double ptGenFrom, double ptGenTo, edm::Service<TFileService>& fs): AnalyserBase(triggerAlgo), ptGenFrom(ptGenFrom), ptGenTo(ptGenTo) {
    const int ptBins = 1000;
    const int etaBins = 240;//96;
    const int phiBins = 360;

    TFileDirectory subDir = fs->mkdir( ("EfficiencyAnalyser_" + triggerAlgo->name + "_ptGenFrom_" + std::to_string( (int)ptGenFrom) + "_ptGenTo_" + std::to_string( (int)ptGenTo)).c_str());
    //versus gen (tracking particle) pt, eta, phi
    muCandGenPtMuons = subDir.make<TH1I>("muCandGenPtMuons", "muCandGenPtMuons; pt gen; #events", ptBins, 0., 500.);
    muCandGenEtaMuons = subDir.make<TH1D>("muCandGenEtaMuons", "muCandGenEtaMuons; eta gen; #events", etaBins, -2.4, 2.4);



    ostringstream title;
    title<<" ptTTTrack >= " <<triggerAlgo->ptCut<<" GeV"<<", "<<ptGenFrom<<" < ptGen < "<<ptGenTo <<" GeV, event 0 ";
    title<<"; eta; #events";

    gpMuonGenEtaMuons_withPtCuts = subDir.make<TH1D>("gpMuonGenEtaMuons_withPtCuts", ("gpMuonGenEtaMuons_withPtCuts " + title.str()).c_str(), etaBins, -2.4, 2.4);
    ttMuonGenEtaMuons_withPtCuts = subDir.make<TH1D>("ttMuonGenEtaMuons_withPtCuts", ("ttMuonGenEtaMuons_withPtCuts " + title.str()).c_str(), etaBins, -2.4, 2.4);
    muCandGenEtaMuons_withPtCuts = subDir.make<TH1D>("muCandGenEtaMuons_withPtCuts", ("muCandGenEtaMuons_withPtCuts " + title.str()).c_str(), etaBins, -2.4, 2.4);

    muCandGenPhiMuons = subDir.make<TH1D>("muCandGenPhiMuons", "muCandGenPhiMuons; phi gen; #events", phiBins, -M_PI, M_PI);

    ptGenPtMuCandMuonsEv0 = subDir.make<TH2I>("ptGenPtMuCandMuonsEv0", "ptGenPtMuCandMuonsEv0; pT gen [GeV]; pT ttTrack [GeV]; #", 100, 0, 100, 100, 0, 100);
    ptGenPtMuCandMuonsPu = subDir.make<TH2I>("ptGenPtMuCandMuonsPu", "ptGenPtMuCandMuonsPu; pT gen [GeV]; pT ttTrack [GeV]; #", 100, 0, 100, 100, 0, 100);


    ptGenPtMuCandMuonsEv0Barrel = subDir.make<TH2I>("ptGenPtMuCandMuonsEv0Barrel", "ptGenPtMuCandMuonsEv0Barrel |#eta| < 0.82; generated p_{T} [GeV]; ttTrack p_{T} [GeV]; #", 100, 0, 100, 100, 0, 100);
    ptGenPtMuCandMuonsEv0Overlap = subDir.make<TH2I>("ptGenPtMuCandMuonsEv0Overlap", "ptGenPtMuCandMuonsEv0Overlap 0.82 < |#eta| < 1.24; generated p_{T} [GeV]; ttTrack p_{T} [GeV]; #", 100, 0, 100, 100, 0, 100);
    ptGenPtMuCandMuonsEv0Endcap = subDir.make<TH2I>("ptGenPtMuCandMuonsEv0Endcap", "ptGenPtMuCandMuonsEv0Endcap 1.24 < |#eta|; generated p_{T} [GeV]; ttTrack p_{T} [GeV]; #", 100, 0, 100, 100, 0, 100);

    muonsPdfSumFiredPlanes = subDir.make<TH2I>("muonsPdfSumFiredPlanes", "muonsPdfSumFiredPlanes; pdfSum; FiredPlanes; #", 100, 0, 10000, 19, -0.5, 18.5);

    betaLikelihoodFiredPlanesMuons = subDir.make<TH2I>("betaLikelihoodFiredPlanesMuons", "betaLikelihoodFiredPlanesMuons; betaLikelihood; firedPlanes; #", 100, 0, 100, 19, -0.5, 18.5);
    muCandBetaMuons = subDir.make<TH1I>("muCandBetaMuons", "muCandBetaMuons; beta measured; #events", 22, 0., 1.1);
    betaGenBetaL1Mu = subDir.make<TH2I>("betaGenBetaL1Mu", "betaGenBetaL1Mu, staus; betaGen; Beta L1MuCand; #", 20, 0., 1., 42, -1., 1.1);

    lostTtMuonPt = subDir.make<TH1D>("lostTtMuonPt", "lostTtMuonPt; ttMuonPt [GeV]; #events", ptBins, 0., 500.);;
    lostTtMuonEta_ptGen20GeV = subDir.make<TH1D>("lostTtMuonEta_ptGen20GeV", "lostTtMuonEta_ptGen20GeV; eta; #events", etaBins, -2.4, 2.4);
    lostTtMuonPhi_ptGen20GeV = subDir.make<TH1D>("lostTtMuonPhi_ptGen20GeV", "lostTtMuonPhi_ptGen20GeV; phi; #events", phiBins, -M_PI, M_PI);

    acceptedCandidatesVsPtGen = subDir.make<TH2I>("acceptedCandidatesVsPtGen", "number of accepted candidates; ptGen [GeV]; accepted candidates", 50, 0., 100., 5, -0.5, 4.5);
  }


  virtual ~EfficiencyAnalyser() {}

  void fillHistos(const edm::Event& event, edm::Ptr< TrackingParticle >& trackParticle, edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > bestMatchedTTTrack);

private:
  double ptGenFrom = 0;
  double ptGenTo = 100000;

  TH1I* muCandGenPtMuons = nullptr;


  TH1D* muCandGenEtaMuons = nullptr;

  TH1D* gpMuonGenEtaMuons_withPtCuts = nullptr;
  TH1D* ttMuonGenEtaMuons_withPtCuts = nullptr;
  TH1D* muCandGenEtaMuons_withPtCuts = nullptr;

  //TH1D* muCandGenEtaMuons_ptGenCutm2_ptTTCut = nullptr;


  TH1D* muCandGenPhiMuons = nullptr;

  TH2I* ptGenPtMuCandMuonsEv0 = nullptr;
  TH2I* ptGenPtMuCandMuonsPu = nullptr;

  TH2I* ptGenPtMuCandMuonsEv0Barrel = nullptr;
  TH2I* ptGenPtMuCandMuonsEv0Overlap = nullptr;
  TH2I* ptGenPtMuCandMuonsEv0Endcap = nullptr;

  TH2I* muonsPdfSumFiredPlanes = nullptr;

  TH2I* betaLikelihoodFiredPlanesMuons = nullptr;
  TH1I* muCandBetaMuons = nullptr;
  TH2I* betaGenBetaL1Mu = nullptr;

  TH1D* lostTtMuonPt = nullptr;
  TH1D* lostTtMuonEta_ptGen20GeV = nullptr;
  TH1D* lostTtMuonPhi_ptGen20GeV = nullptr;

  TH2I* acceptedCandidatesVsPtGen = nullptr;
};

void EfficiencyAnalyser::fillHistos(const edm::Event& event, edm::Ptr< TrackingParticle >& trackParticle, edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > bestMatchedTTTrack) {
  acceptedCandidatesVsPtGen->Fill(trackParticle->pt(), numberOfAcceptedCandidates);

  if(trackParticle->eventId().event() == 0 &&
      trackParticle->pt() >= ptGenFrom && trackParticle->pt() <= ptGenTo)
  {
    gpMuonGenEtaMuons_withPtCuts->Fill(trackParticle->eta());

    if(bestMatchedTTTrack.isNonnull() && bestMatchedTTTrack->getMomentum(4).perp() >= triggerAlgo->ptCut)
      ttMuonGenEtaMuons_withPtCuts->Fill(trackParticle->eta());

    if(ptOfBestL1MuCand > 0 && bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
      muCandGenEtaMuons_withPtCuts->Fill(trackParticle->eta());


      auto& layerHitBits = bestL1MuCand->getFiredLayerBits(30);
      int firedLayers = layerHitBits.count();
      int pdfSum = bestL1MuCand->pdfSum();
      muonsPdfSumFiredPlanes->Fill(pdfSum, firedLayers);

      betaLikelihoodFiredPlanesMuons->Fill(bestL1MuCand->getBetaLikelihood(), firedLayers);
      muCandBetaMuons->Fill(bestL1MuCand->getBeta()); //todo remove as it is redundant
      double l1MuBeta = bestL1MuCand->getBeta();
      if(l1MuBeta >= 1)
        l1MuBeta = 0.99;

      if(abs(trackParticle->pdgId()) == 1000015 ) //TODO only status are selected here
        betaGenBetaL1Mu->Fill(trackParticle->p4().Beta(), l1MuBeta);
    }
    else {//no candidate with pt > ptCut
      if(abs(trackParticle->pdgId()) == 1000015 ) //TODO only status are selected here
        betaGenBetaL1Mu->Fill(trackParticle->p4().Beta(), -1);
    }
  }

  int minMuPt = 3; //3 GeV
  if(ptOfBestL1MuCand > 0) {
    LogTrace("l1tMuBayesEventPrint")<<"\n"<<triggerAlgo->name<<" best muCand track: "<<toString(*bestL1MuCand)<<endl<<endl;

    double muCandPt = bestL1MuCand->getPt();
    if(trackParticle->pt() > minMuPt) {
      muCandGenEtaMuons->Fill(trackParticle->eta());
      muCandGenPhiMuons->Fill(trackParticle->phi());
    }

    muCandGenPtMuons->Fill(trackParticle->pt());

    if(trackParticle->eventId().event() == 0) {
      ptGenPtMuCandMuonsEv0->Fill(trackParticle->pt(), muCandPt);

      if( abs(trackParticle->eta() ) < 0.82)
        ptGenPtMuCandMuonsEv0Barrel->Fill(trackParticle->pt(), muCandPt);
      else if( abs(trackParticle->eta() ) < 1.24)
        ptGenPtMuCandMuonsEv0Overlap->Fill(trackParticle->pt(), muCandPt);
      else if( abs(trackParticle->eta() ) < 2.4)
        ptGenPtMuCandMuonsEv0Endcap->Fill(trackParticle->pt(), muCandPt);
    }
    else {
      ptGenPtMuCandMuonsPu->Fill(trackParticle->pt(), muCandPt);
    }
  }
  else {//no candidate
    if(trackParticle->eventId().event() == 0) {
      ptGenPtMuCandMuonsEv0->Fill(trackParticle->pt(), 0);

      if( abs(trackParticle->eta() ) < 0.82)
        ptGenPtMuCandMuonsEv0Barrel->Fill(trackParticle->pt(), 0);
      else if( abs(trackParticle->eta() ) < 1.24)
        ptGenPtMuCandMuonsEv0Overlap->Fill(trackParticle->pt(), 0);
      else if( abs(trackParticle->eta() ) < 2.4)
        ptGenPtMuCandMuonsEv0Endcap->Fill(trackParticle->pt(), 0);
    }
    else {
      ptGenPtMuCandMuonsPu->Fill(trackParticle->pt(), 0);
    }

    lostTtMuonPt->Fill(trackParticle->pt());
    if(trackParticle->pt() > 20) {
      lostTtMuonEta_ptGen20GeV->Fill(trackParticle->eta());
      lostTtMuonPhi_ptGen20GeV->Fill(trackParticle->phi());
    }

    if(trackParticle->eventId().event() == 0 &&
        trackParticle->pt() >= ptGenFrom && trackParticle->pt() <= ptGenTo) {
      /*if(notAcceptedL1MuCand == nullptr) {
        edm::LogImportant("l1tMuBayesEventPrint")<<"\nrun:lumi:event "<<event.run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<<" "<<triggerAlgo->name<<" no correlator candidate. lost high pT muon"<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<< printTrackigParticleShort(trackParticle)<<endl<<endl; //<<"Line: "<<__LINE__
      }
      else {
        edm::LogImportant("l1tMuBayesEventPrint")<<"\nrun:lumi:event "<<event.run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<<" "<<triggerAlgo->name<<" lost high pT muon due to algo cuts, correlator candidate:\n"<<toString(*notAcceptedL1MuCand)<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<< printTrackigParticleShort(trackParticle)<<endl<<endl; //<<"Line: "<<__LINE__<<". "
      }*/
    }
  }
}


class RateAnalyzer: public AnalyserBase {
public:
  RateAnalyzer(std::shared_ptr<TriggerAlgo>& triggerAlgo, edm::Service<TFileService>& fs): AnalyserBase(triggerAlgo) {
    TFileDirectory subDir = fs->mkdir( ("RateAnalyzer_" + triggerAlgo->name).c_str());

    const int ptBins = 1000;
    const int etaBins = 100;
    const int phiBins = 360;

    muCandPt = subDir.make<TH1D>("muCandPt", "muCandPt; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
    muCandEta = subDir.make<TH1D>("muCandEta", "muCandEta; eta; #events", etaBins, -2.4, 2.4);
    muCandPhi = subDir.make<TH1D>("muCandPhi", "muCandPhi; phi; #events", phiBins, -M_PI, M_PI);

    muCandPtMuons = subDir.make<TH1D>("muCandPtMuons", "muCandPtMuons; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
    //muCandEtaMuons = subDir.make<TH1D>("muCandEtaMuons", "muCandEtaMuons; eta; #events", etaBins, -2.4, 2.4);
    //muCandPhiMuons = subDir.make<TH1D>("muCandPhiMuons", "muCandPhiMuons; phi; #events", phiBins, -M_PI, M_PI);

    muCandPtFakes = subDir.make<TH1D>("muCandPtFakes", "muCandPtFakes; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
    //muCandEtaFakes = subDir.make<TH1D>("muCandEtaFakes", "muCandEtaFakes; eta; #events", etaBins, -2.4, 2.4);
    //muCandPhiFakes = subDir.make<TH1D>("muCandPhiFakes", "muCandPhiFakes; phi; #events", phiBins, -M_PI, M_PI);

    muCandPtWrongTag = subDir.make<TH1D>("muCandPtWrongTag", "muCandPtWrongTag; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
    //muCandEtaWrongTag = subDir.make<TH1D>("muCandEtaWrongTag", "muCandEtaWrongTag; eta; #events", etaBins, -2.4, 2.4);
    //muCandPhiWrongTag = subDir.make<TH1D>("muCandPhiWrongTag", "muCandPhiWrongTag; phi; #events", phiBins, -M_PI, M_PI);

    pdfSumFiredPlanesNotMuons = subDir.make<TH2I>("pdfSumFiredPlanesNotMuons", "pdfSumFiredPlanesNotMuons; pdfSum; firedPlanes; #", 100, 0, 10000, 19, -0.5, 18.5);
    betaLikelihoodFiredPlanesNotMuons = subDir.make<TH2I>("notMuonsBetaLikelihoodFiredPlanes", "notMuonsBetaLikelihoodFiredPlanes; betaLikelihood; FiredPlanes; #", 100, 0, 100, 19, -0.5, 18.5);
    muCandBetaNotMuons = subDir.make<TH1I>("muCandBetaNotMuons", "muCandBetaNotMuons; beta measured; #events", 22, 0., 1.1);

    chi2GenuineTTTrackMuCand = subDir.make<TH1D>("chi2GenuineTTTrackMuCand", "chi2GenuineTTTrackMuCand; chi2; #events", 40, 0., 200.);
    chi2FakeTTTrackMuCand = subDir.make<TH1D>("chi2FakeTTTrackMuCand", "chi2FakeTTTrackMuCand; chi2; #events", 40, 0., 200.);

    chi2DofGenuineTTTrackMuCand = subDir.make<TH1D>("chi2DofGenuineTTTrackMuCand", "chi2DofGenuineTTTrackMuCand; chi2; #events", 40, 0., 200.);
    chi2DofFakeTTTrackMuCand = subDir.make<TH1D>("chi2DofFakeTTTrackMuCand", "chi2DofFakeTTTrackMuCand; chi2; #events", 40, 0., 200.);

    //etaGenPtGenLostBx1 = subDir.make<TH2I>("etaGenPtGenLostBx1", "etaGenPtGenLostBx1; eta; pt [GeV]; ", 50, -2.4, 2.4, 50, 0, 100);
    //etaGenPtGenLostBx2 = subDir.make<TH2I>("etaGenPtGenLostBx2", "etaGenPtGenLostBx2; eta; pt [GeV]; ", 50, -2.4, 2.4, 50, 0, 100);
  }

  virtual ~RateAnalyzer() {}

  void fillHistos(const edm::Event& event, const edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& MCTruthTTTrackHandle,
      edm::Ptr< TrackingParticle > bestMuInBx1, edm::Ptr< TrackingParticle > bestMuInBx2);

private:
  TH1D* muCandPt = nullptr;
  TH1D* muCandEta = nullptr;
  TH1D* muCandPhi = nullptr;

  TH1D* muCandPtMuons = nullptr;
  //TH1D* muCandEtaMuons = nullptr;
  //TH1D* muCandPhiMuons = nullptr;

  //nominator
  TH1D* muCandPtFakes = nullptr;
  //TH1D* muCandEtaFakes = nullptr;
  //TH1D* muCandPhiFakes = nullptr;

  TH1D* muCandPtWrongTag = nullptr;
  //TH1D* muCandEtaWrongTag = nullptr;
  //TH1D* muCandPhiWrongTag = nullptr;

  TH2I* pdfSumFiredPlanesNotMuons = nullptr;
  TH2I* betaLikelihoodFiredPlanesNotMuons = nullptr;
  TH1I* muCandBetaNotMuons = nullptr;

  TH1D* chi2GenuineTTTrackMuCand = nullptr;
  TH1D* chi2FakeTTTrackMuCand = nullptr;

  TH1D* chi2DofGenuineTTTrackMuCand = nullptr;
  TH1D* chi2DofFakeTTTrackMuCand = nullptr;

  //TH2I* etaGenPtGenLostBx1 = nullptr;
  //TH2I* etaGenPtGenLostBx2 = nullptr;
};


void RateAnalyzer::fillHistos(const edm::Event& event, const edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& MCTruthTTTrackHandle,
    edm::Ptr< TrackingParticle > bestMuInBx1, edm::Ptr< TrackingParticle > bestMuInBx2) {
  if(bestL1MuCand != nullptr) {
    muCandPt->Fill(bestL1MuCand->getPt()); //filling the data from the l1track_ptr, but this is principle should be the same as in the omtfMuon
    muCandEta->Fill(bestL1MuCand->getEta());
    muCandPhi->Fill(bestL1MuCand->getPhi());

    auto& ttTrackPtr = bestL1MuCand->getTtTrackPtr();

    if(ttTrackPtr.isNull())
      return;

    edm::Ptr< TrackingParticle > tpMatchedToBestL1MuCand = MCTruthTTTrackHandle->findTrackingParticlePtr(ttTrackPtr);


    auto& layerHitBits = bestL1MuCand->getFiredLayerBits(30);
    int firedLayers = layerHitBits.count();
    int pdfSum = (int)bestL1MuCand->pdfSum();

    int L1Tk_nPar = 4; //todo take from config
    //int ttTrackNstub = ttTrackPtr->getStubRefs().size();
    float chi2dof = ttTrackPtr->getChi2Red(L1Tk_nPar); //(ttTrackPtr->getChi2(L1Tk_nPar)) / (2*ttTrackNstub - L1Tk_nPar);

    if(tpMatchedToBestL1MuCand.isNonnull() ) {
      if(abs(tpMatchedToBestL1MuCand->pdgId()) == 13) {
        muCandPtMuons->Fill(bestL1MuCand->getPt());
        if(bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
          //muCandEtaMuons->Fill(bestL1MuCand->getEta());
          //muCandPhiMuons->Fill(bestL1MuCand->getPhi());
        }
      }
      else {
        muCandPtWrongTag->Fill(bestL1MuCand->getPt());
        if(bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
          //muCandEtaWrongTag->Fill(bestL1MuCand->getEta());
          //muCandPhiWrongTag->Fill(bestL1MuCand->getPhi());

          edm::LogImportant("l1tMuBayesEventPrint")<<"\nrun:lumi:event "<<event.run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
          edm::LogImportant("l1tMuBayesEventPrint")<<" "<<triggerAlgo->name<<" wrongTag muCand track "<<(bestL1MuCand->getPt() > 18 ? " high Pt" : "")<<"\n"<<toString(*bestL1MuCand)<<endl;
          bool tmp_trk_genuine = MCTruthTTTrackHandle->isGenuine(ttTrackPtr);
          bool tmp_trk_loosegenuine = MCTruthTTTrackHandle->isLooselyGenuine(ttTrackPtr);
          edm::LogImportant("l1tMuBayesEventPrint") << printTTTRack(ttTrackPtr, tmp_trk_genuine, tmp_trk_loosegenuine, L1Tk_nPar);
          if(!tpMatchedToBestL1MuCand.isNull() ) {
            edm::LogImportant("l1tMuBayesEventPrint") <<"L:"<<__LINE__<< "TP matched to muCand:";
            edm::LogImportant("l1tMuBayesEventPrint") <<"L:"<<__LINE__<<" "<< printTrackigParticleShort(tpMatchedToBestL1MuCand, MCTruthTTTrackHandle); //<<"Line: "<<__LINE__<<". "
          }
        }
      }

      chi2GenuineTTTrackMuCand->Fill(ttTrackPtr->getChi2(L1Tk_nPar) );
      chi2DofGenuineTTTrackMuCand->Fill(chi2dof);
    }
    else {
      muCandPtFakes->Fill(bestL1MuCand->getPt());
      if(bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
        //muCandEtaFakes->Fill(bestL1MuCand->getEta());
        //muCandPhiFakes->Fill(bestL1MuCand->getPhi());

        pdfSumFiredPlanesNotMuons->Fill(pdfSum, firedLayers);
        betaLikelihoodFiredPlanesNotMuons->Fill(bestL1MuCand->getBetaLikelihood(), firedLayers);
        muCandBetaNotMuons->Fill(bestL1MuCand->getBeta());

        chi2FakeTTTrackMuCand->Fill(ttTrackPtr->getChi2(L1Tk_nPar) );
        chi2DofFakeTTTrackMuCand->Fill(chi2dof);
      }

      if(bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
        bool unknown  = MCTruthTTTrackHandle->isUnknown(ttTrackPtr);
        bool combinatoric  = MCTruthTTTrackHandle->isCombinatoric(ttTrackPtr);
        edm::LogImportant("l1tMuBayesEventPrint")<<"\nrun:lumi:event "<<event.run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<<" "<<triggerAlgo->name<<" fake muCand track "<<(bestL1MuCand->getPt() > 18 ? " high Pt" : "")<<"\n"<<toString(*bestL1MuCand)<<endl;
        edm::LogImportant("l1tMuBayesEventPrint") << " no matched track particle, the ttTrack is fake: "
            << (unknown ? " unknown " : "")
            << (combinatoric ? " combinatoric " : "")
            <<endl;
      }
    }


    //this triggerAlgo fired, because the bestL1MuCand is not null
/*    if(bestL1MuCand->getPt() >= triggerAlgo->ptCut) {
      if(bestMuInBx1.isNonnull())
        etaGenPtGenLostBx1->Fill(bestMuInBx1->eta(), bestMuInBx1->pt() );

      if(bestMuInBx2.isNonnull())
        etaGenPtGenLostBx2->Fill(bestMuInBx2->eta(), bestMuInBx2->pt() );
    }*/
  }
}


class MuCandsMatchingAnalyzer: public AnalyserBase {
public:
  class MatchingCategory {
  public:

    virtual ~MatchingCategory() {};
    std::string name;

    TH1D* candPt = nullptr;
    TH1D* candEta = nullptr;
    //TH1D* muCandPhi = nullptr;

    TH2I* candPt_vertexRho = nullptr; //rho is sqrt(x^2 + y^2)
    TH2I* pdfSumNFiredLayers = nullptr;
    TH2I* candPtFiredMuonLayers = nullptr;
    TH2I* chi2NStubs = nullptr;

    TH2I* ptGenDeltaPt = nullptr;
    TH2I* ptGenDeltaPhi = nullptr;

    MatchingCategory(std::string name, TFileDirectory& subDir, std::shared_ptr<TriggerAlgo>& triggerAlgo): name(name) {
      const int ptBins = 1000;
      const int etaBins = 100;
      const int phiBins = 360;

      candPt = subDir.make<TH1D>( ("candPt_" + name).c_str(), ("candPt_" + name + "; ttTrack pt [GeV]; #events").c_str(), ptBins, 0., 500.);

      std::string ptCutStr = " pt > " + std::to_string( (int)(triggerAlgo->ptCut) ) + " GeV ";

      candEta = subDir.make<TH1D>( ("candEta_" + name).c_str(), ("candEta " + ptCutStr + name + ";  eta; #events").c_str(), etaBins, -2.4, 2.4);
      //candPhi = subDir.make<TH1D>("candPhi", "candPhi; phi; #events", phiBins, -M_PI, M_PI);

      candPt_vertexRho = subDir.make<TH2I>( ("candPt_vertexRho_" + name).c_str(), ("candPt_vertexRho_" + name + "; ttTrack pt [GeV]; #rho=#sqrt{x^{2} + y^{2} } [cm]").c_str(), 20, 0., 100., 20, 0., 200.); //fixme is rho in cm?

      pdfSumNFiredLayers = subDir.make<TH2I>( ("pdfSumNFiredLayers_" + name).c_str(), ("pdfSumNFiredLayers "+ name + ptCutStr + "; pdfSum; NFiredLayers; #").c_str(), 100, 0, 10000, 19, -0.5, 18.5);

      candPtFiredMuonLayers = subDir.make<TH2I>( ("candPtFiredMuonLayers_" + name).c_str(), ("candPtFiredMuonLayers " + name + "; ttTrack pt [GeV]; muon layer; #").c_str(), 50, 0, 100, 23, -0.5, 22.5);

      chi2NStubs = subDir.make<TH2I>( ("chi2NStubs_" + name).c_str(), ("chi2NStubs "+ name + ptCutStr + "; chi2; #nStubs").c_str(), 30, 0., 300., 8, .5, 8.5);

      ptGenDeltaPt = subDir.make<TH2I>( ("ptGenDeltaPt_" + name).c_str(), ("ptGenDeltaPt " + name + "; gen pT [GeV]; (gen pT - ttTrack pT)/(gen pT) [GeV]; #").c_str(), 50, 0, 100, 100, -.5, .5);
      ptGenDeltaPhi = subDir.make<TH2I>( ("ptGenDeltaPhi_" + name).c_str(), ("ptGenDeltaPhi " + name +  "; gen pT [GeV]; (gen phi - ttTrack phi); #").c_str(), 50, 0, 100,  60, -0.3, 0.3);
    }

    virtual void fillHistos(const l1t::BayesMuCorrelatorTrack& l1MuCand, const edm::Ptr< TrackingParticle >& tpMatchedToL1MuCand, bool passesPtCut) {
      int L1Tk_nPar =  4; //TODO take form config
      double pt = l1MuCand.getPt();
      double eta = l1MuCand.getEta();

      candPt->Fill(pt);

      auto& layerHitBits = l1MuCand.getFiredLayerBits(30);
      if( passesPtCut) {
        candEta->Fill(eta);
        int firedLayers = layerHitBits.count();
        int pdfSum = l1MuCand.pdfSum();
        pdfSumNFiredLayers->Fill(pdfSum, firedLayers);

        chi2NStubs->Fill(l1MuCand.getTtTrackPtr()->getChi2(L1Tk_nPar), l1MuCand.getTtTrackPtr()->getStubRefs().size());
      }

      for(unsigned int layer = 0; layer < layerHitBits.size(); layer++)
        if(layerHitBits[layer])
          candPtFiredMuonLayers->Fill(pt, layer);

      if(tpMatchedToL1MuCand.isNonnull()) {
        candPt_vertexRho->Fill(pt, tpMatchedToL1MuCand->vertex().Rho());

        ptGenDeltaPt->Fill(tpMatchedToL1MuCand->pt(), (tpMatchedToL1MuCand->pt() - pt) / tpMatchedToL1MuCand->pt() );
        ptGenDeltaPhi->Fill(tpMatchedToL1MuCand->pt(), tpMatchedToL1MuCand->phi() - l1MuCand.getPhi() );
      }
    }
  };

  class MatchingCategoryDecayedToMuon: public MatchingCategory {
  public:
    virtual ~MatchingCategoryDecayedToMuon() {};

    TH2I* candPt_decayVertexRho = nullptr;

    MatchingCategoryDecayedToMuon(std::string name, TFileDirectory& subDir, std::shared_ptr<TriggerAlgo>& triggerAlgo): MatchingCategory(name, subDir, triggerAlgo) {
      candPt_decayVertexRho = subDir.make<TH2I>( ("candPt_decayVertexRho_" + name).c_str(), ("candPt_decayVertexRho_" + name + "; ttTrack pt [GeV]; #rho=#sqrt{x^{2} + y^{2} } [cm]").c_str(), 10, 0., 50., 50, 0., 500.); //fixme is rho in cm?
    }

    virtual void fillHistos(const l1t::BayesMuCorrelatorTrack& l1MuCand, const edm::Ptr< TrackingParticle >& tpMatchedToL1MuCand, bool passesPtCut) {
      MatchingCategory::fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);

      candPt_decayVertexRho->Fill(l1MuCand.getPt(), tpMatchedToL1MuCand->decayVertices()[0]->position().rho());
    }
  };


  std::unique_ptr<MatchingCategory> muons;
  std::unique_ptr<MatchingCategory> pions;
  std::unique_ptr<MatchingCategoryDecayedToMuon> pionsDecayedToMu;
  std::unique_ptr<MatchingCategory> pionsNotDecayedToMu;
  std::unique_ptr<MatchingCategory> kaons;
  std::unique_ptr<MatchingCategoryDecayedToMuon> kaonsDecayedToMu;
  std::unique_ptr<MatchingCategory> kaonsNotDecayedToMu;
  std::unique_ptr<MatchingCategory> otherParts;
  std::unique_ptr<MatchingCategory> veryLooseMuons;
  std::unique_ptr<MatchingCategory> fakes;

  TH1D* candPt = nullptr;
  TH1D* candEta = nullptr;

  MuCandsMatchingAnalyzer(std::shared_ptr<TriggerAlgo>& triggerAlgo, edm::Service<TFileService>& fs): AnalyserBase(triggerAlgo) {
    TFileDirectory subDir = fs->mkdir( ("MuCandsMatchingAnalyzer_" + triggerAlgo->name).c_str());
    const int ptBins = 1000;
    const int etaBins = 100;

    candPt = subDir.make<TH1D>( ("candPt"), ("candPt; ttTrack pt [GeV]; #events"), ptBins, 0., 500.);;
    candEta = subDir.make<TH1D>( ("candEta"), ("candEta;  eta; #events"), etaBins, -2.4, 2.4);

    muons = std::make_unique<MatchingCategory>("muons", subDir, triggerAlgo);

    pions = std::make_unique<MatchingCategory>("pions", subDir, triggerAlgo);
    pionsDecayedToMu = std::make_unique<MatchingCategoryDecayedToMuon>("pionsDecayedToMu", subDir, triggerAlgo);
    pionsNotDecayedToMu = std::make_unique<MatchingCategory>("pionsNotDecayedToMu", subDir, triggerAlgo);

    kaons = std::make_unique<MatchingCategory>("kaons", subDir, triggerAlgo);
    kaonsDecayedToMu = std::make_unique<MatchingCategoryDecayedToMuon>("kaonsDecayedToMu", subDir, triggerAlgo);
    kaonsNotDecayedToMu = std::make_unique<MatchingCategory>("kaonsNotDecayedToMu", subDir, triggerAlgo);

    otherParts = std::make_unique<MatchingCategory>("otherParts", subDir, triggerAlgo);

    veryLooseMuons = std::make_unique<MatchingCategory>("veryLooseMuons", subDir, triggerAlgo);
    fakes = std::make_unique<MatchingCategory>("fakes", subDir, triggerAlgo);
  }

  virtual ~MuCandsMatchingAnalyzer() {}

  void fillHistos(const edm::Event& event, const edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& MCTruthTTTrackHandle,
      const std::vector<edm::Ptr< TrackingParticle > >& muonTrackingParticles, const l1t::BayesMuCorrelatorTrack& l1MuCand);
};

bool decaysToMuon(edm::Ptr< TrackingParticle >& tpPtr) { //, const std::vector<edm::Ptr< TrackingParticle > >& muonTrackingParticles
  //LogTrace("l1tMuBayesEventPrint")<<std::endl;
  //LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" "<<printTrackigParticleShort(tpPtr);

  for (auto& decayVert : tpPtr->decayVertices()) {
    //LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" decayVertex r "<<decayVert->position().r()<<" phi "<<decayVert->position().phi()<<" z "<<decayVert->position().z()<<" rho "<<decayVert->position().rho()<<std::endl;

    for(auto& daughterTrack : decayVert->daughterTracks() ) {
//      LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" daughterTrack: pdgId "<<daughterTrack->pdgId()
//        << " vertex: r "<<daughterTrack->vertex().r()<<" phi "<<daughterTrack->vertex().phi()<<" x "<<daughterTrack->vertex().x()<<" y "<<daughterTrack->vertex().y()<<" z "<<daughterTrack->vertex().z()
//        <<" rho "<<daughterTrack->vertex().rho()<<" pt "<<daughterTrack->pt()<<std::endl;

      if( abs(daughterTrack->pdgId()) == 13) { //looks that neutrino is not present in the daughterTrack - since it has no track
        return true;
      }
    }

    /*    for(auto& muonTrackingParticle : muonTrackingParticles) {

    }*/
  }

  //LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" dont decays to muon"<<std::endl;
  return false;
}


void MuCandsMatchingAnalyzer::fillHistos(const edm::Event& event, const edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& MCTruthTTTrackHandle,
    const std::vector<edm::Ptr< TrackingParticle > >& muonTrackingParticles, const l1t::BayesMuCorrelatorTrack& l1MuCand) {
  if(triggerAlgo->accept(l1MuCand) == false)
    return;

  candPt->Fill(l1MuCand.getPt());
  candEta->Fill(l1MuCand.getEta());
  //muCandPhi->Fill(itL1MuCand->getPhi());

  auto& ttTrackPtr = l1MuCand.getTtTrackPtr();

  if(ttTrackPtr.isNull()) //possible only if the correlation algorithm was done with something else than ttTracks
    return;

  edm::Ptr< TrackingParticle > tpMatchedToL1MuCand = MCTruthTTTrackHandle->findTrackingParticlePtr(ttTrackPtr);

/*
  auto& layerHitBits = l1MuCand.getFiredLayerBits();
  int firedLayers = layerHitBits.count();
  int pdfSum = (int)l1MuCand.pdfSum();

  int L1Tk_nPar = 4; //todo take from config
  //int ttTrackNstub = ttTrackPtr->getStubRefs().size();
  float chi2dof = ttTrackPtr->getChi2Red(L1Tk_nPar); //(ttTrackPtr->getChi2(L1Tk_nPar)) / (2*ttTrackNstub - L1Tk_nPar);
*/

  bool passesPtCut = (l1MuCand.getPt() >= triggerAlgo->ptCut);

  if(tpMatchedToL1MuCand.isNonnull() ) {
    if(abs(tpMatchedToL1MuCand->pdgId()) == 13) {
      muons->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
    }
    else if(abs(tpMatchedToL1MuCand->pdgId()) == 211) {
      pions->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
      if(decaysToMuon(tpMatchedToL1MuCand)) {
        pionsDecayedToMu->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
      }
      else
        pionsNotDecayedToMu->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
    }
    else if(abs(tpMatchedToL1MuCand->pdgId()) == 321) {
      kaons->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
      if(decaysToMuon(tpMatchedToL1MuCand)) {
        kaonsDecayedToMu->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
      }
      else
        kaonsNotDecayedToMu->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
    }
    else {
      otherParts->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
    }

    //chi2GenuineTTTrackMuCand->Fill(ttTrackPtr->getChi2(L1Tk_nPar) );
    //chi2DofGenuineTTTrackMuCand->Fill(chi2dof);
  }
  else {
    bool isVeryLoose = false;
    for(auto& muonTrackingPart : muonTrackingParticles) {
      //here we have ttTracks tagged as muon by correlator that have no matching genuine/loose genuine tracking particle
      //so we go over all muonTrackingParticles and check if muonTrackingParticle has given ttTrack matched,
      //here, "match" means ttTracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters - so it is very loose match
      std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(muonTrackingPart);
      for(auto& matchedTTTrack : matchedTracks) {
        if(matchedTTTrack == ttTrackPtr) {
          veryLooseMuons->fillHistos(l1MuCand, muonTrackingPart, passesPtCut);
          isVeryLoose = true;

          LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" veryLooseMuon "<<printTrackigParticleShort(muonTrackingPart);
          break;
        }
      }
      if(isVeryLoose)
        break;
    }

    if(!isVeryLoose)
      fakes->fillHistos(l1MuCand, tpMatchedToL1MuCand, passesPtCut);
  }
}

//object definition
class MuCorrelatorAnalyzer : public edm::EDAnalyzer {
public:

  //constructor, function is called when new object is created
  explicit MuCorrelatorAnalyzer(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~MuCorrelatorAnalyzer();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();


private:
  //void analyzeSimTracks(edm::Handle<edm::SimTrackContainer>& simTraksHandle, edm::Handle<l1t::BayesMuCorrTrackBxCollection>& muCorrTracksHandle);
  void analyzeSimTracks(edm::Handle<edm::SimTrackContainer>& simTraksHandle, edm::Handle<reco::GenParticleCollection>& genPartHandle);

  void hscpAnalysis(edm::Handle< std::vector< TrackingParticle > >& trackingParticleHandle, edm::Handle<l1t::BayesMuCorrTrackBxCollection>& muCorrTracksHandle);
  void hscpAnalysis(const edm::Event& event, edm::Ptr< TrackingParticle >& trackPartPtr, edm::Handle<l1t::BayesMuCorrTrackBxCollection>& muCorrTracksHandle);

  std::vector<EfficiencyAnalyser> efficiencyAnalysers;
  std::vector<RateAnalyzer> rateAnalysers;
  std::vector<std::unique_ptr<MuCandsMatchingAnalyzer> > muCandsMatchingAnalyzers;

  std::unique_ptr<MuCandsMatchingAnalyzer> ttTracksMatchingAnalyzer;
  std::unique_ptr<MuCandsMatchingAnalyzer> ttTracksMatchingAnalyzerBarrel;

  std::vector<EfficiencyAnalyser> efficiencyAnalysersCorrelatorWithTrackPart; //used when the correlator works with the Tracking particles and not the ttTracks

  edm::ParameterSet parameterSet;
  unsigned int eventCount;
  //TFile* myRootFile;
  double etaCutFrom = 0;
  double etaCutTo = 2.4;

  int minMuPt = 3; //3 GeV

  MuCorrelatorConfig muCorrConfig;

  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)
  bool LooseMatch;      // use loose MC-matching instead
  int L1Tk_nPar;        // use 4 or 5 parameter track fit?
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt
  double TP_maxEta;     // save TPs with |eta| < maxEta
  double TP_maxZ0;      // save TPs with |z0| < maxZ0

  double TP_maxRho = 30; //[cm] maximum muon vertex rho - to not include in the efficiency analysis the muons from pions, which are not well matched to ttTracks. increase for not pointin muon analysis

  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)

  int muCandQualityCut = 12;

  string analysisType = "rate";

  //bool applyTriggerRules = false; //if true, for the efficiency plots the candidates that has getCanceledByTriggerRules() = true are not counted

  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag TrackingVertexInputTag;

  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > ttClusterToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > ttClusterMCTruthToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > TrackingVertexToken_;


  //tracking particle i.e. generated
  TH1I* gpPerEvent = nullptr;
  TH1I* ttTrackPerEvent = nullptr;

  TH1I* gpMuonPerEvent = nullptr;
  TH1I* ttTrackMuonPerEvent = nullptr;

  TH1D* ttTracksPt = nullptr;

  TH1I* ttTracksFakesPerEvent = nullptr;
  TH1I* ttTracksPt10FakesPerEvent = nullptr;
  TH1I* ttTracksPt20FakesPerEvent = nullptr;
  TH1D* ttTracksFakesPt = nullptr;
  TH1I* ttTracksPt10FakesEta = nullptr;

  TH1I* gpMuonPerBx = nullptr;

  //TH1D* ttTrackPerEventEta2_8 = nullptr;
  //TH1D* ttTrackMuonPerEventEta2_8 = nullptr;

  //TH1D* ttTrackInOmtfRegPerEvent = nullptr;
  //TH1D* ttTrackMuonInOmtfRegPerEvent = nullptr;

  //denominators
  TH1D* gpMuonPt = nullptr;
  TH1D* gpMuonEta = nullptr;
  TH1D* gpMuonEta_ptGen20GeV = nullptr;
  TH1D* gpMuonPhi = nullptr;

  TH1D* ttMuonPt = nullptr;
  TH1D* ttMuonEta = nullptr;
  TH1D* ttMuonEta_ptGen20GeV_ptTT18Gev = nullptr;
  TH1D* ttMuonPhi = nullptr;

  TH1D* ttMuonVeryLoosePt = nullptr;
  TH1D* ttMuonVeryLooseEta_ptGen20GeV_ptTT18Gev = nullptr;

  TH2I* ptGenPtTTMuonEv0 = nullptr;

  TH2I* ptGenPtTTMuonEvPu = nullptr;

  TH2I* ptGenDeltaPtTTMuon = nullptr;
  TH2I* ptGenDeltaPhiTTMuon = nullptr;
  TH2I* ptGenDeltaEtaTTMuon = nullptr;

  TH1D* hRefLayer = nullptr;

  TH2F* gpTrackEta_Pt = nullptr;
  TH2F* gpMuonEta_Pt = nullptr;

  TH2F* ttTrackEta_Pt = nullptr;
  TH2F* ttMuonEta_Pt = nullptr;

  //timing

  TH2I* betaGenBetaL1Mu = nullptr;
  TH2I* betaLikelihoodFiredPlanesMuons = nullptr;
  TH1I* hscpGenEta = nullptr;
  TH1I* hscpGenPt = nullptr;

  TH1I* simTracksBetaPt20 = nullptr;
  TH1I* simTracksPt = nullptr;

  TH1I* tpBetaPt20 = nullptr;

  TH1I* tpBetaHscpPt20 = nullptr;

  TH2I* etaGenPtGenBx1 = nullptr;

  TH2I* etaGenPtGenBx2 = nullptr;

  TH2I* ttTracksPerMuonTPvsPtGen = nullptr;

  edm::EDGetTokenT<l1t::BayesMuCorrTrackBxCollection> inputMuCorr;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> vertexSim;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;
};


MuCorrelatorAnalyzer::MuCorrelatorAnalyzer(const edm::ParameterSet& conf)
: parameterSet(conf), eventCount(0)
{
  inputMuCorr = consumes<l1t::BayesMuCorrTrackBxCollection>(edm::InputTag("simBayesMuCorrelatorTrackProducer", L1TMuonBayesMuCorrelatorTrackProducer::allTracksProductName)); //


  simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //TODO which is correct?
  //simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimTrackSrc"));

  vertexSim = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));

  genParticleToken = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles") );

  etaCutFrom = parameterSet.getParameter< double >("etaCutFrom");
  etaCutTo = parameterSet.getParameter< double >("etaCutTo");

  MyProcess        = parameterSet.getParameter< int >("MyProcess");
  DebugMode        = parameterSet.getParameter< bool >("DebugMode");
  SaveAllTracks    = parameterSet.getParameter< bool >("SaveAllTracks");
  SaveStubs        = parameterSet.getParameter< bool >("SaveStubs");
  LooseMatch       = parameterSet.getParameter< bool >("LooseMatch");
  L1Tk_nPar        = parameterSet.getParameter< int >("L1Tk_nPar");
  TP_minNStub      = parameterSet.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = parameterSet.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = parameterSet.getParameter< double >("TP_minPt");
  TP_maxEta        = parameterSet.getParameter< double >("TP_maxEta");
  TP_maxZ0         = parameterSet.getParameter< double >("TP_maxZ0");
  TP_maxRho        = parameterSet.getParameter< double >("TP_maxRho");
  L1TrackInputTag      = parameterSet.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub    = parameterSet.getParameter< int >("L1Tk_minNStub");

  muCandQualityCut = parameterSet.getParameter< int >("muCandQualityCut");
  analysisType = parameterSet.getParameter< string >("analysisType");

  //applyTriggerRules = parameterSet.getParameter< bool >("applyTriggerRules");

  L1StubInputTag      = parameterSet.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthClusterInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = parameterSet.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingVertexInputTag = parameterSet.getParameter<edm::InputTag>("TrackingVertexInputTag");

  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);

  ttStubToken_ = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  TrackingVertexToken_ = consumes< std::vector< TrackingVertex > >(TrackingVertexInputTag);


  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    throw;
  }

  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    throw;
  }
}


MuCorrelatorAnalyzer::~MuCorrelatorAnalyzer()
{ 
  cout <<" MuCorrelatorAnalyzer end" << endl;
}

void MuCorrelatorAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  std::shared_ptr<TriggerAlgo> singleMuAlgo = std::make_shared<SingleMuAlgo>(20);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts = std::make_shared<SingleMuAlgoSoftCuts>(20);

  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1 = std::make_shared<SingleMuAlgoSoftCuts1>(20);

  std::shared_ptr<TriggerAlgo> singleMuAlgoPtCut10 = std::make_shared<SingleMuAlgo>(10);
  std::shared_ptr<TriggerAlgo> singleMuAlgoBarrelPtCut10 = std::make_shared<SingleMuAlgoBarrel>(10);

  std::shared_ptr<TriggerAlgo> singleMuAlgoPtCut5 = std::make_shared<SingleMuAlgo>(5);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCutsPtCut5 = std::make_shared<SingleMuAlgoSoftCuts>(5);

  std::shared_ptr<TriggerAlgo> singleMuAlgoPtCut0 = std::make_shared<SingleMuAlgo>(0);

  std::shared_ptr<TriggerAlgo> singleMuAlgoPtCut3 = std::make_shared<SingleMuAlgo>(3);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCutsPtCut3 = std::make_shared<SingleMuAlgoSoftCuts>(3);

  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1PtCut0 = std::make_shared<SingleMuAlgoSoftCuts1>(0);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1PtCut3 = std::make_shared<SingleMuAlgoSoftCuts1>(3);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1PtCut5 = std::make_shared<SingleMuAlgoSoftCuts1>(5);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1PtCut10 = std::make_shared<SingleMuAlgoSoftCuts1>(10);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts1PtCut20 = std::make_shared<SingleMuAlgoSoftCuts1>(20);

  std::shared_ptr<TriggerAlgo> singleMuAlgoPdfSumSoftCuts = std::make_shared<SingleMuAlgoPdfSumSoftCuts>(20);

  std::shared_ptr<TriggerAlgo> singleMuAlgoBarrel = std::make_shared<SingleMuAlgoBarrel>(20);
  std::shared_ptr<TriggerAlgo> singleMuAlgoOverlap = std::make_shared<SingleMuAlgoOverlap>(20);
  std::shared_ptr<TriggerAlgo> singleMuAlgoEndcap = std::make_shared<SingleMuAlgoEndcap>(20);

  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCutsOverlap = std::make_shared<SingleMuAlgoSoftCutsOverlap>(20);
  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCutsOverlap1 = std::make_shared<SingleMuAlgoSoftCutsOverlap1>(20);

  //std::shared_ptr<TriggerAlgo> singleMuAlgoHardCuts = std::make_shared<SingleMuAlgoHardCuts>(20);

  std::shared_ptr<TriggerAlgo> hscpAlgo20 = std::make_shared<HscpAlgo>(20);
  std::shared_ptr<TriggerAlgo> hscpAlgo30 = std::make_shared<HscpAlgo>(30);

  std::shared_ptr<TriggerAlgo> hscpAlgoHardCuts20 = std::make_shared<HscpAlgoHardCuts>(20);
  std::shared_ptr<TriggerAlgo> hscpAlgoSoftCuts20 = std::make_shared<HscpAlgoSoftCuts>(20);
  std::shared_ptr<TriggerAlgo> hscpAlgoPdfSumCuts20 = std::make_shared<HscpAlgoPdfSumCuts>(20);



  if(analysisType == "rate") {
    rateAnalysers.emplace_back(singleMuAlgo, fs);

    rateAnalysers.emplace_back(singleMuAlgoSoftCuts, fs);
    rateAnalysers.emplace_back(singleMuAlgoSoftCuts1, fs);

    rateAnalysers.emplace_back(singleMuAlgoPdfSumSoftCuts, fs);

    //rateAnalysers.emplace_back(singleMuAlgoHardCuts, fs);

    rateAnalysers.emplace_back(singleMuAlgoBarrel, fs);
    rateAnalysers.emplace_back(singleMuAlgoOverlap, fs);
    rateAnalysers.emplace_back(singleMuAlgoEndcap, fs);
    rateAnalysers.emplace_back(singleMuAlgoSoftCutsOverlap, fs);

    rateAnalysers.emplace_back(singleMuAlgoSoftCutsOverlap1, fs);

    rateAnalysers.emplace_back(hscpAlgo20, fs);

    //rateAnalysers.emplace_back(hscpAlgo30, fs);

    //rateAnalysers.emplace_back(hscpAlgoHardCuts20, fs);
    rateAnalysers.emplace_back(hscpAlgoSoftCuts20, fs);
    rateAnalysers.emplace_back(hscpAlgoPdfSumCuts20, fs);
  }
  else if(analysisType == "efficiency") {
    efficiencyAnalysers.emplace_back(singleMuAlgo, 25, 10000, fs);

    efficiencyAnalysers.emplace_back(singleMuAlgoSoftCuts, 25, 10000, fs);
    efficiencyAnalysers.emplace_back(singleMuAlgoSoftCuts1, 25, 10000, fs);

    efficiencyAnalysers.emplace_back(singleMuAlgoPdfSumSoftCuts, 25, 10000, fs);

    //efficiencyAnalysers.emplace_back(singleMuAlgoPtCut5, 7, 10, fs);
    //efficiencyAnalysers.emplace_back(singleMuAlgoSoftCutsPtCut5, 7, 15, fs);

    efficiencyAnalysers.emplace_back(singleMuAlgoPtCut5, 7, 15, fs);
    efficiencyAnalysers.emplace_back(singleMuAlgoSoftCutsPtCut5, 7, 15, fs);
    efficiencyAnalysers.emplace_back(singleMuAlgoSoftCuts1PtCut5, 7, 15, fs);

    //efficiencyAnalysers.emplace_back(singleMuAlgoPtCut5, 7, 20, fs);
    //efficiencyAnalysers.emplace_back(singleMuAlgoSoftCutsPtCut5, 7, 20, fs);

    //efficiencyAnalysers.emplace_back(singleMuAlgoPtCut3, 5, 10, fs);
    //efficiencyAnalysers.emplace_back(singleMuAlgoSoftCutsPtCut3, 5, 10, fs);





    //efficiencyAnalysers.emplace_back(singleMuAlgoHardCuts, fs);

    efficiencyAnalysers.emplace_back(hscpAlgo20, 25, 10000, fs);

    //efficiencyAnalysers.emplace_back(hscpAlgo30, fs);

    //efficiencyAnalysers.emplace_back(hscpAlgoHardCuts20, fs);
    efficiencyAnalysers.emplace_back(hscpAlgoSoftCuts20, 25, 10000, fs);
    efficiencyAnalysers.emplace_back(hscpAlgoPdfSumCuts20, 25, 10000, fs);
  }
  else if(analysisType == "withTrackPart") {
    efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(singleMuAlgo, 20, 100000, fs);
    efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(hscpAlgo20, 20, 100000, fs);
    efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(hscpAlgo30, 20, 100000, fs);
    efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(hscpAlgoHardCuts20, 20, 100000,  fs);
    efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(hscpAlgoSoftCuts20, 20, 100000, fs);
    //efficiencyAnalysersCorrelatorWithTrackPart.emplace_back(hscpAlgoPdfSumCuts20, fs);

  }

  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoPtCut10, fs));
  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoBarrelPtCut10, fs));
  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgo, fs));
  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoSoftCuts, fs));
  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoSoftCuts1, fs));
  muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoPdfSumSoftCuts, fs));

  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoHardCuts, fs));
  //muCandsMatchingAnalyzers.emplace_back(hscpAlgo20, fs);


  std::shared_ptr<TriggerAlgo> allTTTRacks = std::make_shared<AllTTTRacks>(10);
  ttTracksMatchingAnalyzer = std::make_unique<MuCandsMatchingAnalyzer>(allTTTRacks, fs);

  std::shared_ptr<TriggerAlgo> allTTTRacksBarrel = std::make_shared<AllTTTRacksBarrel>(10);
  ttTracksMatchingAnalyzerBarrel = std::make_unique<MuCandsMatchingAnalyzer>(allTTTRacksBarrel, fs);

  //make a new Root file
  /*  string  outRootFile = parameterSet.getParameter<std::string>("outRootFile");
  myRootFile = new TFile(outRootFile.c_str(),"RECREATE");
  if(myRootFile == 0 || myRootFile-> IsZombie() ) {
    cout<<"file "<<outRootFile<<" not opened"<<endl;
    throw cms::Exception("MuCorrelatorConfig::getPatternPtRange: patternPts vector not initialized");

  cout<<"file created "<<outRootFile<<endl;*/
  /*  const int nIpt = 27;
  double lower[nIpt + 1] = {0, 4., 4.5, 5., 6., 7., 8., 10., 12., 14., 16., 18.,
      20., 22. , 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120., 140., 500.};*/
  //  for (int i = 0; i <= nIpt; i++){
  //    lower[i] = lower[i] - 0.25;
  //  }

  const int ptBins = 1000;
  const int etaBins = 100;
  const int phiBins = 360;
  //create histograms

  //TODO finish
  //TH1D* gpMuonPerEvent = new TH1D("gpMuonPt", "gpMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);
  //TH1D* ttTrackMuonPerEvent;

  gpPerEvent = fs->make<TH1I>("gpPerEvent", "generated particles per event, |#eta| < 2.4; N_tracks; #events", 100, 0., 5000.);
  ttTrackPerEvent = fs->make<TH1I>("ttTrackPerEvent", "ttTracks per event, |#eta| < 2.4; N_tracks; #events", 500, 0., 1000.);

  gpMuonPerEvent = fs->make<TH1I>("gpMuonPerEvent", ("generated muons per event, |#eta| < 2.4, pt > " + std::to_string(minMuPt) + "; N_tracks; #events").c_str(), 50, -0.5, 50.-0.5);
  ttTrackMuonPerEvent = fs->make<TH1I>("ttTrackMuonPerEvent", ("ttTracks matched to muons per event, |#eta| < 2.4, pt > " + std::to_string(minMuPt) + ";  N_tracks; #events").c_str(), 50, -0.5, 50.-0.5);

  gpMuonPerBx = fs->make<TH1I>("gpMuonPerBx", "generated muons per BX, |#eta| < 2.4; bx; #events", 21, -10.5, 10.5);

  //ttTrackPerEventEta2_8 = fs->make<TH1D>("ttTrackPerEventEta2_8", "ttTracks per event, |#eta| < 2.8; N_tracks; #events", 500, 0., 500.);
  //ttTrackMuonPerEventEta2_8 = fs->make<TH1D>("ttTrackMuonPerEventEta2_8", "ttTracks per event, |#eta| < 2.8; N_tracks; #events", 50, 0., 50.);

  //ttTrackInOmtfRegPerEvent = fs->make<TH1D>("ttTrackInOmtfRegPerEvent", "ttTracks in OMTF region per event; N_tracks; #events", 500, 0., 500.);
  //ttTrackMuonInOmtfRegPerEvent =fs->make<TH1D>("ttTrackMuonInOmtfRegPerEvent", "ttTracks matched to muons in OMTF region per event; N_tracks; #events", 50, 0., 50.);


  //denominators
  gpMuonPt = fs->make<TH1D>("gpMuonPt", "gpMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);
  gpMuonEta = fs->make<TH1D>("gpMuonEta", "gpMuonEta; eta; #events", etaBins, -2.4, 2.4);
  gpMuonEta_ptGen20GeV = fs->make<TH1D>("gpMuonEta_ptGen20GeV", "gpMuonEta_ptGen20GeV; eta; #events", etaBins, -2.4, 2.4);
  gpMuonPhi = fs->make<TH1D>("gpMuonPhi", "gpMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  ttMuonPt = fs->make<TH1D>("ttMuonPt", "ttMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);
  ttMuonEta = fs->make<TH1D>("ttMuonEta", "ttMuonEta; eta; #events", etaBins, -2.4, 2.4);
  ttMuonEta_ptGen20GeV_ptTT18Gev = fs->make<TH1D>("ttMuonEta_ptGen20GeV_ptTT18Gev", "ttMuonEta_ptGen20GeV_ptTT18GeV; eta; #events", etaBins, -2.4, 2.4);
  ttMuonPhi = fs->make<TH1D>("ttMuonPhi", "ttMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  ttMuonVeryLoosePt = fs->make<TH1D>("ttMuonVeryLoosePt", "ttMuonVeryLoosePt; gen pT [GeV]; #events", ptBins, 0., 500.);
  ttMuonVeryLooseEta_ptGen20GeV_ptTT18Gev = fs->make<TH1D>("ttMuonVeryLooseEta_ptGen20GeV_ptTT18Gev", "ttMuonVeryLooseEta_ptGen20GeV_ptTT18Gev; eta; #events", etaBins, -2.4, 2.4);

  ttTracksPt = fs->make<TH1D>("ttTracksPt", "ttTracksPt; ttTrack pT [GeV]; #events", ptBins, 0., 500.);;

  ttTracksFakesPerEvent = fs->make<TH1I>("ttTracksFakesPerEvent", ("fake ttTracks per event, |#eta| < 2.4, pt > " + std::to_string(minMuPt) + " GeV; N_tracks; #events").c_str(), 50, 0., 50.);
  ttTracksPt10FakesPerEvent = fs->make<TH1I>("ttTracksPt10FakesPerEvent", ("fake ttTracks per event, |#eta| < 2.4, pt > " + std::to_string(10) + " GeV; N_tracks; #events").c_str(), 50, 0., 50.);
  ttTracksPt20FakesPerEvent = fs->make<TH1I>("ttTracksPt20FakesPerEvent", ("fake ttTracks per event, |#eta| < 2.4, pt > " + std::to_string(20) + " GeV; N_tracks; #events").c_str(), 50, 0., 50.);
  ttTracksFakesPt = fs->make<TH1D>("ttTracksFakesPt", "ttTracksFakesPt; ttTrack pT [GeV]; #events", ptBins, 0., 500.);;
  ttTracksPt10FakesEta = fs->make<TH1I>("ttTracksPt10FakesEta", "fake ttTracks pt > 10 GeV; eta; #events", 50, -2.4, 2.4);


  ttTrackEta_Pt = fs->make<TH2F>("ttTrackEta_Pt", "ttTrackEta_Pt; ttTRack eta; ttTRack pT [GeV]; #events", etaBins/10, -2.4, 2.4, 20, 0, 40);
  ttMuonEta_Pt = fs->make<TH2F>("ttMuonEta_Pt", "ttMuonEta_Pt ttTracks matched to muon; ttTrack eta; ttTrack pT [GeV]; #events", etaBins/10, -2.4, 2.4, 20, 0, 40);

  gpTrackEta_Pt = fs->make<TH2F>("gpTrackEta_Pt", "gpTrackEta_Pt; eta; pT [GeV]; #events", etaBins/10, -2.4, 2.4, 20, 0, 40);
  gpMuonEta_Pt = fs->make<TH2F>("gpMuonEta_Pt", "gpMuonEta_Pt; eta; pT [GeV]; #events", etaBins/10, -2.4, 2.4, 20, 0, 40);

  //nominators
  ptGenPtTTMuonEv0 = fs->make<TH2I>("ptGenPtTTMuonEv0", "ptGenPtTTMuonEv0; gen pT [GeV]; ttTrack pT [GeV]; #", 100, 0, 100, 100, 0, 100);

  ptGenPtTTMuonEvPu = fs->make<TH2I>("ptGenPtTTMuonEvPu", "ptGenPtTTMuonEvPu; gen pT [GeV]; ttTrack pT [GeV]; #", 100, 0, 100, 100, 0, 100);

  ptGenDeltaPtTTMuon = fs->make<TH2I>("ptGenDeltaPtTTMuon", "ptGenDeltaPtTTMuon; gen pT [GeV]; (gen pT - ttTrack pT)/(gen pT) [GeV]; #", 100, 0, 100, 100, -.5, .5);
  ptGenDeltaPhiTTMuon = fs->make<TH2I>("ptGenDeltaPhiTTMuon", "ptGenDeltaPhiTTMuon; gen pT [GeV]; (gen phi - ttTrack phi); #", 100, 0, 100,  100, -0.5, 0.5);
  ptGenDeltaEtaTTMuon = fs->make<TH2I>("ptGenDeltaEtaTTMuon", "ptGenDeltaEtaTTMuon; gen pT [GeV]; (gen eta - ttTrack eta); #", 100, 0, 100,  100, -0.1, 0.1);

  hRefLayer = fs->make<TH1D>("refLayer", "refLayer; refLayer; #events", 10, -0.5, 10-0.5);

  betaGenBetaL1Mu = fs->make<TH2I>("betaGenBetaL1Mu", "betaGenBetaL1Mu; betaGen; Beta L1MuCand; #", 20, 0., 1., 42, -1., 1.1);
  betaLikelihoodFiredPlanesMuons = fs->make<TH2I>("betaLikelihoodFiredPlanesMuons", "betaLikelihoodFiredPlanesMuons; betaLikelihood; firedPlanes; #", 100, 0, 100, 19, -0.5, 18.5);

  hscpGenEta = fs->make<TH1I>("hscpGenEta", "hscpGenEta; eta; #events", etaBins, -3, 3);
  hscpGenPt  = fs->make<TH1I>("hscpGenPt", "hscpGenPt; ptGen [GeV]; #events", 100, 0, 1000);

  simTracksBetaPt20 = fs->make<TH1I>("simTracksBetaPt20", "simTracks beta, pt > 20 GeV; beta; #events",  20, 0., 1.);
  simTracksPt = fs->make<TH1I>("simTracksPt", "simTracksPt; ptGen [GeV]; #events", 100, 0, 1000);

  tpBetaPt20 = fs->make<TH1I>("tpBetaPt20", "tracking particle beta, pt > 20 GeV; beta; #events",  20, 0., 1.);

  tpBetaHscpPt20 = fs->make<TH1I>("tpBetaHscpPt20", "tracking particle beta, HSCP only, pt > 20 GeV; beta; #events",  20, 0., 1.);

  etaGenPtGenBx1 = fs->make<TH2I>("etaGenPtGenBx1", "etaGenPtGenBx1; eta; pt [GeV]; ", 50, -2.4, 2.4, 50, 0, 100);
  etaGenPtGenBx2 = fs->make<TH2I>("etaGenPtGenBx2", "etaGenPtGenBx2; eta; pt [GeV]; ", 50, -2.4, 2.4, 50, 0, 100);

  ttTracksPerMuonTPvsPtGen = fs->make<TH2I>("ttTracksPerMuonTPvsPtGen", "ttTracks matched to Muon Tracking Particle; ptGen [GeV]; ttTracks cout", 50, 0, 100, 5, -0.5, 4.5);

  edm::LogImportant("MuCorrelatorAnalyzer")<< "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;
}

void MuCorrelatorAnalyzer::endJob()
{
  cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;
}

/*void MuCorrelatorAnalyzer::hscpAnalysis(edm::Handle< std::vector< TrackingParticle > >& trackingParticleHandle, edm::Handle<l1t::BayesMuCorrTrackBxCollection>& muCorrTracksHandle) {
  for (unsigned int iTrackPart = 0; iTrackPart != trackingParticleHandle->size(); iTrackPart++ ) {
    edm::Ptr< TrackingParticle > trackPartPtr(trackingParticleHandle, iTrackPart);
    hscpAnalysis(trackPartPtr, muCorrTracksHandle);
  }
}*/


void MuCorrelatorAnalyzer::hscpAnalysis(const edm::Event& event, edm::Ptr< TrackingParticle >& trackPartPtr, edm::Handle<l1t::BayesMuCorrTrackBxCollection>& muCorrTracksHandle) {
  auto muCorrTracks = muCorrTracksHandle.product();
  int bxNumber = 0;

  if( trackPartPtr->pt() > 10) {
    LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" sim.type() "<<trackPartPtr->pdgId()
                      <<" Beta() "<<trackPartPtr->p4().Beta()<<std::endl;

    if (abs(trackPartPtr->pdgId()) == 1000015) {
      hscpGenEta->Fill(trackPartPtr->momentum().eta());
      hscpGenPt->Fill(trackPartPtr->pt());
    }

    for(auto& effAnalys : efficiencyAnalysersCorrelatorWithTrackPart)
      effAnalys.reset();

    bool wasL1MuCand = false;
    for(auto itL1MuCand = muCorrTracks->begin(bxNumber); itL1MuCand != muCorrTracks->end(bxNumber); ++itL1MuCand) {
      auto l1MuCandSimTrackPtr = itL1MuCand->getTrackPartPtr();

      if(l1MuCandSimTrackPtr.isNonnull() && l1MuCandSimTrackPtr == trackPartPtr) {
        if (abs(trackPartPtr->pdgId()) == 1000015) {
          for(auto& effAnalys : efficiencyAnalysersCorrelatorWithTrackPart)
            effAnalys.takeCanidate(itL1MuCand);
        }

        int quality = itL1MuCand->hwQual();
        if(quality < muCandQualityCut)
          continue;

        float l1MuBeta = itL1MuCand->getBeta();
        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<" l1MuBeta "<<l1MuBeta<<" "<<(l1MuBeta ? "" : "!!!!!!!")<<endl;
        if(l1MuBeta == 1)
          l1MuBeta = 0.99;
        betaGenBetaL1Mu->Fill(l1MuCandSimTrackPtr->p4().Beta(), l1MuBeta);

        betaLikelihoodFiredPlanesMuons->Fill(itL1MuCand->getBetaLikelihood(), itL1MuCand->getFiredLayerCnt());

        wasL1MuCand = true; //by definition there can be only candidate for a given tracking particle
        break;
      }
    }

    if(!wasL1MuCand) {
      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" no L1MuCand !!!!!!!!!!!!!"<<std::endl;
      betaGenBetaL1Mu->Fill(trackPartPtr->p4().Beta(), -1);
    }

    //TODO fix
    if (abs(trackPartPtr->pdgId()) == 1000015) {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > bestMatchedTTTrack;
      for(auto& effAnalys : efficiencyAnalysersCorrelatorWithTrackPart)
        effAnalys.fillHistos(event, trackPartPtr, bestMatchedTTTrack);
    }
  }
}

void MuCorrelatorAnalyzer::analyzeSimTracks(edm::Handle<edm::SimTrackContainer>& simTraksHandle, edm::Handle<reco::GenParticleCollection>& genPartHandle ) {
  LogTrace("l1tMuBayesEventPrint") <<"analyzeSimTracks\n printing simTraks"<<std::endl;
  for (unsigned int iSimTrack = 0; iSimTrack != simTraksHandle->size(); iSimTrack++ ) {
    edm::Ptr< SimTrack > simTrackPtr(simTraksHandle, iSimTrack);
    simTracksPt->Fill(simTrackPtr->momentum().pt());
    if(simTrackPtr->momentum().pt() > 20 ) {
      simTracksBetaPt20->Fill(simTrackPtr->momentum().Beta());
    }

    //if(abs(simTrackPtr->type() ) >= 1000000 && abs(simTrackPtr->type()) < 9999999)
    if(simTrackPtr->momentum().pt() > 20 )
    { //stau is 1000015
    //if(abs(simTrackPtr->type()) == 13 ) {
      LogTrace("l1tMuBayesEventPrint") <<" analyzeSimTracks: event "<<simTrackPtr->eventId().event()<<" pdgId "<<simTrackPtr->type()
              <<" pt "<<simTrackPtr->momentum().pt()<<" Beta "<<simTrackPtr->momentum().Beta()<<" genpartIndex "<<simTrackPtr->genpartIndex()<<std::endl;
    }

  }


  LogTrace("l1tMuBayesEventPrint") <<"\n printing GenParticle"<<std::endl;
  for (unsigned int iGenTrack = 0; iGenTrack != genPartHandle->size(); iGenTrack++ ) {
    edm::Ptr< reco::GenParticle > genPartPtr(genPartHandle, iGenTrack);

    if(abs(genPartPtr->pdgId()) == 1000015) {
      LogTrace("l1tMuBayesEventPrint") <<" analyzeSimTracks "<<"genPartPtr->pdgId() "<<genPartPtr->pdgId()<<" status() "<<genPartPtr->status()<<" charge() "<<genPartPtr->charge()
                      <<" p "<<genPartPtr->p()<<" pt "<<genPartPtr->pt()<<" vertex().Rho() "<<genPartPtr->vertex().Rho()<<" eta "<<genPartPtr->eta()<<" phi "<<genPartPtr->phi()
                      <<" numberOfDaughters() "<<genPartPtr->numberOfDaughters()<<" numberOfMothers "<<genPartPtr->numberOfMothers()<<endl;

      LogTrace("l1tMuBayesEventPrint")<<"mothers" <<endl;
      for(auto& mother : genPartPtr->motherRefVector()) {
        LogTrace("l1tMuBayesEventPrint") <<"      mother->pdgId() "<<mother->pdgId()
                            <<" status() "<<mother->status()
                            <<" p "<<mother->p()<<" pt "<<mother->pt()<<" vertex().Rho() "<<mother->vertex().Rho()
                            <<endl;
      }

      LogTrace("l1tMuBayesEventPrint")<<"daughters" <<endl;
      for(auto& daughter : genPartPtr->daughterRefVector()) {
        LogTrace("l1tMuBayesEventPrint") <<"      daughter->pdgId() "<<daughter->pdgId()
                              <<" status() "<<daughter->status() //<<"  numberOfMothers() "<<daughter->numberOfMothers()<<
                              <<" p "<<daughter->p()<<" pt "<<daughter->pt()<<" vertex().Rho() "<<daughter->vertex().Rho()<<endl;


        for(auto& mother : daughter->motherRefVector()) {
          LogTrace("l1tMuBayesEventPrint") <<"                mother->pdgId() "<<mother->pdgId()
                                            <<" status() "<<mother->status()
                                            <<" p "<<mother->p()<<" pt "<<mother->pt()<<" vertex().Rho() "<<mother->vertex().Rho()
                                            <<endl;
        }

      }
      LogTrace("l1tMuBayesEventPrint") <<endl;
    }
    /*else {
        if(genPartPtr->numberOfMothers() > 0) {
          bool hasStauMother = false;
          for(auto& mother : genPartPtr->motherRefVector()) {
            if( (mother->pdgId()) == 1000015 ) {
              hasStauMother = true;
            }
          }
          if(hasStauMother) {
            LogTrace("l1tMuBayesEventPrint") <<" analyzeSimTracks genPartPtr->pdgId() "<<genPartPtr->pdgId()<<" status() "<<genPartPtr->status()
                  <<" p "<<genPartPtr->p()<<" pt "<<genPartPtr->pt()<<" vertex().Rho() "<<genPartPtr->vertex().Rho()
                  <<" numberOfDaughters() "<<genPartPtr->numberOfDaughters()<<" numberOfMothers "<<genPartPtr->numberOfMothers()<<endl;

            for(auto& mother : genPartPtr->motherRefVector()) {
              LogTrace("l1tMuBayesEventPrint") <<"      mother->pdgId() "<<mother->pdgId()<<endl;
            }
          }
        }
      }*/
  }

}

void analyseMuonTrackingParticles(edm::Ptr< TrackingParticle >& tpPtr, const std::vector<SimTrack>& simTraks, const std::vector<SimVertex>& simVertexes) {
  //tpPtr->decayVertices();
  //tpPtr->parentVertex();
  //tpPtr->vertex();

  //if(tpPtr->eventId().event() != 0 )
  {//&& tpPtr->decayVertices().size() > 0
    //if ( (abs(tpPtr->pdgId()) == 211  && tpPtr->pt() > 4 && tpPtr->vertex().r() < 10) || (abs(tpPtr->pdgId()) == 13 && tpPtr->pt() > 3) )
    { //13
      //if(tpPtr->pt() > 4)
      {
        //if(tpPtr->vertex().Rho() > 1)
        {
          LogTrace("l1tMuBayesEventPrint") <<" analyseMuonTrackingPartivles: event "<<tpPtr->eventId().event()<<" pdgId "<<tpPtr->pdgId()
                << " vertex: r "<<tpPtr->vertex().r()<<" phi "<<tpPtr->vertex().phi()<<" z "<<tpPtr->vertex().z()<<" pt "<<tpPtr->pt()<<std::endl;
          for(auto& simTrack : tpPtr->g4Tracks()) {
            LogTrace("l1tMuBayesEventPrint") <<" simTrack.type() "<<simTrack.type()<<" vertIndex "<<simTrack.vertIndex()<<" decayVertices.size "<<tpPtr->decayVertices().size()<<std::endl;
            //LogTrace("l1tMuBayesEventPrint") <<" simTrack.vertIndex()).parentIndex() "<<simVertexes.at(simTrack.vertIndex()).parentIndex()<<std::endl; does not work, looks that vertices for the pileup events are not stored
          }

          for(auto& genTrack : tpPtr->genParticles() ) {
            LogTrace("l1tMuBayesEventPrint") <<" genTrack->pdgId() "<<genTrack->pdgId()<<std::endl;
          }

          for (auto& decayVert : tpPtr->decayVertices()) {
            LogTrace("l1tMuBayesEventPrint") <<" decayVertex r "<<decayVert->position().r()<<" phi "<<decayVert->position().phi()<<" z "<<decayVert->position().z()<<std::endl;
          }
          LogTrace("l1tMuBayesEventPrint") <<std::endl;
        }
      }
    }
  }
}

void MuCorrelatorAnalyzer::analyze(
    const edm::Event& event, const edm::EventSetup& es)
{
  //const unsigned int omtflayersCnt = muCorrConfig.nLayers();
  LogTrace("l1tMuBayesEventPrint") << "\n\nMuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" new event "<<event.id()<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< endl;


  edm::Handle<edm::SimTrackContainer> simTraksHandle;
  event.getByToken(simTrackToken, simTraksHandle);

  const std::vector<SimTrack>& simTraks = *(simTraksHandle.product());

  edm::Handle<edm::SimVertexContainer> simVx;
  event.getByToken(vertexSim, simVx);
  const std::vector<SimVertex>& simVertexes = *(simVx.product());

  edm::Handle<reco::GenParticleCollection> genPartHandle;
  event.getByToken(genParticleToken, genPartHandle);

  LogTrace("l1tMuBayesEventPrint") << "\nMuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" vertexSim "<<simVx->size()<< endl;

  //analyzeSimTracks(simTraksHandle, genPartHandle); //muCorrTracksHandle


  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  event.getByToken(ttTrackToken_, TTTrackHandle);

  // MC truth association maps

  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  event.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > trackingParticleHandle;
  //edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  event.getByToken(TrackingParticleToken_, trackingParticleHandle);
  //event.getByToken(TrackingVertexToken_, TrackingVertexHandle);

  //std::cout <<" L1 MUONS: "<<std::endl;
  //edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  //event.getByToken(inputOMTF, l1Omtf);
  //cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" omtfCands->size(bxNumber) "<<omtfCands->size()<< endl;

  edm::Handle<l1t::BayesMuCorrTrackBxCollection> muCorrTracksHandle;
  event.getByToken(inputMuCorr, muCorrTracksHandle);
  auto muCorrTracks = muCorrTracksHandle.product();

  int bxNumber = 0;

  //hscpAnalysis(trackingParticleHandle, muCorrTracksHandle);

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------

  //LogTrace("l1tMuBayesEventPrint") << endl << "Loop over tracking particles!" << endl;

  gpPerEvent->Fill(trackingParticleHandle->size() );
  ttTrackPerEvent->Fill(TTTrackHandle->size());

  int allGpMuCnt = 0;
  int ttTrackMuCnt = 0;
  //int ttTrackInOmtfCnt = 0;
  int ttTrackMuInOmtfCnt = 0;

  //boost::dynamic_bitset omtfCandsIsMatchedToMuonTrakPart(omtfCands->size(bxNumber));
  LogTrace("l1tMuBayesEventPrint")<<"trackingParticleHandle->size() "<<trackingParticleHandle->size()<<" isValid "<<trackingParticleHandle.isValid();

  auto printTrackigParticleShort = [&](edm::Ptr< TrackingParticle >& tpPtr) {
    std::ostringstream ostr;
    ostr<< "Tracking particle,    pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
                << " eta: " <<setw(7)<< tpPtr->eta()
                << " phi: " <<setw(7)<< tpPtr->phi()
                << " beta: " <<setw(7)<< tpPtr->p4().Beta()
                << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
                << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr).size()
                <<" bx: "<<tpPtr->eventId().bunchCrossing()
                <<dec<<" key "<<tpPtr.key();
    return ostr.str();
  };

  auto printTTTRack = [&](const edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > >& ttTrackPtr, bool genuine, bool loosegenuine) {
    std::ostringstream ostr;
    ostr << "matched ttTrack       pt: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).perp()
                     <<" charge: "<<setw(2)<<(ttTrackPtr->getRInv() > 0 ? 1 : -1)
                     << " eta: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).eta()
                     << " phi: " <<setw(7)<< ttTrackPtr->getMomentum(L1Tk_nPar).phi()
                     << " chi2: " <<setw(7)<< ttTrackPtr->getChi2(L1Tk_nPar)
                     << " consistency: " <<setw(7)<< ttTrackPtr->getStubPtConsistency(L1Tk_nPar)
                     << " z0: " <<setw(7)<< ttTrackPtr->getPOCA(L1Tk_nPar).z()
                     << " nstub: " <<setw(7)<< ttTrackPtr->getStubRefs().size()
                     << (genuine  ? " genuine " : "")
                     << (loosegenuine ? " loose genuine " : "");
    return ostr.str();
  };


  //for pre-fireing analysis
  edm::Ptr< TrackingParticle > bestMuInBx1;
  edm::Ptr< TrackingParticle > bestMuInBx2;

  std::vector<edm::Ptr< TrackingParticle > > muonTrackingParticles;

  for (unsigned int iTP = 0; iTP < trackingParticleHandle->size(); ++iTP) {
    edm::Ptr< TrackingParticle > tpPtr(trackingParticleHandle, iTP);

    //analyseMuonTrackingParticles(tpPtr, simTraks, simVertexes);

    //int tmp_eventid = tp_ptr->eventId().event();
    //if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tp_d0_prod = -tpPtr->vx()*sin(tpPtr->phi()) + tpPtr->vy()*cos(tpPtr->phi());


    gpTrackEta_Pt->Fill(tpPtr->eta(), tpPtr->pt());

    if (abs(tpPtr->pdgId()) == 13 || abs(tpPtr->pdgId()) == 1000015 ) {  //|| tpPtr->pt() > 20 //todo 1000015 is stau
      //only muons
    }
    else
      continue;

    if (tpPtr->pt() < TP_minPt)
      continue;
    if (fabs(tpPtr->eta()) > TP_maxEta)
      continue;

    if(abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) { //TODO it duplicates the above condition,
    }
    else
      continue;

    if(tpPtr->vertex().rho() > TP_maxRho)
      continue;

    if(tpPtr->pt() > 20) {
      tpBetaPt20->Fill(tpPtr->p4().Beta());

      if(abs(tpPtr->pdgId()) == 1000015 && tpPtr->eventId().event() == 0) {
        tpBetaHscpPt20->Fill(tpPtr->p4().Beta());
      }
    }

    if(analysisType == "withTrackPart")
      hscpAnalysis(event, tpPtr, muCorrTracksHandle);

    if (abs(tpPtr->pdgId()) == 1000015) {
      hscpGenEta->Fill(tpPtr->momentum().eta());
      hscpGenPt->Fill(tpPtr->pt());
    }

    if(tpPtr->eventId().bunchCrossing() == 1) {
      if(bestMuInBx1.isNull() || tpPtr->pt() > bestMuInBx1->pt())
        bestMuInBx1 = tpPtr;
    }
    else if(tpPtr->eventId().bunchCrossing() == 2) {
      if(bestMuInBx2.isNull() || tpPtr->pt() > bestMuInBx2->pt())
        bestMuInBx2 = tpPtr;
    }

    gpMuonPerBx->Fill(tpPtr->eventId().bunchCrossing());

    if(tpPtr->eventId().bunchCrossing() != 0)
      continue;

    muonTrackingParticles.emplace_back(tpPtr);

    // ----------------------------------------------------------------------------------------------
    // get d0/z0 propagated back to the IP
    float tmp_tp_t = tan(2.0*atan(1.0) - 2.0*atan(exp(-tpPtr->eta())));

    float delx = -tpPtr->vx();
    float dely = -tpPtr->vy();

    float A = 0.01*0.5696;
    float Kmagnitude = A / tpPtr->pt();

    float tmp_tp_charge = tpPtr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;

    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tpPtr->phi()));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tpPtr->phi()));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));

    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tpPtr->phi()-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tpPtr->vz() + tmp_tp_t*delphi/(2.0*K);
    // ----------------------------------------------------------------------------------------------

    //if (fabs(tmp_tp_z0) > TP_maxZ0) todo does it matter for muons??????????????????
    //  continue;

    /*    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tpPtr->vx()*tpPtr->vx() + tpPtr->vy()*tpPtr->vy());
    float tmp_tp_dxy = dxy;
    if (MyProcess==6 && (dxy > 1.0)) continue;*/


    if(tpPtr->pt() > minMuPt) {
      allGpMuCnt++;
      gpMuonEta->Fill(tpPtr->eta());
    }
    if(tpPtr->pt() > 20) {
      gpMuonEta_ptGen20GeV->Fill(tpPtr->eta());
    }

    gpMuonPt->Fill(tpPtr->pt());
    gpMuonPhi->Fill(tpPtr->phi());

    gpMuonEta_Pt->Fill(tpPtr->eta(), tpPtr->pt());

    auto printTrackigParticle = [&]() {
      std::ostringstream ostr;
      ostr<< "Tracking particle     pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
                  << " eta: " <<setw(7)<< tpPtr->eta()
                  << " phi: " <<setw(7)<< tpPtr->phi()
                  << " beta: " <<setw(7)<< tpPtr->p4().Beta()
                  << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
                  << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr).size()
                  << " bx: "<<tpPtr->eventId().bunchCrossing()
                  << dec<<" key "<<tpPtr.key()
                  << " z0: " <<setw(7)<< tmp_tp_z0 << " d0: " <<setw(7)<< tmp_tp_d0
                  << " z_prod: " <<setw(7)<< tpPtr->vz() << " d_prod: " <<setw(7)<< tp_d0_prod;
      return ostr.str();
    };

    LogTrace("l1tMuBayesEventPrint")<<"\n\n"<<"L:"<<__LINE__<<" "<<printTrackigParticle();

    // ----------------------------------------------------------------------------------------------
    // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)
    /*

    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() < 1) {
      if (DebugMode) cout << "No matching TTClusters for TP, continuing..." << endl;
      continue;
    }


    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int) theStubRefs.size();

    if (TP_minNStub > 0) {
      if (DebugMode) cout << "Only consider TPs with >= " << TP_minNStub << " stubs" << endl;
      if (nStubTP < TP_minNStub) {
        if (DebugMode) cout << "TP fails minimum nbr stubs requirement! Continuing..." << endl;
        continue;
      }
    }
     */


    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle
    std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr);

    int nMatch = 0;
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > bestMatchedTTTrack; //null at the beginning
    float chi2dofOfBestMatchedTrack = 99999;

    //we reset here, because we are looking for one candidate (highest pt) for each tracking particle
    for(auto& effAnalys : efficiencyAnalysers)
      effAnalys.reset();

    // ----------------------------------------------------------------------------------------------
    // loop over ttTracks matched to this tracking particles
    // here, "match" means ttTracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack
    /* MCTruthTTTrackHandle->findTTTrackPtrs gives for every tracking particle every ttTrack that has at least one cluster generated by this tracking particle
     * so it is very loose
     * see TTTrackAssociator< Ref_Phase2TrackerDigi_ >::produce
     */
    for(auto& matchedTTTrack : matchedTracks) {
      bool tmp_trk_genuine = MCTruthTTTrackHandle->isGenuine(matchedTTTrack);
      bool tmp_trk_loosegenuine = MCTruthTTTrackHandle->isLooselyGenuine(matchedTTTrack);

      /**
       * if the ttTrack is loosegenuine or genuine it means that there is only one tracking particle matched to it
       * i.e. the ttTrack stubs were generated by only one tracking particle
       * additionally all stubs but one were generated by this tracking particle
       * see TTTrackAssociator.cc line 142
       * the question is what exactly is this one stub that is not associated to any tracking particle? (probably unknown? or also combinatoric?)
       * see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTStub
       *
       * but if it is genuine then ALL stubs of the ttTrack must be generated by one tracking particle, see MCTruthTTTrackHandle->isGenuine
       */

/*      int matchingQuality = 0;
      if(tmp_trk_genuine)
        matchingQuality = 2;
      else if(tmp_trk_loosegenuine)
        matchingQuality = 1;
      else
        matchingQuality = 0;*/


      if (LooseMatch && !tmp_trk_loosegenuine)
        continue;
      if (!LooseMatch && !tmp_trk_genuine)
        continue;


      edm::Ptr< TrackingParticle > tpOfMatchedTTTrack = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTTTrack);
      LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTTTRack(matchedTTTrack, tmp_trk_genuine, tmp_trk_loosegenuine);

      // ----------------------------------------------------------------------------------------------
      // number of stubs in this matchedTTTrack
      int ttTrackNstub = matchedTTTrack->getStubRefs().size();

      if (ttTrackNstub < L1Tk_minNStub) {
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << "t tTrackNstub < "<<L1Tk_minNStub<<" - skipping this ttTrack" << endl;
        continue;
      }

      //dmatch_pt  = fabs(my_tp->p4().pt() - tpPtr->pt());
      /*float dmatch_eta = fabs(tpOfMatchedTTTrack->p4().eta() - tpPtr->eta());
      float dmatch_phi = fabs(tpOfMatchedTTTrack->p4().phi() - tpPtr->phi());
      float dmatch_id = tpOfMatchedTTTrack->pdgId();*/

      float tmp_trk_chi2dof = (matchedTTTrack->getChi2(L1Tk_nPar)) / (2*ttTrackNstub - L1Tk_nPar);

      // ensure that track is uniquely matched to the TP we are looking at!
      if (tpOfMatchedTTTrack.isNull()) { //if the ttTrack is at least loose genuine this cannot happened
        edm::LogImportant("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << " ttTrack matched to tracking particle is NOT matched to any tracking particle  (i.e. is not genuine nor loose) - not possible!!!!!!!!!!! " << endl;
      }
      else if(tpPtr != tpOfMatchedTTTrack) {//is this possible at all? - yes if the ttTrack was matched to another particle, but with the current has one cluster in common (possible if loose genuine)
        LogTrace("l1tMuBayesEventPrint")<<std::endl;
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" tpPtr != tpOfMatchedTTTrack !!!!!!!!!!!!";
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" current tracking particle";
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTrackigParticleShort(tpPtr);

        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTTTRack(matchedTTTrack, tmp_trk_genuine, tmp_trk_loosegenuine);
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" tpOfMatchedTTTrack";
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTrackigParticleShort(tpOfMatchedTTTrack);
        LogTrace("l1tMuBayesEventPrint")<<std::endl;
      }
      else {//the ttTrack is uniquely  matched to the current tracking particle
        //so we assume that the matched ttTrack (and thus omtf track) is at least loose genuine
        nMatch++;

        if(bestMatchedTTTrack.isNonnull()) {
          edm::LogImportant("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << "WARNING *** 2 or more (loose) genuine ttTrack match to tracking particle !!!!!\n"<<printTrackigParticleShort(tpPtr) << endl;
          edm::LogImportant("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << " current bestMatchedTTTrack\n"<<printTTTRack(bestMatchedTTTrack, false, false);
          edm::LogImportant("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << " matchedTTTrack\n"<<printTTTRack(matchedTTTrack, tmp_trk_genuine, tmp_trk_loosegenuine)<<"\n";
        }

        if (bestMatchedTTTrack.isNull() || tmp_trk_chi2dof < chi2dofOfBestMatchedTrack) {
          bestMatchedTTTrack = matchedTTTrack;
          chi2dofOfBestMatchedTrack = tmp_trk_chi2dof;
        }

        //finding omtf candidate corresponding to the matchedTTTrack
        //LogTrace("l1tMuBayesEventPrint")<<" omtfCands->size "<<omtfCands->size(bxNumber)<<endl;

        //we should find one candidate with max pt and use it for filling the histograms
        for(auto itL1MuCand = muCorrTracks->begin(bxNumber); itL1MuCand != muCorrTracks->end(bxNumber); ++itL1MuCand) {
          auto& muCandTtTrackPtr = itL1MuCand->getTtTrackPtr();
          //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > omtfTtTrackPtr(TTTrackHandle, omtfTTTrackIndex);

          if(muCandTtTrackPtr.isNonnull() && muCandTtTrackPtr.get() == matchedTTTrack.get() ) { //it should be only one omtfTrack for given ttTRack
            LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<toString(*itL1MuCand)<<endl;

            for(auto& effAnalys : efficiencyAnalysers)
              effAnalys.takeCanidate(itL1MuCand);

            break;
          }
          //in principle for a given tracking particle there can be two (or more) matched ttTracks, but then most probably the OMTF will select only one, the second will be ghost busted
        }
      }
    }// end loop over matched ttTracks

    // ----------------------------------------------------------------------------------------------
    ttTracksPerMuonTPvsPtGen->Fill(tpPtr->pt(), nMatch);

    for(auto& effAnalys : efficiencyAnalysers)
      effAnalys.fillHistos(event, tpPtr, bestMatchedTTTrack);

    if (nMatch > 1) {
      //edm::LogImportant("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << "WARNING *** 2 or more (loose) genuine ttTrack match to tracking particle !!!!! "<<printTrackigParticleShort(tpPtr) << endl;
    }

    if (nMatch > 0) {//there is ttTrack matching to the trackingParticle
      float matchTTTrackPt   = bestMatchedTTTrack->getMomentum(L1Tk_nPar).perp();
      float matchTTTrackEta  = bestMatchedTTTrack->getMomentum(L1Tk_nPar).eta();
      float matchTTTrackPhi  = bestMatchedTTTrack->getMomentum(L1Tk_nPar).phi();

      if(tpPtr->pt() > minMuPt) {
        ttTrackMuCnt++;
        ttMuonEta->Fill(tpPtr->eta());
        ttTrackMuInOmtfCnt++;
      }
      if(tpPtr->pt() > 20 && matchTTTrackPt >= 18) {
        ttMuonEta_ptGen20GeV_ptTT18Gev->Fill(tpPtr->eta());
      }

      ttMuonPt->Fill(tpPtr->pt());
      ttMuonPhi->Fill(tpPtr->phi());

      ttMuonEta_Pt->Fill(matchTTTrackEta, matchTTTrackPt);

      if(tpPtr->eventId().event() == 0) {
        ptGenPtTTMuonEv0->Fill(tpPtr->pt(), matchTTTrackPt);
      }
      else {
        ptGenPtTTMuonEvPu->Fill(tpPtr->pt(), matchTTTrackPt);
      }

      ptGenDeltaPtTTMuon->Fill(tpPtr->pt(), (tpPtr->pt() - matchTTTrackPt) / tpPtr->pt());
      ptGenDeltaPhiTTMuon->Fill(tpPtr->pt(), tpPtr->phi() - matchTTTrackPhi);
      ptGenDeltaEtaTTMuon->Fill(tpPtr->pt(), tpPtr->eta() - matchTTTrackEta);

    }
    else { //no matched ttTrack
      if(tpPtr->eventId().event() == 0) {
        ptGenPtTTMuonEv0->Fill(tpPtr->pt(), 0);
      }
      else {
        ptGenPtTTMuonEvPu->Fill(tpPtr->pt(), 0);
      }

      double tpMinPt = 5;
      if (tpPtr->pt() > tpMinPt) {
        LogTrace("l1tMuBayesEventPrint") <<"\n"<<"L:"<<__LINE__<<" no ttTrack matched to tracking particle:";
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTrackigParticle();
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" Printing all ttTracks matching this tracking particle" << endl;

        for(auto& matchedTTTrack : matchedTracks) {
          LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTTTRack(matchedTTTrack, false, false) << endl;
          if( abs(matchedTTTrack->getMomentum(L1Tk_nPar).phi() - tpPtr->phi()) < 0.01 &&
              abs(matchedTTTrack->getMomentum(L1Tk_nPar).eta() - tpPtr->eta()) < 0.01 &&
              abs( (matchedTTTrack->getMomentum(L1Tk_nPar).perp() - tpPtr->pt() ) / tpPtr->pt()) < 0.1)
          {
            if(MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTTTrack).isNull() ) { //we require that this ttTrack is not well matched to other particle
              ttMuonVeryLoosePt->Fill(tpPtr->pt());
              if(tpPtr->pt() > 20 && matchedTTTrack->getMomentum(L1Tk_nPar).perp() >= 18 ) {
                ttMuonVeryLooseEta_ptGen20GeV_ptTT18Gev->Fill(tpPtr->eta());
              }
              LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" selected as VeryLoose"<< endl;
            }
          }
        }
/*
        LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" Printing all ttTracks around this tracking particle" << endl;

        //looking for ttTrack not matched by the MCTruthTTTrackHandle
        int l1TrackIndx = 0;
        for (auto iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
          edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, l1TrackIndx++);

          float ttTrkPt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
          float ttTrkEta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
          float ttTrkPhi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
          float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm

          float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
          int tmp_trk_nstub  = (int) iterL1Track->getStubRefs().size();

          //cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;

          int tmp_trk_genuine = 0;
          int tmp_trk_loose = 0;
          int tmp_trk_unknown = 0;
          int tmp_trk_combinatoric = 0;
          if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
          if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
          if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
          if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;

          //TODO for phi  around -pi and +pi the problem is with the delta, ignoring for the moment
          if(tpPtr->pt() > tpMinPt && abs(ttTrkPhi -  tpPtr->phi()) < 0.2 && abs(ttTrkEta - tpPtr->eta()) < 0.05) {
            edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

            //if (DebugMode )
            {//&& abs(my_tp->pdgId()) == 13
              LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__ << " L1 track, pt: " << ttTrkPt //<<" RInv "<<l1track_ptr->getRInv()
                  << " eta: " << ttTrkEta << " phi: " << ttTrkPhi
                  << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub<<", "
                  <<(tmp_trk_genuine ? "genuine, " : "")
                  <<(tmp_trk_loose ? "loose, " : "")
                  <<(tmp_trk_unknown ? "unknown, " : "")
                  <<(tmp_trk_combinatoric ? "combinatoric, " : "");

              if( !my_tp.isNull() ) {
                LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__ <<" is matched to tracking particle: ";
                LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<< " Tracking particle,    pt: " <<setw(7)<< my_tp->pt() << " eta: " <<setw(7)<< my_tp->eta() << " phi: " <<setw(7)<< my_tp->phi()
                                      << " pdgid: " <<setw(7)<< my_tp->pdgId() << " eventID: " <<setw(7)<< my_tp->eventId().event()
                                      << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(my_tp).size()
                                      <<" bx: "<<my_tp->eventId().bunchCrossing()<<" address "<<my_tp.get()
                                      <<" key "<<my_tp.key()
                                      <<" "<<(my_tp == tpPtr ? " my_tp == tp_ptr " : "my_tp != tp_ptr")<< endl;
              }
              else {
                LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" no matching tracking particle"<<endl;
              }
              LogTrace("l1tMuBayesEventPrint")<<"";
            }
          }
        }
        */
      }
    }
  } //end loop tracking particles

  gpMuonPerEvent->Fill(allGpMuCnt);
  ttTrackMuonPerEvent->Fill(ttTrackMuCnt);
  //ttTrackMuonInOmtfRegPerEvent->Fill(ttTrackMuInOmtfCnt);

  /*
   * in principle only running with LooseMatch= True has sense, otherwise how to count the ttTracks that are loose genuine matched to muon tracking particle
   * and are tagged as muon by omtf?
   */

  for(auto& rateAnalyser : rateAnalysers)
    rateAnalyser.reset();

  for(auto itL1MuCand = muCorrTracks->begin(bxNumber);
        itL1MuCand != muCorrTracks->end(bxNumber); ++itL1MuCand) {

    //int omtfTTTrackIndex = (int)itL1MuCand->trackAddress().at(3);
    //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ttTrackPtr(TTTrackHandle, omtfTTTrackIndex);

/*    auto& ttTrackPtr = itL1MuCand->getTtTrackPtr();
    if(ttTrackPtr.isNull() ) {
      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<"ttTrackPtr.isNull() !!!!!!!!!!!!!"<<endl;
      //continue;
    }
    else {
      int ttTrackNstub = ttTrackPtr->getStubRefs().size();
      if (ttTrackNstub < L1Tk_minNStub) this cut is applied in the TTTracksInputMaker::loadTTTracks
        continue;
    }*/

    if(abs(itL1MuCand->getEta()) >= etaCutFrom && abs(itL1MuCand->getEta() ) <= etaCutTo) {
    }
    else
      continue;

    for(auto& rateAnalyser : rateAnalysers)
      rateAnalyser.takeCanidate(itL1MuCand);

    for(auto& muCandsMatchingAnalyzer : muCandsMatchingAnalyzers)
      muCandsMatchingAnalyzer->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, *itL1MuCand);
  }

  if(bestMuInBx1.isNonnull())
    etaGenPtGenBx1->Fill(bestMuInBx1->eta(), bestMuInBx1->pt() );
  if(bestMuInBx2.isNonnull())
    etaGenPtGenBx2->Fill(bestMuInBx2->eta(), bestMuInBx2->pt() );

  for(auto& rateAnalyser : rateAnalysers)
    rateAnalyser.fillHistos(event, MCTruthTTTrackHandle, bestMuInBx1, bestMuInBx2);

  int ttTracksFakesPerEventCnt = 0;
  int ttTracksPt10FakesPerEventCnt = 0;
  int ttTracksPt20FakesPerEventCnt = 0;
  for(unsigned int ttTrackIndx = 0; ttTrackIndx < TTTrackHandle->size(); ttTrackIndx++) {
    //the ttTracks that were not matched to any tracking particle, even loosely
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ttTrackPtr(TTTrackHandle, ttTrackIndx);
    //edm::Ptr< TrackingParticle > trackingParticle = MCTruthTTTrackHandle->findTrackingParticlePtr(ttTrackPtr);

    if ( (int)ttTrackPtr->getStubRefs().size() < L1Tk_minNStub) {
      continue;
    }

    float ttTrkPt = ttTrackPtr->getMomentum(L1Tk_nPar).perp();
    float ttTrkEta  = ttTrackPtr->getMomentum(L1Tk_nPar).eta();

    if(ttTrkPt > 20)
      LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<" "<<printTTTRack(ttTrackPtr, false, false);

    ttTrackEta_Pt->Fill(ttTrkEta, ttTrkPt);

    //do we really need the cut on the eta?
    if(abs(ttTrkEta) >= etaCutFrom && abs(ttTrkEta) <= etaCutTo) {
    }
    else
      continue;

    ttTracksPt->Fill(ttTrkPt);



    if(MCTruthTTTrackHandle->isLooselyGenuine(ttTrackPtr) == false) { //genuine is also loosely genuine,

      //float ttTrkPhi  = ttTrackPtr->getMomentum(L1Tk_nPar).phi();

      if(ttTrkPt > minMuPt )
        ttTracksFakesPerEventCnt++;

      if(ttTrkPt > 10) {
        ttTracksPt10FakesPerEventCnt++;
        ttTracksPt10FakesEta->Fill(ttTrkEta);
      }

      if(ttTrkPt > 20)
        ttTracksPt20FakesPerEventCnt++;

      ttTracksFakesPt->Fill(ttTrkPt);
    }

    l1t::BayesMuCorrelatorTrack dummy(ttTrackPtr); //very ugly, but should work - the idea is to reuse the MatchingAnalyzer
    ttTracksMatchingAnalyzer->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, dummy);
    ttTracksMatchingAnalyzerBarrel->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, dummy);
  }
  ttTracksFakesPerEvent->Fill(ttTracksFakesPerEventCnt);
  ttTracksPt10FakesPerEvent->Fill(ttTracksPt10FakesPerEventCnt);
  ttTracksPt20FakesPerEvent->Fill(ttTracksPt20FakesPerEventCnt);
}




DEFINE_FWK_MODULE(MuCorrelatorAnalyzer);

