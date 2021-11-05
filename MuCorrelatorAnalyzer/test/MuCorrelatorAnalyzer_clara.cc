//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
//#include "DataFormats/L1TMuon/interface/TkMuonBayesTrack.h"
//#include "L1Trigger/L1TkMuonBayes/plugins/L1TkMuonBayesTrackProducer.h"
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
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
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

//#include "L1Trigger/L1TkMuonBayes/interface/TkMuBayesProcConfig.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>
#include <stddef.h>

#define EDM_ML_LOGDEBUG
class TriggerAlgo {
public:
  TriggerAlgo(std::string name, double ptCut): name(name), ptCut(ptCut) {};
  virtual ~TriggerAlgo() {};
  
  virtual bool accept(const l1t::TrackerMuon& muCorrelatorTrack) = 0;

  std::string name;

  double ptCut = 0;
};
class SingleMuAlgo: public TriggerAlgo {
public:
  SingleMuAlgo(double ptCut): TriggerAlgo("SingleMuAlgo" + std::to_string((int)ptCut), ptCut) {};
  virtual ~SingleMuAlgo() {};

  virtual bool accept(const l1t::TrackerMuon& muCorrelatorTrack){
    return true;
  }
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

  //  virtual void takeCanidate(const edm::Ref<l1t::RegionalMuonCandBxCollection> itL1MuCand);
  virtual void takeCanidate(const l1t::TrackerMuon *itL1MuCand);

protected:
  std::shared_ptr<TriggerAlgo> triggerAlgo;

  l1t::TrackerMuon const* bestL1MuCand = nullptr;
  l1t::TrackerMuon const* notAcceptedL1MuCand = nullptr;

  double ptOfBestL1MuCand = -1;

  int numberOfAcceptedCandidates = 0;
};

void AnalyserBase::takeCanidate(const l1t::TrackerMuon *itL1MuCand) {
  if(triggerAlgo->accept(*itL1MuCand) ) {
    if(ptOfBestL1MuCand < itL1MuCand->pt() ) {
      if(ptOfBestL1MuCand > 0)
        edm::LogImportant("l1tMuBayesEventPrint") <<"L:"<<__LINE__<<" AnalyserBase::takeCanidate: another track already exists with pt "<<ptOfBestL1MuCand<<" new track has pt "<<itL1MuCand->pt() <<std::endl  ;
      numberOfAcceptedCandidates++;
      ptOfBestL1MuCand = itL1MuCand->pt();
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

    ptGenPtMuCandMuonsEv0HighPt = subDir.make<TH2I>("ptGenPtMuCandMuonsEv0HighPt", "ptGenPtMuCandMuonsEv0; pT gen [GeV]; pT ttTrack [GeV]; #", 100, 0, 1000, 100, 0, 1000);

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

  TH2I* ptGenPtMuCandMuonsEv0HighPt = nullptr;

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
      
      if(bestMatchedTTTrack.isNonnull() && bestMatchedTTTrack->momentum().perp() >= triggerAlgo->ptCut)
	ttMuonGenEtaMuons_withPtCuts->Fill(trackParticle->eta());
      
      if(ptOfBestL1MuCand > 0 && bestL1MuCand->pt() >= triggerAlgo->ptCut) {
	muCandGenEtaMuons_withPtCuts->Fill(trackParticle->eta());
	

      /*      auto& layerHitBits = bestL1MuCand->getFiredLayerBits(30);
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
      */
      }
    }
  
  int minMuPt = 3; //3 GeV
  if(ptOfBestL1MuCand > 0) {
    //    LogTrace("l1tMuBayesEventPrint")<<"\n"<<std::setw(32)<<triggerAlgo->name<<" best muCand track: "<<toString(*bestL1MuCand)<<endl<<endl;

    double muCandPt = bestL1MuCand->pt();
    if(trackParticle->pt() > minMuPt) {
      muCandGenEtaMuons->Fill(trackParticle->eta());
      muCandGenPhiMuons->Fill(trackParticle->phi());
    }

    muCandGenPtMuons->Fill(trackParticle->pt());

    if(trackParticle->eventId().event() == 0) {
      ptGenPtMuCandMuonsEv0->Fill(trackParticle->pt(), muCandPt);
      ptGenPtMuCandMuonsEv0HighPt->Fill(trackParticle->pt(), muCandPt);

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
      ptGenPtMuCandMuonsEv0HighPt->Fill(trackParticle->pt(), 0);

      if( abs(trackParticle->eta() ) < 0.82)
        ptGenPtMuCandMuonsEv0Barrel->Fill(trackParticle->pt(), 0);
      else if( abs(trackParticle->eta() ) < 1.24)
        ptGenPtMuCandMuonsEv0Overlap->Fill(trackParticle->pt(), 0);
      else if( abs(trackParticle->eta() ) < 2.4)
        ptGenPtMuCandMuonsEv0Endcap->Fill(trackParticle->pt(), 0);

      if(trackParticle->pt() > ptGenFrom) {
        LogTrace("l1tMuBayesEventPrint")<<"\n"<<std::setw(32)<<triggerAlgo->name<<" no muCand track !!!!"<<endl;
      }
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

  std::vector<EfficiencyAnalyser> efficiencyAnalysers;
 
  //std::vector<EfficiencyAnalyser> efficiencyAnalysersCorrelatorWithTrackPart; //used when the correlator works with the Tracking particles and not the ttTracks

  edm::ParameterSet parameterSet;
  unsigned int eventCount;
  //TFile* myRootFile;
  double etaCutFrom = 0;
  double etaCutTo = 2.4;

  int minMuPt = 3; //3 GeV



  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)
  bool LooseMatch;      // use loose MC-matching instead
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt
  double TP_maxEta;     // save TPs with |eta| < maxEta
  double TP_maxZ0;      // save TPs with |z0| < maxZ0

  double TP_maxRho = 30; //[cm] maximum muon vertex rho - to not include in the efficiency analysis the muons from pions, which are not well matched to ttTracks. increase for not pointin muon analysis

  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)

  int muCandQualityCut = 12;

  string analysisType = "efficiency";

  //bool applyTriggerRules = false; //if true, for the efficiency plots the candidates that has getCanceledByTriggerRules() = true are not counted

  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag L1TkMuon;
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
  edm::EDGetTokenT< std::vector< l1t::TrackerMuon > > MuTrackToken_;
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

  //TH1I* simTracksBetaPt20 = nullptr;
  //TH1I* simTracksPt = nullptr;

  TH1I* tpBetaPt20 = nullptr;

  TH1I* tpBetaHscpPt20 = nullptr;

  TH2I* etaGenPtGenBx1 = nullptr;

  TH2I* etaGenPtGenBx2 = nullptr;

  TH2I* ttTracksPerMuonTPvsPtGen = nullptr;

  edm::EDGetTokenT<l1t::TrackerMuonCollection> inputMuCorr;
  //edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> vertexSim;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;
};

MuCorrelatorAnalyzer::MuCorrelatorAnalyzer(const edm::ParameterSet& conf)
: parameterSet(conf), eventCount(0)
{
  inputMuCorr = consumes<l1t::TrackerMuonCollection>(edm::InputTag("L1TkMuon")); //


  //simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //TODO which is correct?
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
  TP_minNStub      = parameterSet.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = parameterSet.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = parameterSet.getParameter< double >("TP_minPt");
  TP_maxEta        = parameterSet.getParameter< double >("TP_maxEta");
  TP_maxZ0         = parameterSet.getParameter< double >("TP_maxZ0");
  TP_maxRho        = parameterSet.getParameter< double >("TP_maxRho");
  L1TrackInputTag      = parameterSet.getParameter<edm::InputTag>("L1TrackInputTag");
  //L1TkMuon      = parameterSet.getParameter<edm::InputTag>("L1TkMuon");
  MCTruthTrackInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub    = parameterSet.getParameter< int >("L1Tk_minNStub");

  muCandQualityCut = parameterSet.getParameter< int >("muCandQualityCut");
  analysisType = parameterSet.getParameter< string >("analysisType");

  //applyTriggerRules = parameterSet.getParameter< bool >("applyTriggerRules");

  L1StubInputTag      = parameterSet.getParameter<edm::InputTag>("L1StubInputTag");
  //MCTruthClusterInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  //MCTruthStubInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthStubInputTag");
  //TrackingParticleInputTag = parameterSet.getParameter<edm::InputTag>("TrackingParticleInputTag");
  //TrackingVertexInputTag = parameterSet.getParameter<edm::InputTag>("TrackingVertexInputTag");

  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);

  ttStubToken_ = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  //  MuTrackToken_ = consumes< std::vector< l1t::TrackerMuon > >(L1TkMuon);
  TrackingVertexToken_ = consumes< std::vector< TrackingVertex > >(TrackingVertexInputTag);


  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
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
  /*  std::shared_ptr<TriggerAlgo> singleMuAlgoSoftCuts = std::make_shared<SingleMuAlgoSoftCuts>(20);

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
  */

  if(analysisType == "rate") {
    /*    rateAnalysers.emplace_back(singleMuAlgo, fs);

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
    */
  }
  else if(analysisType == "efficiency") {
    efficiencyAnalysers.emplace_back(singleMuAlgo, 25, 10000, fs);

    /*    efficiencyAnalysers.emplace_back(singleMuAlgoSoftCuts, 25, 10000, fs);
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

    //efficiencyAnalysers.emplace_back(hscpAlgo20, 25, 10000, fs);

    //efficiencyAnalysers.emplace_back(hscpAlgo30, fs);

    //efficiencyAnalysers.emplace_back(hscpAlgoHardCuts20, fs);
    //efficiencyAnalysers.emplace_back(hscpAlgoSoftCuts20, 25, 10000, fs);
    //efficiencyAnalysers.emplace_back(hscpAlgoPdfSumCuts20, 25, 10000, fs);
    */
  }
  
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoPtCut10, fs));
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoBarrelPtCut10, fs));
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgo, fs));
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoSoftCuts, fs));
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoSoftCuts1, fs));
  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoPdfSumSoftCuts, fs));

  //muCandsMatchingAnalyzers.emplace_back(std::make_unique<MuCandsMatchingAnalyzer>(singleMuAlgoHardCuts, fs));
  //muCandsMatchingAnalyzers.emplace_back(hscpAlgo20, fs);


  
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

  
  //simTracksBetaPt20 = fs->make<TH1I>("simTracksBetaPt20", "simTracks beta, pt > 20 GeV; beta; #events",  20, 0., 1.);
  //simTracksPt = fs->make<TH1I>("simTracksPt", "simTracksPt; ptGen [GeV]; #events", 100, 0, 1000);

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

/*void MuCorrelatorAnalyzer::hscpAnalysis(edm::Handle< std::vector< TrackingParticle > >& trackingParticleHandle, edm::Handle<l1t::TrackerMuonBxCollection>& muCorrTracksHandle) {
  for (unsigned int iTrackPart = 0; iTrackPart != trackingParticleHandle->size(); iTrackPart++ ) {
    edm::Ptr< TrackingParticle > trackPartPtr(trackingParticleHandle, iTrackPart);
    hscpAnalysis(trackPartPtr, muCorrTracksHandle);
  }
}*/

void MuCorrelatorAnalyzer::analyze(
    const edm::Event& event, const edm::EventSetup& es)
{
  //const unsigned int omtflayersCnt = muCorrConfig.nLayers();
  LogTrace("l1tMuBayesEventPrint") << "\n\nMuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" new event "<<event.id()<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< endl;


  //edm::Handle<edm::SimTrackContainer> simTraksHandle;
  //event.getByToken(simTrackToken, simTraksHandle);

  //const std::vector<SimTrack>& simTraks = *(simTraksHandle.product());

  //edm::Handle<edm::SimVertexContainer> simVx;
  //event.getByToken(vertexSim, simVx);
  //const std::vector<SimVertex>& simVertexes = *(simVx.product());

  edm::Handle<reco::GenParticleCollection> genPartHandle;
  event.getByToken(genParticleToken, genPartHandle);

  //LogTrace("l1tMuBayesEventPrint") << "\nMuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" vertexSim "<<simVx->size()<< endl;

  //analyzeSimTracks(simTraksHandle, genPartHandle); //muCorrTracksHandle


  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  event.getByToken(ttTrackToken_, TTTrackHandle);

  // MC truth association maps

  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  event.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > trackingParticleHandle;
  //  edm::Handle< std::vector< l1t::TrackerMuon > > tkmuoneHandle;
  //edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  event.getByToken(TrackingParticleToken_, trackingParticleHandle);
  //  event.getByToken(MuTrackToken_, tkmuoneHandle);
  //event.getByToken(TrackingVertexToken_, TrackingVertexHandle);

  //std::cout <<" L1 MUONS: "<<std::endl;
  //edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  //event.getByToken(inputOMTF, l1Omtf);
  //cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" omtfCands->size(bxNumber) "<<omtfCands->size()<< endl;

  edm::Handle<l1t::TrackerMuonCollection> muCorrTracksHandle;
  event.getByToken(inputMuCorr, muCorrTracksHandle);
  auto muCorrTracks = muCorrTracksHandle.product();

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
    ostr << "matched ttTrack       pt: " <<setw(7)<< ttTrackPtr->momentum().perp()
                     <<" charge: "<<setw(2)<<(ttTrackPtr->rInv() > 0 ? 1 : -1)
                     << " eta: " <<setw(7)<< ttTrackPtr->momentum().eta()
                     << " phi: " <<setw(7)<< ttTrackPtr->momentum().phi()
                     << " chi2: " <<setw(7)<< ttTrackPtr->chi2()
                     << " consistency: " <<setw(7)<< ttTrackPtr->stubPtConsistency()
                     << " z0: " <<setw(7)<< ttTrackPtr->POCA().z()
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
        //tpBetaHscpPt20->Fill(tpPtr->p4().Beta());
      }
    }

    if(analysisType == "withTrackPart")
      //hscpAnalysis(event, tpPtr, muCorrTracksHandle);

    if (abs(tpPtr->pdgId()) == 1000015) {
      //hscpGenEta->Fill(tpPtr->momentum().eta());
      //hscpGenPt->Fill(tpPtr->pt());
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

        if (bestMatchedTTTrack.isNull() || matchedTTTrack->chi2Red() < chi2dofOfBestMatchedTrack) {
          bestMatchedTTTrack = matchedTTTrack;
          chi2dofOfBestMatchedTrack = matchedTTTrack->chi2Red();
        }

        //finding omtf candidate corresponding to the matchedTTTrack
        //LogTrace("l1tMuBayesEventPrint")<<" omtfCands->size "<<omtfCands->size(bxNumber)<<endl;

        //we should find one candidate with max pt and use it for filling the histograms
        for(auto itL1MuCand = muCorrTracks->begin(); itL1MuCand != muCorrTracks->end(); ++itL1MuCand) {
          auto& muCandTtTrackPtr = itL1MuCand->trkPtr();
          //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > omtfTtTrackPtr(TTTrackHandle, omtfTTTrackIndex);

          if(muCandTtTrackPtr.isNonnull() && muCandTtTrackPtr.get() == matchedTTTrack.get() ) { //it should be only one omtfTrack for given ttTRack
	    //            LogTrace("l1tMuBayesEventPrint")<<"L:"<<__LINE__<<toString(*itL1MuCand)<<endl;

            for(auto& effAnalys : efficiencyAnalysers)
              effAnalys.takeCanidate(&(*itL1MuCand));

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
      float matchTTTrackPt   = bestMatchedTTTrack->momentum().perp();
      float matchTTTrackEta  = bestMatchedTTTrack->momentum().eta();
      float matchTTTrackPhi  = bestMatchedTTTrack->momentum().phi();

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
          if( abs(matchedTTTrack->momentum().phi() - tpPtr->phi()) < 0.01 &&
              abs(matchedTTTrack->momentum().eta() - tpPtr->eta()) < 0.01 &&
              abs( (matchedTTTrack->momentum().perp() - tpPtr->pt() ) / tpPtr->pt()) < 0.1)
          {
            if(MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTTTrack).isNull() ) { //we require that this ttTrack is not well matched to other particle
              ttMuonVeryLoosePt->Fill(tpPtr->pt());
              if(tpPtr->pt() > 20 && matchedTTTrack->momentum().perp() >= 18 ) {
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

          float ttTrkPt   = iterL1Track->getMomentum().perp();
          float ttTrkEta  = iterL1Track->getMomentum().eta();
          float ttTrkPhi  = iterL1Track->getMomentum().phi();
          float tmp_trk_z0   = iterL1Track->getPOCA().z(); //cm

          float tmp_trk_chi2 = iterL1Track->getChi2();
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

  //  for(auto& rateAnalyser : rateAnalysers)
  //  rateAnalyser.reset();

  for(auto itL1MuCand = muCorrTracks->begin();
        itL1MuCand != muCorrTracks->end(); ++itL1MuCand) {

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

    if(abs(itL1MuCand->eta()) >= etaCutFrom && abs(itL1MuCand->eta() ) <= etaCutTo) {
    }
    else
      continue;

    /*    for(auto& rateAnalyser : rateAnalysers)
      rateAnalyser.takeCanidate(itL1MuCand);
    */
    //    for(auto& muCandsMatchingAnalyzer : muCandsMatchingAnalyzers)
    //      muCandsMatchingAnalyzer->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, *itL1MuCand);
  }

  if(bestMuInBx1.isNonnull())
    etaGenPtGenBx1->Fill(bestMuInBx1->eta(), bestMuInBx1->pt() );
  if(bestMuInBx2.isNonnull())
    etaGenPtGenBx2->Fill(bestMuInBx2->eta(), bestMuInBx2->pt() );

  //  for(auto& rateAnalyser : rateAnalysers)
  // rateAnalyser.fillHistos(event, MCTruthTTTrackHandle, bestMuInBx1, bestMuInBx2);

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

    float ttTrkPt = ttTrackPtr->momentum().perp();
    float ttTrkEta  = ttTrackPtr->momentum().eta();

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

      //float ttTrkPhi  = ttTrackPtr->getMomentum().phi();

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

    //l1t::TrackerMuon dummy&(*ttTrackPtr); //very ugly, but should work - the idea is to reuse the MatchingAnalyzer
    //ttTracksMatchingAnalyzer->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, dummy);
    //ttTracksMatchingAnalyzerBarrel->fillHistos(event, MCTruthTTTrackHandle, muonTrackingParticles, dummy);
  }
  ttTracksFakesPerEvent->Fill(ttTracksFakesPerEventCnt);
  ttTracksPt10FakesPerEvent->Fill(ttTracksPt10FakesPerEventCnt);
  ttTracksPt20FakesPerEvent->Fill(ttTracksPt20FakesPerEventCnt);
}




DEFINE_FWK_MODULE(MuCorrelatorAnalyzer);
