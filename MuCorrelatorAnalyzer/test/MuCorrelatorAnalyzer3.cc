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

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>
#include <stddef.h>

#define EDM_ML_LOGDEBUG

using namespace std;

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

  edm::ParameterSet parameterSet;
  unsigned int eventCount;
  //TFile* myRootFile;
  double etaCutFrom = 0;
  double etaCutTo = 2.35;

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
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)

  int muCandQualityCut = 12;

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
  TH1D* gpPerEvent = nullptr;
  TH1D* ttTrackPerEvent = nullptr;

  TH1D* gpMuonPerEvent = nullptr;
  TH1D* ttTrackMuonPerEvent = nullptr;

  TH1D* ttTrackPerEventEta2_8 = nullptr;
  TH1D* ttTrackMuonPerEventEta2_8 = nullptr;

  TH1D* ttTrackInOmtfRegPerEvent = nullptr;
  TH1D* ttTrackMuonInOmtfRegPerEvent = nullptr;

//denominators
  TH1D* gpMuonPt = nullptr;
  TH1D* gpMuonEta = nullptr;
  TH1D* gpMuonEta_ptGen20GeV = nullptr;
  TH1D* gpMuonPhi = nullptr;

  TH1D* ttMuonPt = nullptr;
  TH1D* ttMuonEta = nullptr;
  TH1D* ttMuonEta_ptGen20GeV_ptTT18Gev = nullptr;
  TH1D* ttMuonPhi = nullptr;

  TH1D* omtfTTTrackPt = nullptr;
  TH1D* omtfTTTrackEta = nullptr;
  TH1D* omtfTTTrackPhi = nullptr;

  TH1D* omtfTTTrackMuonPt = nullptr;
  TH1D* omtfTTTrackMuonEta = nullptr;
  TH1D* omtfTTTrackMuonPhi = nullptr;

 //nominator
  TH1D* omtfAndTtMuonPt = nullptr;
  TH1D* omtfAndTtMuonEta = nullptr;
  TH1D* omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev = nullptr;
  TH1D* omtfAndTtMuonPhi = nullptr;


  TH1D* fakeTrackOmtfPt = nullptr;
  //TH1D* fakeTrackOmtfEta = nullptr;
  TH1D* fakeTrackOmtfEta_ptTT18Gev = nullptr;
  TH1D* fakeTrackOmtfPhi_ptTT18Gev = nullptr;

  TH1D* wrongTagOmtfPt = nullptr;
  //TH1D* wrongTagOmtfEta = nullptr;
  TH1D* wrongTagOmtfEta_ptTT18Gev = nullptr;
  TH1D* wrongTagOmtfPhi_ptTT18Gev = nullptr;

  TH1D* lostTtMuonPt = nullptr;
  TH1D* lostTtMuonEta_ptGen20GeV = nullptr;
  TH1D* lostTtMuonPhi_ptGen20GeV = nullptr;

  TH2I* ptGenPtTTMuon = nullptr;
  TH2I* ptGenPtOMtfMuon = nullptr;

  TH2I* ptGenPtTTMuonEv0 = nullptr;
  TH2I* ptGenPtOMtfMuonEv0 = nullptr;

  TH2I* ptGenPtTTMuonEvPu = nullptr;
  TH2I* ptGenPtOMtfMuonEvPu = nullptr;

  TH2I* ptGenDeltaPtTTMuon = nullptr;
  TH2I* ptGenDeltaPhiTTMuon = nullptr;
  TH2I* ptGenDeltaEtaTTMuon = nullptr;

  TH1D* hRefLayer = nullptr;

  TH2I* muonsPdfSumFiredPlanes = nullptr;
  TH2I* notMuonsPdfSumFiredPlanes = nullptr;

  TH1D* chi2GenuineTTTrackOmtf = nullptr;
  TH1D* chi2FakeTTTrackOmtf = nullptr;

  TH1D* chi2DofGenuineTTTrackOmtf = nullptr;
  TH1D* chi2DofFakeTTTrackOmtf = nullptr;


  /*  TH1D* hRefLayerCut;
  TH1D* hRefLayerCut2;
  TH1D* hLayerHits;
  TH1D* hLayerHitsCut;
  TH1D* hLayerHitsCut2;
  TH1D* hQuality;
  TH1D* hNLayers;
  TH1D* hRefL1Layers;
  TH1D* hRefL3Layers;

  TH1D* hEtaSimCut;
  TH1D* hEffVsEtaSim;
  TH1D* hEffVsEtaSimOmtfCut;

  vector<TH2D*> hPtSimPtVsRefLay;*/

  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> inputOMTF;
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
};


MuCorrelatorAnalyzer::MuCorrelatorAnalyzer(const edm::ParameterSet& conf)
: parameterSet(conf), eventCount(0)
{
  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simBayesMuCorrelatorTrackProducer", "MuCorr")); //simOmtfDigis
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));

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
  L1TrackInputTag      = parameterSet.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = parameterSet.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub    = parameterSet.getParameter< int >("L1Tk_minNStub");

  muCandQualityCut = parameterSet.getParameter< int >("muCandQualityCut");

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

  //make a new Root file
/*  string  outRootFile = parameterSet.getParameter<std::string>("outRootFile");
  myRootFile = new TFile(outRootFile.c_str(),"RECREATE");
  if(myRootFile == 0 || myRootFile-> IsZombie() ) {
    cout<<"file "<<outRootFile<<" not opened"<<endl;
    throw cms::Exception("OMTFConfiguration::getPatternPtRange: patternPts vector not initialized");

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

  gpPerEvent = fs->make<TH1D>("gpPerEvent", "generated particles per event, |#eta| < 2.4; N_tracks; #events", 500, 0., 2000.);
  ttTrackPerEvent = fs->make<TH1D>("ttTrackPerEvent", "ttTracks per event, |#eta| < 2.4; N_tracks; #events", 500, 0., 1000.);

  gpMuonPerEvent = fs->make<TH1D>("gpMuonPerEvent", "generated muons per event, |#eta| < 2.4; N_tracks; #events", 50, 0., 50.);
  ttTrackMuonPerEvent = fs->make<TH1D>("ttTrackMuonPerEvent", "ttTracks matched to muons per event; |#eta| < 2.4, N_tracks; #events", 50, 0., 50.);

  ttTrackPerEventEta2_8 = fs->make<TH1D>("ttTrackPerEventEta2_8", "ttTracks per event, |#eta| < 2.8; N_tracks; #events", 500, 0., 500.);
  ttTrackMuonPerEventEta2_8 = fs->make<TH1D>("ttTrackMuonPerEventEta2_8", "ttTracks per event, |#eta| < 2.8; N_tracks; #events", 50, 0., 50.);

  ttTrackInOmtfRegPerEvent = fs->make<TH1D>("ttTrackInOmtfRegPerEvent", "ttTracks in OMTF region per event; N_tracks; #events", 500, 0., 500.);
  ttTrackMuonInOmtfRegPerEvent =fs->make<TH1D>("ttTrackMuonInOmtfRegPerEvent", "ttTracks matched to muons in OMTF region per event; N_tracks; #events", 50, 0., 50.);


  //denominators
  gpMuonPt = fs->make<TH1D>("gpMuonPt", "gpMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);
  gpMuonEta = fs->make<TH1D>("gpMuonEta", "gpMuonEta; eta; #events", etaBins, -2.4, 2.4);
  gpMuonEta_ptGen20GeV = fs->make<TH1D>("gpMuonEta_ptGen20GeV", "gpMuonEta_ptGen20GeV; eta; #events", etaBins, -2.4, 2.4);
  gpMuonPhi = fs->make<TH1D>("gpMuonPhi", "gpMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  ttMuonPt = fs->make<TH1D>("ttMuonPt", "ttMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);
  ttMuonEta = fs->make<TH1D>("ttMuonEta", "ttMuonEta; eta; #events", etaBins, -2.4, 2.4);
  ttMuonEta_ptGen20GeV_ptTT18Gev = fs->make<TH1D>("ttMuonEta_ptGen20GeV_ptTT18Gev", "ttMuonEta_ptGen20GeV_ptTT18Gev; eta; #events", etaBins, -2.4, 2.4);
  ttMuonPhi = fs->make<TH1D>("ttMuonPhi", "ttMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  omtfTTTrackPt = fs->make<TH1D>("omtfTTTrackPt", "omtfTTTrackPt; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
  omtfTTTrackEta = fs->make<TH1D>("omtfTTTrackEta", "omtfTTTrackEta; eta; #events", etaBins, -2.4, 2.4);
  omtfTTTrackPhi = fs->make<TH1D>("omtfTTTrackPhi", "omtfTTTrackPhi; phi; #events", phiBins, -M_PI, M_PI);

  omtfTTTrackMuonPt = fs->make<TH1D>("omtfTTTrackMuonPt", "omtfTTTrackMuonPt; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
  omtfTTTrackMuonEta = fs->make<TH1D>("omtfTTTrackMuonEta", "omtfTTTrackMuonEta; eta; #events", etaBins, -2.4, 2.4);
  omtfTTTrackMuonPhi = fs->make<TH1D>("omtfTTTrackMuonPhi", "omtfTTTrackMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  fakeTrackOmtfPt = fs->make<TH1D>("fakeTrackOmtfPt", "fakeTrackOmtfPt; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
  fakeTrackOmtfEta_ptTT18Gev = fs->make<TH1D>("fakeTrackOmtfEta_ptTT18Gev", "fakeTrackOmtfEta_ptTT18Gev; eta; #events", etaBins, -2.4, 2.4);
  fakeTrackOmtfPhi_ptTT18Gev = fs->make<TH1D>("fakeTrackOmtfPhi_ptTT18Gev", "fakeTrackOmtfPhi_ptTT18Gev; phi; #events", phiBins, -M_PI, M_PI);

  wrongTagOmtfPt = fs->make<TH1D>("wrongTagOmtfPt", "wrongTagOmtfPt; ttTrack pt [GeV]; #events", ptBins, 0., 500.);;
  wrongTagOmtfEta_ptTT18Gev = fs->make<TH1D>("wrongTagOmtfEta_ptTT18Gev", "wrongTagOmtfEta_ptTT18Gev; eta; #events", etaBins, -2.4, 2.4);
  wrongTagOmtfPhi_ptTT18Gev = fs->make<TH1D>("wrongTagOmtfPhi_ptTT18Gev", "wrongTagOmtfPhi_ptTT18Gev; phi; #events", phiBins, -M_PI, M_PI);

  //nominators
  //versus gen (tracking particle) pt, eta, phi
  omtfAndTtMuonPt = fs->make<TH1D>("omtfAndTtMuonPt", "omtfAndTtMuonPt; gen pT [GeV]; #events", ptBins, 0., 500.);;
  omtfAndTtMuonEta = fs->make<TH1D>("omtfAndTtMuonEta", "omtfAndTtMuonEta; eta; #events", etaBins, -2.4, 2.4);
  omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev = fs->make<TH1D>("omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev", "omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev; eta; #events", etaBins, -2.4, 2.4);
  omtfAndTtMuonPhi = fs->make<TH1D>("omtfAndTtMuonPhi", "omtfAndTtMuonPhi; phi; #events", phiBins, -M_PI, M_PI);

  lostTtMuonPt = fs->make<TH1D>("lostTtMuonPt", "lostTtMuonPt; ttMuonPt [GeV]; #events", ptBins, 0., 500.);;
  lostTtMuonEta_ptGen20GeV = fs->make<TH1D>("lostTtMuonEta_ptGen20GeV", "lostTtMuonEta_ptGen20GeV; eta; #events", etaBins, -2.4, 2.4);
  lostTtMuonPhi_ptGen20GeV = fs->make<TH1D>("lostTtMuonPhi_ptGen20GeV", "lostTtMuonPhi_ptGen20GeV; phi; #events", phiBins, -M_PI, M_PI);


  ptGenPtTTMuon = fs->make<TH2I>("ptGenPtTTMuon", "ptGenPtTTMuon; gen pT [GeV]; ttTrack pT [GeV]; #", 100, 0, 100, 100, 0, 100);
  ptGenPtOMtfMuon = fs->make<TH2I>("ptGenPtOMtfMuon", "ptGenPtOMtfMuon; gen pT [GeV]; OMTF Track pT [GeV]; #", 100, 0, 100, 100, 0, 100);

  ptGenPtTTMuonEv0 = fs->make<TH2I>("ptGenPtTTMuonEv0", "ptGenPtTTMuonEv0; gen pT [GeV]; ttTrack pT [GeV]; #", 100, 0, 100, 100, 0, 100);
  ptGenPtOMtfMuonEv0 = fs->make<TH2I>("ptGenPtOMtfMuonEv0", "ptGenPtOMtfMuonEv0; gen pT [GeV]; OMTF Track pT [GeV]; #", 100, 0, 100, 100, 0, 100);

  ptGenPtTTMuonEvPu = fs->make<TH2I>("ptGenPtTTMuonEvPu", "ptGenPtTTMuonEvPu; gen pT [GeV]; ttTrack pT [GeV]; #", 100, 0, 100, 100, 0, 100);
  ptGenPtOMtfMuonEvPu = fs->make<TH2I>("ptGenPtOMtfMuonEvPu", "ptGenPtOMtfMuonEvPu; gen pT [GeV]; OMTF Track pT [GeV]; #", 100, 0, 100, 100, 0, 100);


  ptGenDeltaPtTTMuon = fs->make<TH2I>("ptGenDeltaPtTTMuon", "ptGenDeltaPtTTMuon; gen pT [GeV]; (gen pT - ttTrack pT) [GeV]; #", 100, 0, 100, 100, -5, 5);
  ptGenDeltaPhiTTMuon = fs->make<TH2I>("ptGenDeltaPhiTTMuon", "ptGenDeltaPhiTTMuon; gen pT [GeV]; (gen phi - ttTrack phi); #", 100, 0, 100,  100, -0.5, 0.5);
  ptGenDeltaEtaTTMuon = fs->make<TH2I>("ptGenDeltaEtaTTMuon", "ptGenDeltaEtaTTMuon; gen pT [GeV]; (gen eta - ttTrack eta); #", 100, 0, 100,  100, -0.1, 0.1);

  muonsPdfSumFiredPlanes = fs->make<TH2I>("muonsPdfSumFiredPlanes", "muonsPdfSumFiredPlanes; pdfSum; FiredPlanes; #", 100, 0, 10000, 19, -0.5, 18.5);
  notMuonsPdfSumFiredPlanes = fs->make<TH2I>("notMuonsPdfSumFiredPlanes", "notMuonsPdfSumFiredPlanes; pdfSum; FiredPlanes; #", 100, 0, 10000, 19, -0.5, 18.5);

  chi2GenuineTTTrackOmtf = fs->make<TH1D>("chi2GenuineTTTrackOmtf", "chi2GenuineTTTrackOmtf; chi2; #events", 40, 0., 200.);
  chi2FakeTTTrackOmtf = fs->make<TH1D>("chi2FakeTTTrackOmtf", "chi2FakeTTTrackOmtf; chi2; #events", 40, 0., 200.);

  chi2DofGenuineTTTrackOmtf = fs->make<TH1D>("chi2DofGenuineTTTrackOmtf", "chi2DofGenuineTTTrackOmtf; chi2; #events", 40, 0., 200.);
  chi2DofFakeTTTrackOmtf = fs->make<TH1D>("chi2DofFakeTTTrackOmtf", "chi2DofFakeTTTrackOmtf; chi2; #events", 40, 0., 200.);

  hRefLayer = fs->make<TH1D>("refLayer", "refLayer; refLayer; #events", 10, -0.5, 10-0.5);

  edm::LogImportant("MuCorrelatorAnalyzer")<< "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;
}

void MuCorrelatorAnalyzer::endJob()
{
  cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;
}

void MuCorrelatorAnalyzer::analyze(
    const edm::Event& event, const edm::EventSetup& es)
{
  LogTrace("l1tMuBayesEventPrint") << "\nOmtfTTAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" new event "<<event.id()<< endl;

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
  edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  event.getByToken(inputOMTF, l1Omtf);
  auto omtfCands = l1Omtf.product();
  //cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ <<" omtfCands->size(bxNumber) "<<omtfCands->size()<< endl;

  int bxNumber = 0;



  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------

  //LogTrace("l1tMuBayesEventPrint") << endl << "Loop over tracking particles!" << endl;

  gpPerEvent->Fill(trackingParticleHandle->size() );
  ttTrackPerEvent->Fill(TTTrackHandle->size());

  int minMuPt = 3; //3 GeV
  int allGpMuCnt = 0;
  int ttTrackMuCnt = 0;
  int ttTrackInOmtfCnt = 0;
  int ttTrackMuInOmtfCnt = 0;

  //boost::dynamic_bitset omtfCandsIsMatchedToMuonTrakPart(omtfCands->size(bxNumber));
  LogTrace("l1tMuBayesEventPrint")<<"trackingParticleHandle->size() "<<trackingParticleHandle->size()<<" isValid "<<trackingParticleHandle.isValid();

  auto printTrackigParticleShort = [&](edm::Ptr< TrackingParticle >& tpPtr) {
    std::ostringstream ostr;
    ostr<< "Tracking particle,    pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
        << " eta: " <<setw(7)<< tpPtr->eta()
        << " phi: " <<setw(7)<< tpPtr->phi()
        << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
        << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr).size()
        <<" bx: "<<tpPtr->eventId().bunchCrossing()
        <<dec<<" key "<<tpPtr.key()<< endl;
    return ostr.str();
  };

  auto printTTTRack = [&](edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > >& ttTrackPtr, bool genuine, bool loosegenuine) {
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

  for (unsigned int iTP = 0; iTP < trackingParticleHandle->size(); ++iTP) {
    edm::Ptr< TrackingParticle > tpPtr(trackingParticleHandle, iTP);

    if(tpPtr->eventId().bunchCrossing() != 0)
      continue;

    //int tmp_eventid = tp_ptr->eventId().event();
    //if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tp_d0_prod = -tpPtr->vx()*sin(tpPtr->phi()) + tpPtr->vy()*cos(tpPtr->phi());

    if (abs(tpPtr->pdgId()) == 13 ) {  //|| tpPtr->pt() > 20
      //only muons
    }
    else
      continue;

    if (tpPtr->pt() < TP_minPt)
      continue;
    if (fabs(tpPtr->eta()) > TP_maxEta)
      continue;

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

    if (fabs(tmp_tp_z0) > TP_maxZ0)
      continue;

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


    auto printTrackigParticle = [&]() {
      std::ostringstream ostr;
      ostr<< "Tracking particle,    pt: " <<setw(7)<< tpPtr->pt() <<" charge: "<<setw(2)<<tpPtr->charge()
          << " eta: " <<setw(7)<< tpPtr->eta()
          << " phi: " <<setw(7)<< tpPtr->phi()
          << " pdgid: " <<setw(7)<< tpPtr->pdgId() << " eventID: " <<setw(7)<< tpPtr->eventId().event()
          << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(tpPtr).size()
          << " bx: "<<tpPtr->eventId().bunchCrossing()
          << dec<<" key "<<tpPtr.key()
          << " z0: " <<setw(7)<< tmp_tp_z0 << " d0: " <<setw(7)<< tmp_tp_d0
          << " z_prod: " <<setw(7)<< tpPtr->vz() << " d_prod: " <<setw(7)<< tp_d0_prod
          << endl;
      return ostr.str();
    };

    if(abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
      gpMuonPt->Fill(tpPtr->pt());
      gpMuonPhi->Fill(tpPtr->phi());

      LogTrace("l1tMuBayesEventPrint")<<"\n\n"<<printTrackigParticle();
    }

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
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > bestMatchedTTTrack; //null at the begining
    float chi2dofOfBestMatchedTrack = 99999;

    bool isOmtf = false;
    // ----------------------------------------------------------------------------------------------
    // loop over matched L1 tracks
    // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters
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
      if (LooseMatch && !tmp_trk_loosegenuine)
        continue;
      if (!LooseMatch && !tmp_trk_genuine)
        continue;

      edm::Ptr< TrackingParticle > tpOfMatchedTTTrack = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTTTrack);
      if (abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
        LogTrace("l1tMuBayesEventPrint")<<printTTTRack(matchedTTTrack, tmp_trk_genuine, tmp_trk_loosegenuine);
        if (tpOfMatchedTTTrack.isNull()) { //if the ttTrack is at least loosegenuine this cannot happened
          LogTrace("l1tMuBayesEventPrint") << "track matched to TP is NOT uniquely matched to a TP" << endl;
        }
        else {
          /*LogTrace("l1tMuBayesEventPrint") << "TP matched to ttTrack pt: " <<setw(7)<< tpOfMatchedTTTrack->p4().pt() << " eta: " <<setw(7)<< tpOfMatchedTTTrack->momentum().eta()
                     << " phi: " <<setw(7)<< tpOfMatchedTTTrack->momentum().phi() << " z0: " <<setw(7)<< tpOfMatchedTTTrack->vertex().z() << endl;*/
          LogTrace("l1tMuBayesEventPrint")<<printTrackigParticleShort(tpOfMatchedTTTrack);
        }
      }
      // ----------------------------------------------------------------------------------------------
      // number of stubs in this matchedTTTrack
      int ttTrackNstub = matchedTTTrack->getStubRefs().size();

      if (ttTrackNstub < L1Tk_minNStub) {
        LogTrace("l1tMuBayesEventPrint") << "ttTrackNstub < "<<L1Tk_minNStub<<" - skipping this ttTrack" << endl;
        continue;
      }

      //dmatch_pt  = fabs(my_tp->p4().pt() - tpPtr->pt());
      /*float dmatch_eta = fabs(tpOfMatchedTTTrack->p4().eta() - tpPtr->eta());
      float dmatch_phi = fabs(tpOfMatchedTTTrack->p4().phi() - tpPtr->phi());
      float dmatch_id = tpOfMatchedTTTrack->pdgId();*/

      float tmp_trk_chi2dof = (matchedTTTrack->getChi2(L1Tk_nPar)) / (2*ttTrackNstub - L1Tk_nPar);

      // ensure that track is uniquely matched to the TP we are looking at!
      if(tpPtr != tpOfMatchedTTTrack) {
        LogTrace("l1tMuBayesEventPrint")<<"\ntpPtr != tpOfMatchedTTTrack !!!!!!!!!!!!";
        LogTrace("l1tMuBayesEventPrint")<<"current tracking particle";
        LogTrace("l1tMuBayesEventPrint")<<printTrackigParticleShort(tpPtr);

        LogTrace("l1tMuBayesEventPrint")<<printTTTRack(matchedTTTrack, tmp_trk_genuine, tmp_trk_loosegenuine);
        LogTrace("l1tMuBayesEventPrint")<<"tpOfMatchedTTTrack";
        LogTrace("l1tMuBayesEventPrint")<<printTrackigParticleShort(tpOfMatchedTTTrack);
      }
      else {
        //so we assume that the matched ttTrack (and thus omtf track) is at least loose genuine
        nMatch++;
        if (bestMatchedTTTrack.isNull() || tmp_trk_chi2dof < chi2dofOfBestMatchedTrack) {
          bestMatchedTTTrack = matchedTTTrack;
          chi2dofOfBestMatchedTrack = tmp_trk_chi2dof;
        }

        //finding omtf candidate corresponding to the matchedTTTrack
        //LogTrace("l1tMuBayesEventPrint")<<" omtfCands->size "<<omtfCands->size(bxNumber)<<endl;
        for(l1t::RegionalMuonCandBxCollection::const_iterator itOmtfCand = omtfCands->begin(bxNumber);
            itOmtfCand != l1Omtf.product()->end(bxNumber); ++itOmtfCand) {

          /* LogTrace("l1tMuBayesEventPrint")<<" itOmtfCand->trackAddress().size "<<itOmtfCand->trackAddress().size()<<endl;
            LogTrace("l1tMuBayesEventPrint")<<" itOmtfCand->trackAddress().at(0) "<<itOmtfCand->trackAddress().at(0)<<endl;
            LogTrace("l1tMuBayesEventPrint")<<" itOmtfCand->trackAddress().at(1) "<<itOmtfCand->trackAddress().at(1)<<endl;
            LogTrace("l1tMuBayesEventPrint")<<" itOmtfCand->trackAddress().at(2) "<<itOmtfCand->trackAddress().at(2)<<endl;
            LogTrace("l1tMuBayesEventPrint")<<" itOmtfCand->trackAddress().at(3) "<<itOmtfCand->trackAddress().at(3)<<endl;*/

          if(itOmtfCand->hwQual() < muCandQualityCut)
            continue;

          int omtfTTTrackIndex = (int)itOmtfCand->trackAddress().at(3);
          edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, omtfTTTrackIndex);

          if(l1track_ptr.get() == matchedTTTrack.get() ) {
            double ptTTTrackOmtf = l1track_ptr->getMomentum(L1Tk_nPar).perp();
            if(tpPtr->pt() > minMuPt) {
              omtfAndTtMuonEta->Fill(tpPtr->eta());
            }
            if(tpPtr->pt() > 20 && ptTTTrackOmtf >= 18) {
              omtfAndTtMuonEta_ptGen20GeV_ptTT18Gev->Fill(tpPtr->eta());
            }

            if (abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
              omtfAndTtMuonPt->Fill(tpPtr->pt());
              omtfAndTtMuonPhi->Fill(tpPtr->phi());

              ptGenPtOMtfMuon->Fill(tpPtr->pt(), ptTTTrackOmtf);

              if(tpPtr->eventId().event() == 0) {
                ptGenPtOMtfMuonEv0->Fill(tpPtr->pt(), ptTTTrackOmtf);
              }
              else {
                ptGenPtOMtfMuonEvPu->Fill(tpPtr->pt(), ptTTTrackOmtf);
              }

              int layerHits = (int)itOmtfCand->trackAddress().at(0);
              std::bitset<29> layerHitBits(layerHits);
              int firedLayers = layerHitBits.count();
              int pdfSum = (int)itOmtfCand->trackAddress().at(2);
              if(ptTTTrackOmtf >= 18)
                muonsPdfSumFiredPlanes->Fill(pdfSum, firedLayers);

              isOmtf = true;
              //if (DebugMode)
              {
                //int refLayer = (int)it->trackAddress().at(1);
                int omtfTTTrackIndex = (int)itOmtfCand->trackAddress().at(3);
                LogTrace("l1tMuBayesEventPrint")<<" omtf track: omtfTTTrackIndex "<<omtfTTTrackIndex<<" pt "<<itOmtfCand->hwPt()<<" = "<<(itOmtfCand->hwPt()-1)/2.<<" GeV "
                    <<" eta "<<OMTFConfiguration::hwEtaToEta(itOmtfCand->hwEta())
                <<" phi "<<OMTFConfiguration::hwPhiToGlobalPhi(itOmtfCand->hwPhi() )
                <<" hwSign "<<itOmtfCand->hwSign()
                <<" "<<layerHitBits<<endl<<endl;
              }
              break; //can be only one correlator muon for one ttTrack
            }
          }
          /*
           * in principle for a given tracking particle there can be two (or more) matched ttTracks, but then most probably the OMTF will select only one, the second will be ghost busted
           */
        }
      }
    }// end loop over matched L1 tracks

    // ----------------------------------------------------------------------------------------------
    if(!isOmtf && nMatch && abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
      lostTtMuonPt->Fill(tpPtr->pt());
      if(tpPtr->pt() > 20) {
        lostTtMuonEta_ptGen20GeV->Fill(tpPtr->eta());
        lostTtMuonPhi_ptGen20GeV->Fill(tpPtr->phi());
      }

      ptGenPtOMtfMuon->Fill(tpPtr->pt(), 0);
      if(tpPtr->eventId().event() == 0) {
        ptGenPtOMtfMuonEv0->Fill(tpPtr->pt(), 0);
      }
      else {
        ptGenPtOMtfMuonEvPu->Fill(tpPtr->pt(), 0);
      }

      LogTrace("l1tMuBayesEventPrint")<<" no omtf candidate. omtfCands->size(): "<<omtfCands->size()<<(tpPtr->pt() > 10 ? " missing high Pt" : "")<<endl<<endl;
    }
    if (nMatch > 1)
      LogTrace("l1tMuBayesEventPrint") << "WARNING *** 2 or more matches to genuine L1 tracks ***" << endl;

    if (nMatch > 0) {
      float matchTTTrackPt   = bestMatchedTTTrack->getMomentum(L1Tk_nPar).perp();
      float matchTTTrackEta  = bestMatchedTTTrack->getMomentum(L1Tk_nPar).eta();
      float matchTTTrackPhi  = bestMatchedTTTrack->getMomentum(L1Tk_nPar).phi();

      if(tpPtr->pt() > minMuPt) {
        ttTrackMuCnt++;
        ttMuonEta->Fill(tpPtr->eta());
      }
      if(tpPtr->pt() > 20 && matchTTTrackPt >= 18) {
        ttMuonEta_ptGen20GeV_ptTT18Gev->Fill(tpPtr->eta());
      }
      if(abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
        if(tpPtr->pt() > minMuPt) {
          ttTrackMuInOmtfCnt++;
        }
        ttMuonPt->Fill(tpPtr->pt());
        ttMuonPhi->Fill(tpPtr->phi());

        ptGenPtTTMuon->Fill(tpPtr->pt(), matchTTTrackPt);

        if(tpPtr->eventId().event() == 0) {
          ptGenPtTTMuonEv0->Fill(tpPtr->pt(), matchTTTrackPt);
        }
        else {
          ptGenPtTTMuonEvPu->Fill(tpPtr->pt(), matchTTTrackPt);
        }

        ptGenDeltaPtTTMuon->Fill(tpPtr->pt(), tpPtr->pt() - matchTTTrackPt);
        ptGenDeltaPhiTTMuon->Fill(tpPtr->pt(), tpPtr->phi() - matchTTTrackPhi);
        ptGenDeltaEtaTTMuon->Fill(tpPtr->pt(), tpPtr->eta() - matchTTTrackEta);
      }
    }
    else { //no matched ttTrack
      if(abs(tpPtr->eta()) >= etaCutFrom && abs(tpPtr->eta()) <= etaCutTo) {
        ptGenPtTTMuon->Fill(tpPtr->pt(), 0);
        ptGenPtOMtfMuon->Fill(tpPtr->pt(), 0);

        if(tpPtr->eventId().event() == 0) {
          ptGenPtTTMuonEv0->Fill(tpPtr->pt(), 0);
          ptGenPtOMtfMuonEv0->Fill(tpPtr->pt(), 0);
        }
        else {
          ptGenPtTTMuonEvPu->Fill(tpPtr->pt(), 0);
          ptGenPtOMtfMuonEvPu->Fill(tpPtr->pt(), 0);
        }

        double tpMinPt = 5;
        if (tpPtr->pt() > tpMinPt) {// && DebugMode
          LogTrace("l1tMuBayesEventPrint") << "\nno ttTrack matched to tracking particle:";
          LogTrace("l1tMuBayesEventPrint")<<printTrackigParticle();
          LogTrace("l1tMuBayesEventPrint")<<"Printing all ttTracks around tracking particle" << endl;
        }

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
              LogTrace("l1tMuBayesEventPrint") << "L1 track, pt: " << ttTrkPt //<<" RInv "<<l1track_ptr->getRInv()
                  << " eta: " << ttTrkEta << " phi: " << ttTrkPhi
                  << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub<<", "
                  <<(tmp_trk_genuine ? "genuine, " : "")
                  <<(tmp_trk_loose ? "loose, " : "")
                  <<(tmp_trk_unknown ? "unknown, " : "")
                  <<(tmp_trk_combinatoric ? "combinatoric, " : "");

              if( !my_tp.isNull() ) {
                LogTrace("l1tMuBayesEventPrint") <<" is matched to tracking particle: ";
                LogTrace("l1tMuBayesEventPrint")<< "Tracking particle,    pt: " <<setw(7)<< my_tp->pt() << " eta: " <<setw(7)<< my_tp->eta() << " phi: " <<setw(7)<< my_tp->phi()
                          << " pdgid: " <<setw(7)<< my_tp->pdgId() << " eventID: " <<setw(7)<< my_tp->eventId().event()
                          << " ttTracks Cnt " << MCTruthTTTrackHandle->findTTTrackPtrs(my_tp).size()
                          <<" bx: "<<my_tp->eventId().bunchCrossing()<<" address "<<my_tp.get()
                          <<" key "<<my_tp.key()
                          <<" "<<(my_tp == tpPtr ? " my_tp == tp_ptr " : "my_tp != tp_ptr")<< endl;
              }
              else {
                LogTrace("l1tMuBayesEventPrint")<<" no matching tracking particle"<<endl;
              }
              LogTrace("l1tMuBayesEventPrint")<<"";
            }
          }
        }
      }
    }
  } //end loop tracking particles

  gpMuonPerEvent->Fill(allGpMuCnt);
  ttTrackMuonPerEvent->Fill(ttTrackMuCnt);
  ttTrackMuonInOmtfRegPerEvent->Fill(ttTrackMuInOmtfCnt);

  /*
   * in principle only running with LooseMatch= True has sense, otherwise how to count the ttTracks that are loose genuine matched to muon tracking particle
   * and are tagged as muon by omtf?
   */
  //TODO in principle one should find one candidate with max pt and use it for filling the histograms
  l1t::RegionalMuonCandBxCollection::const_iterator itBestOmtfCand = l1Omtf.product()->end(bxNumber);
  float bestPt = 0;
  for(l1t::RegionalMuonCandBxCollection::const_iterator itOmtfCand = omtfCands->begin(bxNumber);
        itOmtfCand != l1Omtf.product()->end(bxNumber); ++itOmtfCand) {

    int omtfTTTrackIndex = (int)itOmtfCand->trackAddress().at(3);
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ttTrackPtr(TTTrackHandle, omtfTTTrackIndex);

    int ttTrackNstub = ttTrackPtr->getStubRefs().size();
    if (ttTrackNstub < L1Tk_minNStub)
      continue;

    int quality = itBestOmtfCand->hwQual();
    if(quality < muCandQualityCut)
      continue;
    /*
     * so only the omtf ttTRacks with L1Tk_minNStub or more stubs are counted further, including counting as fakes
     */

    if(bestPt < ttTrackPtr->getMomentum(L1Tk_nPar).perp() ) {
      bestPt = ttTrackPtr->getMomentum(L1Tk_nPar).perp();
      itBestOmtfCand = itOmtfCand;
    }
  }


  if(itBestOmtfCand != l1Omtf.product()->end(bxNumber))
  {
    int layerHits = (int)itBestOmtfCand->trackAddress().at(0);
    //int refLayer = (int)it->trackAddress().at(1);
    int omtfTTTrackIndex = (int)itBestOmtfCand->trackAddress().at(3);

    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ttTrackPtr(TTTrackHandle, omtfTTTrackIndex);

    int ttTrackNstub = ttTrackPtr->getStubRefs().size();

    //std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
    float ttTrkPt   = ttTrackPtr->getMomentum(L1Tk_nPar).perp();
    float ttTrkEta  = ttTrackPtr->getMomentum(L1Tk_nPar).eta();
    float ttTrkPhi  = ttTrackPtr->getMomentum(L1Tk_nPar).phi();

    omtfTTTrackPt->Fill(ttTrkPt); //filling the data from the l1track_ptr, but this is principle should be the same as in the omtfMuon
    omtfTTTrackEta->Fill(ttTrkEta);
    omtfTTTrackPhi->Fill(ttTrkPhi);

       //cout << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__ << endl;

    /*int tmp_trk_genuine = 0;
    int tmp_trk_loose = 0;
    int tmp_trk_unknown = 0;
    int tmp_trk_combinatoric = 0;
    if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
    if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
    if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
    if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;
*/

    edm::Ptr< TrackingParticle > tpMatchedToOmtfTTTrack = MCTruthTTTrackHandle->findTrackingParticlePtr(ttTrackPtr);
    if(tpMatchedToOmtfTTTrack.isNonnull() && abs(tpMatchedToOmtfTTTrack->pdgId()) == 13) {
      omtfTTTrackMuonPt->Fill(ttTrkPt);
      omtfTTTrackMuonEta->Fill(ttTrkEta);
      omtfTTTrackMuonPhi->Fill(ttTrkPhi);
    }

    float chi2dof = (ttTrackPtr->getChi2(L1Tk_nPar)) / (2*ttTrackNstub - L1Tk_nPar);
    if(tpMatchedToOmtfTTTrack.isNull()) {
      chi2FakeTTTrackOmtf->Fill(ttTrackPtr->getChi2(L1Tk_nPar) );
      chi2DofFakeTTTrackOmtf->Fill(chi2dof);
    }
    else {
      chi2GenuineTTTrackOmtf->Fill(ttTrackPtr->getChi2(L1Tk_nPar) );
      chi2DofGenuineTTTrackOmtf->Fill(chi2dof);
    }

    //tpMatchedToOmtfTTTrack.isNull() means that the ttTRack is not loose genuine nor genuine
    if(tpMatchedToOmtfTTTrack.isNull() || abs(tpMatchedToOmtfTTTrack->pdgId()) != 13) {
      std::bitset<29> layerHitBits(layerHits);
      int firedLayers = layerHitBits.count();
      int pdfSum = (int)itBestOmtfCand->trackAddress().at(2);

      //ptGenPtTTMuonEvPu->Fill(0., ttTrkPt);
      ptGenPtOMtfMuonEvPu->Fill(0., ttTrkPt);

      if(tpMatchedToOmtfTTTrack.isNull() ) {
        fakeTrackOmtfPt->Fill(ttTrkPt);
        if(ttTrkPt >= 18) {
          fakeTrackOmtfEta_ptTT18Gev->Fill(ttTrkEta);
          fakeTrackOmtfPhi_ptTT18Gev->Fill(ttTrkPhi);
        }
      }
      else if(abs(tpMatchedToOmtfTTTrack->pdgId()) != 13) {
        wrongTagOmtfPt->Fill(ttTrkPt);
        if(ttTrkPt >= 18) {
          wrongTagOmtfEta_ptTT18Gev->Fill(ttTrkEta);
          wrongTagOmtfPhi_ptTT18Gev->Fill(ttTrkPhi);
        }
      }

      if(ttTrkPt >= 18)
        notMuonsPdfSumFiredPlanes->Fill(pdfSum, firedLayers);

      if (ttTrkPt > 10.) //DebugMode ||
      {
      	edm::LogImportant("l1tMuBayesEventPrint")<<"\nrun:lumi:event "<<event.run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
        edm::LogImportant("l1tMuBayesEventPrint")<<"fake omtf track: omtfTTTrackIndex "<<omtfTTTrackIndex<<" pt "<<itBestOmtfCand->hwPt()<<" = "<<(itBestOmtfCand->hwPt()-1)/2.<<" GeV "
            <<" eta "<<OMTFConfiguration::hwEtaToEta(itBestOmtfCand->hwEta())
        <<" phi "<<OMTFConfiguration::hwPhiToGlobalPhi(itBestOmtfCand->hwPhi() )
        <<" hwSign "<<itBestOmtfCand->hwSign()
        <<" pdfSum "<<pdfSum<<" firedLayers "<<firedLayers
        <<" "<<layerHitBits<<(itBestOmtfCand->hwPt() > 21 ? " fake high Pt" : "")<<endl;

        bool tmp_trk_genuine = MCTruthTTTrackHandle->isGenuine(ttTrackPtr);
        bool tmp_trk_loosegenuine = MCTruthTTTrackHandle->isLooselyGenuine(ttTrackPtr);
        edm::LogImportant("l1tMuBayesEventPrint") << printTTTRack(ttTrackPtr, tmp_trk_genuine, tmp_trk_loosegenuine);

        if(!tpMatchedToOmtfTTTrack.isNull() ) {
          edm::LogImportant("l1tMuBayesEventPrint") << "TP matched to omtf:";
          edm::LogImportant("l1tMuBayesEventPrint") << printTrackigParticleShort(tpMatchedToOmtfTTTrack);
        }
        else {
          bool unknown  = MCTruthTTTrackHandle->isUnknown(ttTrackPtr);
          bool combinatoric  = MCTruthTTTrackHandle->isCombinatoric(ttTrackPtr);
          edm::LogImportant("l1tMuBayesEventPrint") << " no matched track particle, the ttTrack is fake: "
               << (unknown ? " unknown " : "")
               << (combinatoric ? " combinatoric " : "")
               <<endl;
        }
        edm::LogImportant("l1tMuBayesEventPrint")<<"";
      }
    }
  }

}




DEFINE_FWK_MODULE(MuCorrelatorAnalyzer);

