#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include <string>


using namespace std;

//object definition
class OmtfAnalyzer : public edm::EDAnalyzer {
public:

  //constructor, function is called when new object is created
  explicit OmtfAnalyzer(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~OmtfAnalyzer();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  double omtfPtCut1 = 20;
  double simPtCut1 = 25;
  double simPtCut2 = omtfPtCut1;
  double simPtCut3 = 10;

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  //added variables
  TFile* myRootFile;
  TH1D* hPhiOmtf;
  TH2D* hPtSimPt;
  TH1D* hPtSim;
  TH1D* hPtOmtf;
  TH1D* hEtaSim;
  TH1D* hEtaOmtf;
  TH1D* hEtaOmtfSim;
  TH1D* hEtaOmtfSimCut;
  TH1D* hPhiSim;
  TH1D* hRefLayer;
  TH1D* hRefLayerCut;
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

  vector<TH2D*> hPtSimPtVsRefLay;

  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> inputOMTF;
  edm::EDGetTokenT<edm::SimTrackContainer> inputSim;
};


OmtfAnalyzer::OmtfAnalyzer(const edm::ParameterSet& conf) 
  : theConfig(conf), theEventCount(0) 
{
  cout <<" CTORXX" << endl;
  inputOMTF = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simOmtfDigis", "OMTF")); //simOmtfDigis
  inputSim =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
}


OmtfAnalyzer::~OmtfAnalyzer() 
{ 
  cout <<" DTOR" << endl;
}

void OmtfAnalyzer::beginJob()
{
  //make a new Root file

  string  outRootFile = theConfig.getParameter<std::string>("outRootFile");
  myRootFile = new TFile(outRootFile.c_str(),"RECREATE");
  if(myRootFile == 0 || myRootFile-> IsZombie() ) {
    cout<<"file "<<outRootFile<<" not opened"<<endl;
    throw cms::Exception("OMTFConfiguration::getPatternPtRange: patternPts vector not initialized");
  }

  cout<<"file created "<<outRootFile<<endl;
  const int nIpt = 27;
  double lower[nIpt + 1] = {0, 4., 4.5, 5., 6., 7., 8., 10., 12., 14., 16., 18.,
      20., 22. , 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120., 140., 500.};
//  for (int i = 0; i <= nIpt; i++){
//    lower[i] = lower[i] - 0.25;
//  }
  const int bins = 1000;
  //create histograms

  string omtfPtCut1Str= std::to_string((int)omtfPtCut1);
  string simPtCut1Str = std::to_string((int)simPtCut1);
  string simPtCut2Str = std::to_string((int)simPtCut2);
  string simPtCut3Str = std::to_string((int)simPtCut3);

  hPhiOmtf =new TH1D("hPhiOmtf","test; phi omtf; #events",360, -2*M_PI, 2*M_PI);
  hPtSimPt = new TH2D("hPtSimPt", "pt_sim - pt_omtf; pt_sim; pt_omtf", bins, 0., 500., nIpt, lower);
  hPtSim = new TH1D("pt_sim", "pt_sim; pt_sim; #events", bins, 0., 500.);
  hPtOmtf = new TH1D("pt", "pt_omtf; pt_omtf; #events", bins, 0., 500.);
  hEtaSim = new TH1D("eta_sim", "eta_sim; eta_sim; #events", bins, -1., 3.);
  hEtaOmtf = new TH1D("eta_omtf", "eta_omtf; eta_omtf; #events", bins, -1., 3.);
  hEtaOmtfSim = new TH1D("eta_omtfsim", "EtaSim, if OMTF fired; eta_sim; #events", bins, -1., 3.);
  hEtaOmtfSimCut = new TH1D("eta_omtfsimcut", ("EtaSim, if ptSim > " + simPtCut1Str + " GeV and OMTF fired; eta_sim; #events").c_str(), bins, -1., 3.);
  hPhiSim = new TH1D("phi_sim", "phi_sim; phi_sim; #events", bins, 0., 2*M_PI);
  hRefLayer = new TH1D("refLayer", "refLayer; refLayer; #events", 10, -0.5, 10-0.5);

  hRefLayerCut = new TH1D("refLayerCut", ("RefLayer pt_omtf >= " + omtfPtCut1Str + " GeV; refLayer; #events").c_str(), 10, -0.5, 10-0.5);
  hRefLayerCut2 = new TH1D("refLayerCut2", ("RefLayer pt_omtf >= " + omtfPtCut1Str + ", sim < "  + simPtCut3Str + " GeV; refLayer; #events").c_str(), 10, -0.5, 10-0.5);

  hLayerHits = new TH1D("LayerHits", "LayerHits; Layer; #events", 18, -0.5, 18-0.5);
  hLayerHitsCut = new TH1D("LayerHitsCut", ("LayerHitsCut pt_omtf >= " + omtfPtCut1Str + " GeV; Layer; #events").c_str(), 18, -0.5, 18-0.5);
  hLayerHitsCut2 = new TH1D("LayerHitsCut2", ("LayerHitsCut2 pt_omtf >= " + omtfPtCut1Str + ", sim < " + simPtCut3Str + " GeV; Layer; #events").c_str(), 18, -0.5, 18-0.5);

  hRefL1Layers = new TH1D("RefL1Layers", ("RefLayer 1, fired layers ptOmtf >= " + omtfPtCut1Str + " GeV, sim < " + simPtCut3Str + " GeV; Layer; #events").c_str(), 18, -0.5, 18-0.5);
  hRefL3Layers = new TH1D("RefL3Layers", ("RefLayer 3, fired layers ptOmtf >= " + omtfPtCut1Str + " GeV, sim < " + simPtCut3Str + " GeV; Layer; #events").c_str(), 18, -0.5, 18-0.5);

  hNLayers = new TH1D("NLayers;", ("Fired layers count; pt_omtf >= " + omtfPtCut1Str + " GeV, sim < " + simPtCut3Str + " GeV; Fired layers; #events").c_str(), 18, -0.5, 18-0.5);
  hQuality = new TH1D("Quality", "Quality; quality; #events", 32, -0.5, 32-0.5);

  hEtaSimCut =   new TH1D("hEtaSimCut", ("eta_sim if ptSim > " + simPtCut1Str + " GeV; eta_sim; #events").c_str(), bins, -1., 3.);
  hEffVsEtaSim = new TH1D("hEffVsEtaSim", ("Efficiency vs EtaSim, ptSim > " + simPtCut1Str + " GeV; eta_sim; efficiency").c_str(), bins, -1., 3.);
  hEffVsEtaSimOmtfCut = new TH1D("hEffVsEtaSimOmtfCut", ("Efficiency vs EtaSim, ptOmtf > " + omtfPtCut1Str + " GeV ptSim > " + simPtCut1Str + " GeV; eta_sim; efficiency").c_str(), bins, -1., 3.);

  cout << "HERE Cwiczenie::beginJob()" << endl;

  for(int i = 0; i < 8; i++) //in the 8 the events without omtf cand are going
  {
    hPtSimPtVsRefLay.push_back(new TH2D(Form("hPtSimPt_ref%i", i), "test; pt_sim; pt_omtf", bins, 0., 500., nIpt, lower));
  }
  cout << "HERE Cwiczenie::beginJob()" << endl;
}

void OmtfAnalyzer::endJob()
{
  cout << "HERE Cwiczenie::endJob():"<<__LINE__ << endl;
  myRootFile->cd();
  //write histogram data
  hPhiOmtf->Write();
  cout << "HERE Cwiczenie::endJob():"<<__LINE__ << endl;
  hPtSimPt->Write();
  hPtSim->Write();
  hPtOmtf->Write();
  hEtaSim->Write();
  hEtaOmtf->Write();
  hEtaOmtfSim->Write();
  hEtaOmtfSimCut->Write();
  hPhiSim->Write();
  hRefLayer->Write();
  hRefLayerCut->Write();
  hRefLayerCut2->Write();
  hLayerHits->Write();
  hLayerHitsCut->Write();
  hLayerHitsCut2->Write();
  hQuality->Write();
  hRefL1Layers->Write();
  hRefL3Layers->Write();
  hNLayers->Write();

  hEffVsEtaSim->Divide(hEtaOmtfSimCut, hEtaSimCut);
  hEffVsEtaSim->Write();

  hEffVsEtaSimOmtfCut->Divide(hEtaSimCut);
  hEffVsEtaSimOmtfCut->Write();

  hEtaSimCut->Write();

  for(unsigned int i = 0; i < hPtSimPtVsRefLay.size(); i++)
  {
    hPtSimPtVsRefLay.at(i)->Write();
    //delete hPtSimPt.at(i);
  }

  myRootFile->Close();
/*  delete histo;
  delete hPtSimPt;
  delete myRootFile;
  delete hPt;
  delete hPtSim;
  delete hEtaSim;
  delete hEtaOmtf;
  delete hEtaOmtfSim;
  delete hEtaOmtfSimCut;
  delete hPhiSim;
  delete hRefLayer;
  delete hRefLayerCut;
  delete hRefLayerCut2;
  delete hLayerHits;
  delete hLayerHitsCut;
  delete hLayerHitsCut2;
  delete hQuality;
  delete hRefL1Layers;
  delete hRefL3Layers;
  delete hNLayers;*/
  cout << "HERE Cwiczenie::endJob()" << endl;
}

void OmtfAnalyzer::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  //std::cout << " HERE Cwiczenie::analyze "<< std::endl;

  float pt_sim = 0, phi_sim = 0, eta_sim = 0;
  float  phi_omtf = -1, eta_omtf = 0, pt_omtf_max = 0;
  std::vector<float> pt_omtf;
  int refLayer = 9;
  int layerHit = 0;
  unsigned int simMuonCount=0; 
  unsigned int omtfMuonCount=0;
  int quality = 0;

  //std::cout <<" SIMULATED MUONS: "<<std::endl;
  edm::Handle<edm::SimTrackContainer> simTk;
  ev.getByToken(inputSim, simTk);
  std::vector<SimTrack> mySimTracks = *(simTk.product());
  for (std::vector<SimTrack>::const_iterator it=mySimTracks.begin(); it<mySimTracks.end(); it++) {
    const SimTrack & track = *it;
    if ( track.type() == -99) continue;
    if ( track.vertIndex() != 0) continue;

    //sucessful muon, add to count
    simMuonCount++;

    phi_sim = track.momentum().phi(); //momentum azimutal angle
    pt_sim = track.momentum().pt(); //transverse momentum
    eta_sim = track.momentum().eta(); //pseudorapidity
    //std::cout <<" trackId: " <<track.trackId() 
     //     << " pt_sim: " << pt_sim <<" eta_sim: "<<eta_sim<<" phi_sim: "<<phi_sim
     //     <<" vtx: "<<track.vertIndex()<<" type: "<<track.type() // 13 or -13 is a muon
     //     << std::endl; 
  }
  if(simMuonCount != 1) {
     //cout<<"    Simulated muon count != 1"<<endl;
     return;
  }
  hPtSim->Fill(pt_sim);
  hEtaSim->Fill(eta_sim);
  hPhiSim->Fill(phi_sim);
  if(pt_sim > simPtCut1) {
    hEtaSimCut->Fill(eta_sim);
  }

  //std::cout <<" L1 MUONS: "<<std::endl;
  edm::Handle<l1t::RegionalMuonCandBxCollection> l1Omtf;
  ev.getByToken(inputOMTF, l1Omtf);
  int bxNumber = 0;
  l1t::RegionalMuonCandBxCollection::const_iterator itBestOmtfCand;

  for(l1t::RegionalMuonCandBxCollection::const_iterator it = l1Omtf.product()->begin(bxNumber);
       it != l1Omtf.product()->end(bxNumber); ++it) {
    omtfMuonCount++;

    quality = it->hwQual();
    //QUALITY CUT
    if(quality > 8) {
      if((it->hwPt()-1.)/2. > pt_omtf_max)
      {
        itBestOmtfCand = it;
        pt_omtf_max = (it->hwPt()-1.)/2.;
      }
    }
  }

  if(pt_omtf_max)
  {
    refLayer = (int)itBestOmtfCand->trackAddress().at(1);
    layerHit = (int)itBestOmtfCand->trackAddress().at(0);

    phi_omtf = ( (15.+itBestOmtfCand->processor()*60.)/360. + itBestOmtfCand->hwPhi()/576. ) *2*M_PI;
    if (phi_omtf > 2*M_PI)
      phi_omtf -=  2*M_PI;

    eta_omtf = itBestOmtfCand->hwEta()/240.*2.26;
  }

  if(simMuonCount==1) {//does not matter, return before
    if(eta_sim > 0.82 && eta_sim < 1.24)
    {
      hPtSimPt->Fill(pt_sim, pt_omtf_max);
      hPtOmtf->Fill(pt_omtf_max);
    }
  }

  //all values used for hist filling are for the itBestOmtfCand
  if(omtfMuonCount >= 1 &&eta_sim > 0.82 && eta_sim < 1.24) {
    hQuality->Fill(quality);

    if(pt_omtf_max) { //after quality cut
      hPhiOmtf->Fill(phi_omtf);
      hEtaOmtf->Fill(eta_omtf);
      hEtaOmtfSim->Fill(eta_sim);
      hRefLayer->Fill(refLayer);

      hPtSimPtVsRefLay.at(refLayer)->Fill(pt_sim, pt_omtf_max);

      for(int i = 0; i < 18; i ++) {
        if((layerHit >> i)%2 == 1) {
          hLayerHits->Fill(i);
        }
      }

    }

    if(pt_sim > simPtCut1)
      hEtaOmtfSimCut->Fill(eta_sim);

    if(pt_omtf_max >= omtfPtCut1) {
      if(pt_sim > simPtCut1) {
        hEffVsEtaSimOmtfCut->Fill(eta_sim);
      }

      for(int i = 0; i < 18; i ++) {
        if((layerHit >> i)%2 == 1)
          hLayerHitsCut->Fill(i);
      }

      hRefLayerCut->Fill(refLayer);
      if (pt_sim < simPtCut3) {
        hRefLayerCut2->Fill(refLayer);
        int nLayers = 0;
        for(int i = 0; i < 18; i ++) {
          if((layerHit >> i)%2 == 1) {
            hLayerHitsCut2->Fill(i);
            nLayers++;
          }
        }
        hNLayers->Fill(nLayers);

        //RefLayer 3 HIT
        if((layerHit >> 6)%2 == 1) {
          for(int i = 0; i < 18; i ++)
          {
            if((layerHit >> i)%2 == 1)
              hRefL3Layers->Fill(i);
          }
        }

        //RefLayer 1 HIT
        if((layerHit >> 7)%2 == 1) {
          for(int i = 0; i < 18; i ++) {
            if((layerHit >> i)%2 == 1)
              hRefL1Layers->Fill(i);
          }
        }
      }
    }
  }
  
  //write std io
  //cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" useful event count:"<<++theEventCount << endl;
}




DEFINE_FWK_MODULE(OmtfAnalyzer);

