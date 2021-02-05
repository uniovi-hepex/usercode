/*
 * L1MuonAnalyzerOmtf.cc
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "usercode/L1MuonAnalyzer/plugins/L1MuonAnalyzerOmtf.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

namespace L1MuAn {


L1MuonAnalyzerOmtf::L1MuonAnalyzerOmtf(const edm::ParameterSet& edmCfg): muonMatcher(edmCfg) {
  fillMatcherHists = !edmCfg.exists("muonMatcherFile");
  edm::LogImportant("l1tOmtfEventPrint") <<" L1MuonAnalyzerOmtf: line "<<__LINE__<<" fillMatcherHists "<<fillMatcherHists<<std::endl;

  if(edmCfg.exists("matchUsingPropagation") )
    matchUsingPropagation = edmCfg.getParameter<bool>("matchUsingPropagation");

  edm::LogImportant("l1tOmtfEventPrint") <<" L1MuonAnalyzerOmtf: line "<<__LINE__<<" matchUsingPropagation "<<matchUsingPropagation<<std::endl;

  omtfToken = consumes<l1t::RegionalMuonCandBxCollection >(edmCfg.getParameter<edm::InputTag>("L1OMTFInputTag"));

  if(edmCfg.exists("simVertexesTag") ) {
    simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //TODO which is correct?

    simVertexesToken =  consumes<edm::SimVertexContainer>(edmCfg.getParameter<edm::InputTag>("simVertexesTag"));
  }
  if(edmCfg.exists("trackingParticleTag") )
    trackingParticleToken = consumes<TrackingParticleCollection>(edmCfg.getParameter<edm::InputTag>("trackingParticleTag"));

  edm::Service<TFileService> fs;

  candPerEvent = fs->make<TH1D>("candPerEvent", "candPerEvent", 21, -0.5, 20.5);

  analysisType = edmCfg.getParameter< string >("analysisType");

  int binsCnt = (1<<18)-1;
  firedLayersEventCntOmtf = fs->make<TH1S>("firedLayersEventCntOmtf", "firedLayersEventCntOmtf", binsCnt, 0, binsCnt);
  firedLayersEventCntNN = fs->make<TH1S>("firedLayersEventCntNN", "firedLayersEventCntNN", binsCnt, 0, binsCnt);

  if(edmCfg.exists("nn_pThresholds") )
    nn_pThresholds = edmCfg.getParameter<vector<double> >("nn_pThresholds");

  if(analysisType == "efficiency") {
    ptGenHist = fs->make<TH1I>("ptGenHist", "ptGenHist", 800, 0, 4000);

    double ptGenCut = 25;
    double ptL1Cut = 18;
    TFileDirectory subDir = fs->mkdir("efficiency");
    auto addOmtfAnalysers = [&](TFileDirectory& subDir, std::string name, int qualityCut) {
      omtfEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, name, 0.82, 1.24, qualityCut, 800, 0, 800));
      omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, name, 0.82, 1.24, qualityCut, ptGenCut, ptL1Cut, 100));
      omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, name, qualityCut, ptGenCut, ptL1Cut, 210));
      omtfEfficiencyAnalysers.emplace_back(new LikelihoodDistribution(subDir, name, qualityCut, ptGenCut, ptL1Cut, 120));
    };

    addOmtfAnalysers(subDir, "omtf_q12", 12);
    addOmtfAnalysers(subDir, "omtf_q8", 8);
    //addOmtfAnalysers(subDir, "omtf_q4", 4);
    //addOmtfAnalysers(subDir, "omtf_q1", 1);

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"nn_omtf_q1_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 1, 200, 0, 200));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 1, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 1, ptGenCut, ptL1Cut, 210));

      ostr.str("");
      ostr<<"nn_omtf_q4_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 4, 200, 0, 200));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 4, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 4, ptGenCut, ptL1Cut, 210));

      ostr.str("");
      ostr<<"nn_omtf_q12_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 12, 200, 0, 200));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 12, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 12, ptGenCut, ptL1Cut, 210));
    }
  }

  else if(analysisType == "rate") {
    TFileDirectory subDirRate = fs->mkdir("rate");

    TFileDirectory subDir = subDirRate.mkdir("omtf_q12");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 12, 200, 0, 100));
    omtfCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 12, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q8");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 8, 200, 0, 100));
    omtfCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 8, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q4");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 4, 200, 0, 100));
    omtfCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 4, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q1");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 1, 200, 0, 100));
    omtfCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 1, 200, 0, 100));

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"nn_omtf_q12_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 12, 200, 0, 100));
      omtfNNCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 12, 200, 0, 100));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"nn_omtf_q4_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 4, 200, 0, 100));
      omtfNNCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 4, 200, 0, 100));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"nn_omtf_q1_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 1, 200, 0, 100));
      omtfNNCandsMatchingAnalysers.emplace_back(new CandsMatchingAnalyser(subDir, "", 1, 200, 0, 100));
      edm::LogImportant("l1tOmtfEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
    }
  }
}


L1MuonAnalyzerOmtf::~L1MuonAnalyzerOmtf() {
  // TODO Auto-generated destructor stub
}

void L1MuonAnalyzerOmtf::beginJob() {

}

std::vector<const l1t::RegionalMuonCand*> L1MuonAnalyzerOmtf::ghostBust(const l1t::RegionalMuonCandBxCollection* mtfCands) {
  boost::dynamic_bitset<> isKilled(mtfCands->size(0), false);

  for(unsigned int  i1 = 0; i1 < mtfCands->size(0); ++i1) {
    for(unsigned int  i2 = i1 +1; i2 < mtfCands->size(0); ++i2) {
      auto& mtfCand1 = mtfCands->at( 0, i1 );
      auto& mtfCand2 = mtfCands->at( 0, i2 );

      if( abs( mtfCand1.hwEta() - mtfCand2.hwEta()) < (0.3 / 0.010875) ) {
        int gloablHwPhi1 = l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCand1.hwPhi(), mtfCand1.trackFinderType(), mtfCand1.processor() ) ;
        int gloablHwPhi2 = l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCand2.hwPhi(), mtfCand2.trackFinderType(), mtfCand2.processor() ) ;
        //double globalPhi1 = hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCand1.hwPhi(), mtfCand1.trackFinderType(), mtfCand1.processor() ) );
        //double globalPhi2 = hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCand2.hwPhi(), mtfCand2.trackFinderType(), mtfCand2.processor() ) );

        if(abs(gloablHwPhi1 - gloablHwPhi2) < 8) {//0.0872664626 = 5 deg, i.e. the same window as in the OMTF ghost buster
          if(mtfCand1.hwQual() > mtfCand2.hwQual()) {
            isKilled[i2] = true;
          }
          else
            isKilled[i1] = true;
        }
      }
    }
  }

  std::vector<const l1t::RegionalMuonCand*> resultCands;

  for(unsigned int  i1 = 0; i1 < mtfCands->size(0); ++i1) {
    if(!isKilled[i1] && mtfCands->at( 0, i1 ).hwQual()) { //dropping candidates with quality 0 !!!!!!!!!!!!!!!!!!!! fixme if not needed
      resultCands.push_back(&(mtfCands->at( 0, i1 ) ));
    }

    LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::ghostBust pt "<<std::setw(3)<<mtfCands->at( 0, i1 ).hwPt()<<" qual "<<std::setw(2)<<mtfCands->at( 0, i1 ).hwQual()
        <<" proc "<<std::setw(2)<<mtfCands->at( 0, i1 ).processor()<<" eta "<<std::setw(4)<<mtfCands->at( 0, i1 ).hwEta()<<" gloablEta "<<std::setw(8)<<mtfCands->at( 0, i1 ).hwEta() * 0.010875
        <<" hwPhi "<<std::setw(3)<<mtfCands->at( 0, i1 ).hwPhi()
        <<" globalPhi "<<std::setw(8)<<hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCands->at( 0, i1 ).hwPhi(), mtfCands->at( 0, i1 ).trackFinderType(), mtfCands->at( 0, i1 ).processor() ) )
        <<" fireadLayers "<<std::bitset<18>(mtfCands->at( 0, i1 ).trackAddress().at(0) )
        <<" isKilled "<<isKilled.test(i1)<<std::endl;
  }

  if(resultCands.size() >= 3)
    LogTrace("l1tOmtfEventPrint")<<" ghost !!!!!! " <<std::endl;
  LogTrace("l1tOmtfEventPrint") <<std::endl;

  return resultCands;
}


bool simTrackIsMuonInOmtf(const SimTrack& simTrack) {
  if (abs(simTrack.type()) == 13 || abs(simTrack.type()) == 1000015 ) {  //|| tpPtr->pt() > 20 //todo 1000015 is stau
    //only muons
  }
  else
    return false;

  if(simTrack.momentum().pt() < 2.5) //in the overlap, the propagation of muons with pt less then ~3.2 fails - the actual threshold depends slightly on eta,
    return false;

  LogTrace("l1tOmtfEventPrint") <<"simTrackIsMuonInOmtf, simTrack type "<<std::setw(3)<<simTrack.type()<<" pt "<<std::setw(9)<<simTrack.momentum().pt()<<" eta "<<std::setw(9)<<simTrack.momentum().eta()<<" phi "<<std::setw(9)<<simTrack.momentum().phi()<<std::endl;


  //if( (fabs(simTrack.momentum().eta()) >= 0.82 ) && (fabs(simTrack.momentum().eta()) <= 1.24) ) {
  if( (fabs(simTrack.momentum().eta()) >= 0.72 ) && (fabs(simTrack.momentum().eta()) <= 1.3) ) { //higher margin for matching, otherwise many candidates are marked as ghosts
  }
  else
    return false;;

  /*
  if(tpPtr->vertex().rho() > TP_maxRho)
    continue;
   */

  return true;
}

bool simTrackIsMuonInOmtfBx0(const SimTrack& simTrack) {
  if(simTrack.eventId().bunchCrossing() != 0)
    return false;

  return simTrackIsMuonInOmtf(simTrack);
}

bool trackingParticleIsMuonInOmtfBx0(const TrackingParticle& trackingParticle) {
  //if(trackingParticle.eventId().event() != 0)
/*    LogTrace("l1tOmtfEventPrint") <<"trackingParticleIsMuonInOmtfBx0, pdgId "<<std::setw(3)<<trackingParticle.pdgId()<<" pt "<<std::setw(9)<<trackingParticle.pt()
              <<" eta "<<std::setw(9)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(9)<<trackingParticle.momentum().phi()<<" event "<<trackingParticle.eventId().event()
              <<" bx "<<trackingParticle.eventId().bunchCrossing()<<" eventNot0"<<std::endl;*/

  if(trackingParticle.eventId().bunchCrossing() != 0)
    return false;

  if (abs(trackingParticle.pdgId()) == 13  || abs(trackingParticle.pdgId()) == 1000015 ) {  // || abs(simTrack.pdgId()) == 1000015 || tpPtr->pt() > 20 //todo 1000015 is stau
    //only muons
  }
  else
    return false;

  if(trackingParticle.pt() < 2.5) //in the overlap, the propagation of muons with pt less then ~3.2 fails - the actual threshold depends slightly on eta,
    return false;

  if(trackingParticle.parentVertex().isNonnull() )
    LogTrace("l1tOmtfEventPrint") <<"trackingParticleIsMuonInOmtfBx0, pdgId "<<std::setw(3)<<trackingParticle.pdgId()<<" pt "<<std::setw(9)<<trackingParticle.pt()
          <<" eta "<<std::setw(9)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(9)<<trackingParticle.momentum().phi()<<" event "<<trackingParticle.eventId().event()
          <<" parentVertex Rho "<<trackingParticle.parentVertex()->position().Rho()<<" eta "<<trackingParticle.parentVertex()->position().eta()<<" phi "<<trackingParticle.parentVertex()->position().phi()<<std::endl;
  else
    LogTrace("l1tOmtfEventPrint") <<"trackingParticleIsMuonInOmtfBx0, pdgId "<<std::setw(3)<<trackingParticle.pdgId()<<" pt "<<std::setw(9)<<trackingParticle.pt()
        <<" eta "<<std::setw(9)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(9)<<trackingParticle.momentum().phi();


  //if( (fabs(simTrack.momentum().eta()) >= 0.82 ) && (fabs(trackingParticle.momentum().eta()) <= 1.24) ) {
  if( (fabs(trackingParticle.momentum().eta()) >= 0.72 ) && (fabs(trackingParticle.momentum().eta()) <= 1.3) ) { //higher margin for matching, otherwise many candidates re marked as ghosts
  }
  else
    return false;;



  /*
  if(tpPtr->vertex().rho() > TP_maxRho)
    continue;
   */

  return true;
}

bool trackingParticleIsMuonInOmtfEvent0(const TrackingParticle& trackingParticle) {
  if(trackingParticle.eventId().event() != 0)
    return false;

  return trackingParticleIsMuonInOmtfBx0(trackingParticle);
}

void L1MuonAnalyzerOmtf::analyze(const edm::Event& event, const edm::EventSetup& es) {

  edm::Handle<l1t::RegionalMuonCandBxCollection> l1omtfHandle;
  event.getByToken(omtfToken, l1omtfHandle);

  LogTrace("l1tOmtfEventPrint")<<std::endl<<"\nL1MuonAnalyzerOmtf::analyze. event "<<event.id().event() <<" number of candidates "<<l1omtfHandle.product()->size(0)<<" ---------------------------"<<std::endl;

  std::vector<const l1t::RegionalMuonCand*> ghostBustedCands = ghostBust(l1omtfHandle.product());

  edm::Handle<edm::SimTrackContainer> simTraksHandle;
  edm::Handle<edm::SimVertexContainer> simVertices;
  if(!simTrackToken.isUninitialized()) {
    event.getByToken(simTrackToken, simTraksHandle);
    event.getByToken(simVertexesToken, simVertices);

    LogTrace("l1tOmtfEventPrint")<<"simTraksHandle size "<<simTraksHandle.product()->size()<<std::endl;;
  }

  edm::Handle<TrackingParticleCollection> trackingParticleHandle;

  if(!trackingParticleToken.isUninitialized()) {
    event.getByToken(trackingParticleToken, trackingParticleHandle);
    LogTrace("l1tOmtfEventPrint")<<"trackingParticleHandle size "<<trackingParticleHandle.product()->size()<<std::endl;;
  }



  //todo do little better, move this assignment to constructor
  muonMatcher.setup(es);


  std::function<bool(const SimTrack& )> const& simTrackFilter = simTrackIsMuonInOmtfBx0;
  if(analysisType == "rate") {
    std::function<bool(const SimTrack& )> const& simTrackFilter = simTrackIsMuonInOmtf;
  }

  std::function<bool(const TrackingParticle& )> trackParticleFilter = trackingParticleIsMuonInOmtfEvent0;
  if(analysisType == "rate") {
    trackParticleFilter = trackingParticleIsMuonInOmtfBx0;
  }

  std::vector<MatchingResult> matchingResults;
  if(matchUsingPropagation) {
    if(!trackingParticleToken.isUninitialized())
      matchingResults = muonMatcher.match(ghostBustedCands, trackingParticleHandle.product(), trackParticleFilter);
    else if(!simTrackToken.isUninitialized())
      matchingResults = muonMatcher.match(ghostBustedCands, simTraksHandle.product(), simVertices.product(), simTrackFilter);
  }
  else {

    if(fillMatcherHists) {
      std::function<bool(const SimTrack& )> simTrackFilter = simTrackIsMuonInOmtfBx0;
      muonMatcher.fillHists(ghostBustedCands, simTraksHandle.product(), simTrackFilter);
    }
    else
      matchingResults = muonMatcher.matchWithoutPorpagation(ghostBustedCands, trackingParticleHandle.product(), trackParticleFilter);
  }

  candPerEvent->Fill(l1omtfHandle.product()->size(0));

  if(analysisType == "efficiency") {
    LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyze:"<<__LINE__<<std::endl;
    analyzeEfficiency(event, matchingResults);
  }
  else if(analysisType == "rate") {
    //analyzeRate(event, es);
    analyzeRate(event, matchingResults);
  }

}

void L1MuonAnalyzerOmtf::analyzeEfficiency(const edm::Event& event, std::vector<MatchingResult>& matchingResults) {
  LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyzeEfficiency"<<std::endl;

  for (auto& matchingResult: matchingResults ) {
    if(matchingResult.trackingParticle || matchingResult.simTrack) {
      ptGenHist->Fill(matchingResult.genPt);

      LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyze, sim track type "<<matchingResult.pdgId<<" simTrack pt "<<matchingResult.genPt<<std::endl;

      //const l1t::RegionalMuonCand* bestOmtfCand = nullptr;
      //unsigned int bestCandFiredLayersCnt = 0;

      L1MuonCand l1MuonCand;  //empty cand
      if(matchingResult.muonCand) {
        L1MuonCand l1MuonCand1(*(matchingResult.muonCand));
        l1MuonCand = l1MuonCand1;
        l1MuonCand.ptGev = hwPtToPtGeV(matchingResult.muonCand->hwPt() ); //TODO

        if(matchingResult.genPt >= 25) {
          if(l1MuonCand.hwQual > 0) {//removes candidates with abs(eta) > 1.24
            int firedLayers = matchingResult.muonCand->trackAddress().at(0);
            if( abs(l1MuonCand.ptGev) >= 18) {
              firedLayersEventCntOmtf->AddBinContent(firedLayers +1);
            }

            if(omtfNNEfficiencyAnalysers.size() && abs(hwPtToPtGeV(matchingResult.muonCand->trackAddress().at(10 + 2) ) ) >= 22) {//TODO nn pt for the p threshold 0.45, change if other is needed
              firedLayersEventCntNN->AddBinContent(firedLayers +1);
            }
          }
        }
      }
      else {
        LogTrace("l1tOmtfEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
      }

      for(auto& efficiencyAnalyser : omtfEfficiencyAnalysers) {
        efficiencyAnalyser->fill(matchingResult.genPt, matchingResult.genEta, matchingResult.genPhi, l1MuonCand);
      }

      for(unsigned int i = 0; i < omtfNNEfficiencyAnalysers.size(); i++) {
        if(matchingResult.muonCand) {
          l1MuonCand.ptGev = fabs(hwPtToPtGeV(matchingResult.muonCand->trackAddress().at(10 + i/9) ) ); //TODO check if abs is in the proper place, TODO i/2 because there are 2 Analysers, change if toehre addeds
        }
        omtfNNEfficiencyAnalysers[i]->fill(matchingResult.genPt, matchingResult.genEta, matchingResult.genPhi, l1MuonCand);
      }
    }
  }
}

double L1MuonAnalyzerOmtf::getDeltaR(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand) {

  return 0;
}


bool L1MuonAnalyzerOmtf::matched(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand) {
  if( (simTrackPtr->momentum().eta() * omtfCand.hwEta() ) > 0 ) { //eta has the same sign
    double globalPhi = l1t::MicroGMTConfiguration::calcGlobalPhi( omtfCand.hwPhi(), omtfCand.trackFinderType(), omtfCand.processor() );
    globalPhi = hwGmtPhiToGlobalPhi(globalPhi );

    if(globalPhi > M_PI)
      globalPhi = globalPhi -(2.*M_PI);

    double deltaPhi = abs(foldPhi(simTrackPtr->momentum().phi() - globalPhi) );
    if(simTrackPtr->momentum().pt() < 5) {
      if( deltaPhi < 1.2) {
        return true;
      }
    }
    if(simTrackPtr->momentum().pt() < 10) {
      if( deltaPhi < 0.7) {
        return true;
      }
    }
    if(simTrackPtr->momentum().pt() < 20) {
      if( deltaPhi < 0.3) {
        return true;
      }
    }
    if( deltaPhi < 0.2) {
      return true;
    }
  }
  return false;
}

void L1MuonAnalyzerOmtf::analyzeRate(const edm::Event& event, std::vector<MatchingResult>& matchingResults) { //, const edm::SimVertexContainer* simVertices
  LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyzeRate"<<std::endl;

  MatchingResult* bestOmtfCand = nullptr;
  for (auto& matchingResult: matchingResults ) {
    if(matchingResult.muonCand) {
      if(!bestOmtfCand || bestOmtfCand->muonCand->hwPt() < matchingResult.muonCand->hwPt()) {
        bestOmtfCand = &matchingResult;
      }
    }
    else {
      if(matchingResult.trackingParticle) {
        LogTrace("l1tOmtfEventPrint") <<" analyzeRate: no candidate for tracking particle pdgId"<<matchingResult.trackingParticle->pdgId()<<" pt "<<matchingResult.trackingParticle->pt() <<std::endl;
      }
    }
  }


  L1MuonCand l1MuonCand;  //empty cand
  if(bestOmtfCand) {
    LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyzeRate bestOmtfCand muonCand->hwPt "<<bestOmtfCand->muonCand->hwPt()<<" trackingParticle "<<bestOmtfCand->trackingParticle <<std::endl;

    L1MuonCand l1MuonCand1(*(bestOmtfCand->muonCand) );
    l1MuonCand = l1MuonCand1;
    l1MuonCand.ptGev = hwPtToPtGeV(bestOmtfCand->muonCand->hwPt() ); //TODO

    if(l1MuonCand.hwQual > 0) {//removes candidates with abs(eta) > 1.24
      int firedLayers = bestOmtfCand->muonCand->trackAddress().at(0);
      if(abs(l1MuonCand.ptGev) >= 18) {
        firedLayersEventCntOmtf->AddBinContent(firedLayers +1);
      }

      if(omtfNNRateAnalysers.size() &&  abs(hwPtToPtGeV(bestOmtfCand->muonCand->trackAddress().at(10 + 2) ) ) >= 22) {//TODO nn pt for the p threshold 0.45, change if other is needed
        firedLayersEventCntNN->AddBinContent(firedLayers +1);
      }
    }

    /*if(bestOmtfCand->muonCand->hwQual() >= 1 && bestOmtfCand->muonCand->hwPt() >= 41) { //41
      edm::LogVerbatim("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyzeRate bestOmtfCand hwPt "<<std::setw(3)<<bestOmtfCand->muonCand->hwPt()<<" pt GeV "<<l1MuonCand.ptGev
              <<" qual "<<std::setw(2)<<bestOmtfCand->muonCand->hwQual()
              <<" proc "<<std::setw(2)<<bestOmtfCand->muonCand->processor()<<" eta "<<std::setw(4)<<bestOmtfCand->muonCand->hwEta()<<" gloablEta "<<std::setw(8)<<bestOmtfCand->muonCand->hwEta() * 0.010875
              <<" hwPhi "<<std::setw(3)<<bestOmtfCand->muonCand->hwPhi()
              <<" globalPhi "<<std::setw(8)<<hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( bestOmtfCand->muonCand->hwPhi(), bestOmtfCand->muonCand->trackFinderType(), bestOmtfCand->muonCand->processor() ) )
              <<" fireadLayers "<<std::bitset<18>(bestOmtfCand->muonCand->trackAddress().at(0) )<<std::endl;

      if(bestOmtfCand->trackingParticle) {
        edm::LogVerbatim("l1tOmtfEventPrint") <<"trackingParticleIsMuonInOmtfBx0, pdgId "<<std::setw(3)<<bestOmtfCand->trackingParticle->pdgId()<<" pt "<<std::setw(9)<<bestOmtfCand->trackingParticle->pt()
                  <<" eta "<<std::setw(9)<<bestOmtfCand->trackingParticle->momentum().eta()<<" phi "<<std::setw(9)<<bestOmtfCand->trackingParticle->momentum().phi()<<" event "<<bestOmtfCand->trackingParticle->eventId().event()
                  <<"\nmatching: deltaEta "<<std::setw(8)<<bestOmtfCand->deltaEta<<" deltaPhi "<<std::setw(8)<<bestOmtfCand->deltaPhi<<" matchingLL "<<std::setw(8)<<bestOmtfCand->matchingLikelihood<<std::endl;


        if(bestOmtfCand->trackingParticle->parentVertex().isNonnull() ) {
          edm::LogVerbatim("l1tOmtfEventPrint")<<"parentVertex Rho "<<bestOmtfCand->trackingParticle->parentVertex()->position().Rho()
                  <<" eta "<<bestOmtfCand->trackingParticle->parentVertex()->position().eta()<<" phi "<<bestOmtfCand->trackingParticle->parentVertex()->position().phi()<<std::endl;

          for(auto& parentTrack : bestOmtfCand->trackingParticle->parentVertex()->sourceTracks() ) {
            edm::LogVerbatim("l1tOmtfEventPrint")<<"parentTrackPdgId "<<parentTrack->pdgId()<<std::endl;
          }
        }
      }
      else {
        edm::LogVerbatim("l1tOmtfEventPrint")<<"no matched trackingParticle"<<std::endl;
      }
      edm::LogVerbatim("l1tOmtfEventPrint")<<std::endl;
    }*/


    for(auto& omtfCandsMatchingAnalyser : omtfCandsMatchingAnalysers) {
      omtfCandsMatchingAnalyser->fill(l1MuonCand, bestOmtfCand->trackingParticle);
    }


    for(auto& omtfRateAnalyser : omtfRateAnalysers) {
      omtfRateAnalyser->fill(l1MuonCand);
    }

    /*
     * here is little bit not good, because the  bestOmtfCand is the one with the highest omtf pt, and not the nn pt,
     * so in some rare case the nn pt of the chose cand might be smaller than the nn pt of some other cand
     */
    for(unsigned int i = 0; i < omtfNNRateAnalysers.size(); i++) {
      if(bestOmtfCand) {
        l1MuonCand.ptGev = fabs(hwPtToPtGeV(bestOmtfCand->muonCand->trackAddress().at(10 + i/3) ) ); //TODO check if abs is in the proper place, TODO watch aout the i
      }
      omtfNNRateAnalysers[i]->fill( l1MuonCand);
      omtfNNCandsMatchingAnalysers.at(i)->fill(l1MuonCand, bestOmtfCand->trackingParticle);
    }
  }
  else {
    //LogTrace("l1tOmtfEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
    //if there is no candidates, nothing is filled
  }

}

void L1MuonAnalyzerOmtf::analyzeRate(const edm::Event& event, const edm::EventSetup& es) {
  LogTrace("l1tOmtfEventPrint") <<"L1MuonAnalyzerOmtf::analyzeRate"<<std::endl;

  /*edm::Handle<edm::SimTrackContainer> simTraksHandle;
  event.getByToken(simTrackToken, simTraksHandle); */

  edm::Handle<l1t::RegionalMuonCandBxCollection> l1omtfHandle;
  event.getByToken(omtfToken, l1omtfHandle);

  candPerEvent->Fill(l1omtfHandle.product()->size(0));

  const l1t::RegionalMuonCand* bestOmtfCand = nullptr;
  unsigned int bestCandFiredLayersCnt = 0;

  for(l1t::RegionalMuonCandBxCollection::const_iterator omtfCand = l1omtfHandle.product()->begin(0);
      omtfCand != l1omtfHandle.product()->end(0); ++omtfCand) {
    int refLayer = (int)omtfCand->trackAddress().at(1);
    int layerHits = (int)omtfCand->trackAddress().at(0);
    std::bitset<18> layerHitBits(layerHits);

    double globalPhi = l1t::MicroGMTConfiguration::calcGlobalPhi( omtfCand->hwPhi(), omtfCand->trackFinderType(), omtfCand->processor() )* 2. * M_PI / 576;
    if(globalPhi > M_PI)
      globalPhi = globalPhi -(2.*M_PI);
    LogTrace("l1tOmtfEventPrint") << "L1MuonAnalyzerOmtf::"<<__FUNCTION__<<":"<<__LINE__
        <<" omtf pt "<<omtfCand->hwPt()<<" omtf qual "<<omtfCand->hwQual()<<" omtf hwEta "<<omtfCand->hwEta()<<" omtf hwPhi "<<omtfCand->hwPhi()
        <<" eta "<< (omtfCand->hwEta()*0.010875)
        <<" phi "<<globalPhi
        <<" refLayer "<<refLayer<<" "<<layerHitBits<<endl;

    if(!bestOmtfCand) {
      bestOmtfCand = &(*omtfCand);
      bestCandFiredLayersCnt = layerHitBits.count();
    }
    else {
      if(omtfCand->hwQual() > bestOmtfCand->hwQual() ) {
        bestOmtfCand = &(*omtfCand);
        bestCandFiredLayersCnt = layerHitBits.count();
      }
      else if(omtfCand->hwQual() == bestOmtfCand->hwQual() ) {
        if(layerHitBits.count() > bestCandFiredLayersCnt ) {
          bestOmtfCand = &(*omtfCand);
          bestCandFiredLayersCnt = layerHitBits.count();
        }
        else if(layerHitBits.count() == bestCandFiredLayersCnt ) {
          if(omtfCand->hwPt() > bestOmtfCand->hwPt() ) {
            bestOmtfCand = &(*omtfCand);
            bestCandFiredLayersCnt = layerHitBits.count();
          }
        }
      }
    }
  }

  L1MuonCand l1MuonCand;  //empty cand
  if(bestOmtfCand) {
    L1MuonCand l1MuonCand1(*bestOmtfCand);
    l1MuonCand = l1MuonCand1;
    l1MuonCand.ptGev = hwPtToPtGeV(bestOmtfCand->hwPt() ); //TODO

    if(l1MuonCand.hwQual > 0) {//removes candidates with abs(eta) > 1.24
      int firedLayers = bestOmtfCand->trackAddress().at(0);
      if(abs(l1MuonCand.ptGev) >= 18) {
        firedLayersEventCntOmtf->AddBinContent(firedLayers +1);
      }

      if(omtfNNRateAnalysers.size() && abs(hwPtToPtGeV(bestOmtfCand->trackAddress().at(10 + 2) ) ) >= 22) {//TODO nn pt for the p threshold 0.45, change if other is needed
        firedLayersEventCntNN->AddBinContent(firedLayers +1);
      }
    }

  }
  else {
    LogTrace("l1tOmtfEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
  }

  for(auto& omtfRateAnalyser : omtfRateAnalysers) {
    omtfRateAnalyser->fill(l1MuonCand);
  }

  for(unsigned int i = 0; i < omtfNNRateAnalysers.size(); i++) {
    if(bestOmtfCand) {
      l1MuonCand.ptGev = fabs(hwPtToPtGeV(bestOmtfCand->trackAddress().at(10 + i/3) ) ); //TODO check if abs is in the proper place, TODO watch aout the i
    }
    omtfNNRateAnalysers[i]->fill( l1MuonCand);
  }


}




void L1MuonAnalyzerOmtf::endJob() {
  muonMatcher.saveHists();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MuonAnalyzerOmtf);

}
