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
  omtfToken = consumes<l1t::RegionalMuonCandBxCollection >(edmCfg.getParameter<edm::InputTag>("L1OMTFInputTag"));

  simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //TODO which is correct?

  simVertexesToken =  consumes<edm::SimVertexContainer>(edmCfg.getParameter<edm::InputTag>("simVertexesTag"));

  edm::Service<TFileService> fs;

  candPerEvent = fs->make<TH1D>("candPerEvent", "candPerEvent", 21, -0.5, 20.5);

  analysisType = edmCfg.getParameter< string >("analysisType");

  int binsCnt = (1<<18)-1;
  firedLayersEventCntOmtf = fs->make<TH1S>("firedLayersEventCntOmtf", "firedLayersEventCntOmtf", binsCnt, 0, binsCnt);
  firedLayersEventCntNN = fs->make<TH1S>("firedLayersEventCntNN", "firedLayersEventCntNN", binsCnt, 0, binsCnt);

  if(edmCfg.exists("nn_pThresholds") )
    nn_pThresholds = edmCfg.getParameter<vector<double> >("nn_pThresholds");

  if(analysisType == "efficiency") {
    double ptGenCut = 25;
    double ptL1Cut = 20;
    TFileDirectory subDir = fs->mkdir("efficiency");
    auto addOmtfAnalysers = [&](TFileDirectory& subDir, std::string name, int qualityCut) {
      omtfEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, name, 0.82, 1.24, qualityCut, 200, 0, 200));
      omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, name, 0.82, 1.24, qualityCut, ptGenCut, ptL1Cut, 100));
      omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, name, qualityCut, ptGenCut, ptL1Cut, 210));
      omtfEfficiencyAnalysers.emplace_back(new LikelihoodDistribution(subDir, name, qualityCut, ptGenCut, ptL1Cut, 120));
    };

    addOmtfAnalysers(subDir, "omtf_q12", 12);
    addOmtfAnalysers(subDir, "omtf_q8", 8);
    addOmtfAnalysers(subDir, "omtf_q4", 4);
    addOmtfAnalysers(subDir, "omtf_q1", 1);

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"nn_omtf_q1_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 1, 200, 0, 200));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 1, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 1, ptGenCut, ptL1Cut, 210));

      ostr.str("");
      ostr<<"nn_omtf_q4_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 4, 200, 0, 200));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 4, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 4, ptGenCut, ptL1Cut, 210));

      ostr.str("");
      ostr<<"nn_omtf_q12_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 12, 200, 0, 200));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 12, ptGenCut, ptL1Cut, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 12, ptGenCut, ptL1Cut, 210));
    }
  }

  else if(analysisType == "rate") {
    TFileDirectory subDirRate = fs->mkdir("rate");

    TFileDirectory subDir = subDirRate.mkdir("omtf_q12");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 12, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q8");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 8, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q4");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 4, 200, 0, 100));

    subDir = subDirRate.mkdir("omtf_q1");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 1, 200, 0, 100));

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"nn_omtf_q1_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 1, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"nn_omtf_q4_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 4, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"nn_omtf_q12_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 12, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNRateAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
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

    LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::ghostBust pt "<<mtfCands->at( 0, i1 ).hwPt()<<" qual "<<mtfCands->at( 0, i1 ).hwQual()
        <<" proc "<<mtfCands->at( 0, i1 ).processor()<<" eta "<<mtfCands->at( 0, i1 ).hwEta()<<" gloablEta "<<mtfCands->at( 0, i1 ).hwEta() * 0.010875
        <<" hwPhi "<<mtfCands->at( 0, i1 ).hwPhi()
        <<" globalPhi "<<hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCands->at( 0, i1 ).hwPhi(), mtfCands->at( 0, i1 ).trackFinderType(), mtfCands->at( 0, i1 ).processor() ) )
        <<" fireadLayers "<<std::bitset<18>(mtfCands->at( 0, i1 ).trackAddress().at(0) )
        <<" isKilled "<<isKilled.test(i1)<<std::endl;
  }

  if(resultCands.size() >= 3)
    LogTrace("l1tMuBayesEventPrint")<<" ghost !!!!!! " <<std::endl;
  LogTrace("l1tMuBayesEventPrint") <<std::endl;

  return resultCands;
}


bool simTrackIsMuonInOmtf(const SimTrack& simTrack) {
  if (abs(simTrack.type()) == 13 || abs(simTrack.type()) == 1000015 ) {  //|| tpPtr->pt() > 20 //todo 1000015 is stau
    //only muons
  }
  else
    return false;

  if(simTrack.momentum().pt() < 3) //in the overlap, the propagation of muons with pt less then ~3.2 fails - the actual threshold depends slightly on eta,
    return false;

  if( (fabs(simTrack.momentum().eta()) >= 0.82 ) && (fabs(simTrack.momentum().eta()) <= 1.24) ) {
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

void L1MuonAnalyzerOmtf::analyze(const edm::Event& event, const edm::EventSetup& es) {

  edm::Handle<l1t::RegionalMuonCandBxCollection> l1omtfHandle;
  event.getByToken(omtfToken, l1omtfHandle);

  LogTrace("l1tMuBayesEventPrint")<<std::endl<<"\nL1MuonAnalyzerOmtf::analyze. event "<<event.id().event() <<" number of candidates "<<l1omtfHandle.product()->size(0)<<" ---------------------------"<<std::endl;

  std::vector<const l1t::RegionalMuonCand*> ghostBustedCands = ghostBust(l1omtfHandle.product());


  edm::Handle<edm::SimTrackContainer> simTraksHandle;
  event.getByToken(simTrackToken, simTraksHandle);

  edm::Handle<edm::SimVertexContainer> simVertices;
  event.getByToken(simVertexesToken, simVertices);


  //todo do little better, move this assignment to constructor
  muonMatcher.setup(es);
  std::function<bool(const SimTrack& )> const& simTrackFilter = simTrackIsMuonInOmtfBx0;

  if(analysisType == "rate") {
    std::function<bool(const SimTrack& )> const& simTrackFilter = simTrackIsMuonInOmtf;
  }


  std::vector<MatchingResult> matchingResults = muonMatcher.match(ghostBustedCands, simTraksHandle.product(), simVertices.product(), simTrackFilter);

  candPerEvent->Fill(l1omtfHandle.product()->size(0));

  if(analysisType == "efficiency") {
    analyzeEfficiency(event, matchingResults);
  }
  else if(analysisType == "rate") {
    analyzeRate(event, es);
  }

}

void L1MuonAnalyzerOmtf::analyzeEfficiency(const edm::Event& event, std::vector<MatchingResult>& matchingResults) {
  LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::analyzeEfficiency"<<std::endl;

  for (auto& matchingResult: matchingResults ) {
    if(matchingResult.simTrack) {

      LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::analyze, sim track type "<<matchingResult.simTrack->type()<<" simTrack pt "<<matchingResult.simTrack->momentum().pt()<<std::endl;

      //const l1t::RegionalMuonCand* bestOmtfCand = nullptr;
      //unsigned int bestCandFiredLayersCnt = 0;

      L1MuonCand l1MuonCand;  //empty cand
      if(matchingResult.muonCand) {
        L1MuonCand l1MuonCand1(*(matchingResult.muonCand));
        l1MuonCand = l1MuonCand1;
        l1MuonCand.ptGev = hwPtToPtGeV(matchingResult.muonCand->hwPt() ); //TODO

        if(matchingResult.simTrack->momentum().pt() >= 22) {
          if(l1MuonCand.hwQual > 0) {//removes candidates with abs(eta) > 1.24
            int firedLayers = matchingResult.muonCand->trackAddress().at(0);
            if( abs(l1MuonCand.ptGev) >= 20) {
              firedLayersEventCntOmtf->AddBinContent(firedLayers +1);
            }

            if( abs(hwPtToPtGeV(matchingResult.muonCand->trackAddress().at(10 + 2) ) ) >= 20) {//TODO nn pt for the p threshold 0.45, change if other is needed
              firedLayersEventCntNN->AddBinContent(firedLayers +1);
            }
          }
        }
      }
      else {
        LogTrace("l1tMuBayesEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
      }

      for(auto& efficiencyAnalyser : omtfEfficiencyAnalysers) {
        efficiencyAnalyser->fill(matchingResult.simTrack->momentum().pt(), matchingResult.simTrack->momentum().eta(), matchingResult.simTrack->momentum().phi(), l1MuonCand);
      }

      for(unsigned int i = 0; i < omtfNNEfficiencyAnalysers.size(); i++) {
        if(matchingResult.muonCand) {
          l1MuonCand.ptGev = fabs(hwPtToPtGeV(matchingResult.muonCand->trackAddress().at(10 + i/9) ) ); //TODO check if abs is in the proper place, TODO i/2 because there are 2 Analysers, change if toehre addeds
        }
        omtfNNEfficiencyAnalysers[i]->fill(matchingResult.simTrack->momentum().pt(), matchingResult.simTrack->momentum().eta(), matchingResult.simTrack->momentum().phi(), l1MuonCand);
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



void L1MuonAnalyzerOmtf::analyzeRate(const edm::Event& event, const edm::EventSetup& es) {
  LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::analyzeRate"<<std::endl;

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
    LogTrace("l1tMuBayesEventPrint") << "MuCorrelatorAnalyzer::"<<__FUNCTION__<<":"<<__LINE__
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
      if(abs(l1MuonCand.ptGev) >= 20) {
        firedLayersEventCntOmtf->AddBinContent(firedLayers +1);
      }

      if( abs(hwPtToPtGeV(bestOmtfCand->trackAddress().at(10 + 2) ) ) >= 20) {//TODO nn pt for the p threshold 0.45, change if other is needed
        firedLayersEventCntNN->AddBinContent(firedLayers +1);
      }
    }

  }
  else {
    LogTrace("l1tMuBayesEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
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
