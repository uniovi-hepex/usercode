/*
 * L1MuonAnalyzerOmtf.cc
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/plugins/L1MuonAnalyzerOmtf.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

namespace L1MuAn {


L1MuonAnalyzerOmtf::L1MuonAnalyzerOmtf(const edm::ParameterSet& edmCfg) {
  omtfToken = consumes<l1t::RegionalMuonCandBxCollection >(edmCfg.getParameter<edm::InputTag>("L1OMTFInputTag"));

  simTrackToken =  consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits")); //TODO which is correct?

  edm::Service<TFileService> fs;

  candPerEvent = fs->make<TH1D>("candPerEvent", "candPerEvent", 21, -0.5, 20.5);

  analysisType = edmCfg.getParameter< string >("analysisType");

  if(edmCfg.exists("nn_pThresholds") )
    nn_pThresholds = edmCfg.getParameter<vector<double> >("nn_pThresholds");

  if(analysisType == "efficiency") {
    TFileDirectory subDir = fs->mkdir("efficiency");
    omtfEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, "omtfHighQ", 0.82, 1.24, 12, 200, 0, 200));
    omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, "omtfHighQ", 0.82, 1.24, 12, 10., 1, 100));
    omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, "omtfHighQ", 12, 10., 1, 100));

    omtfEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, "omtfLowQ",  0.82, 1.24, 4, 200, 0, 200));
    omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, "omtfLowQ", 0.82, 1.24, 4, 10., 1, 100));
    omtfEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, "omtfLowQ", 4, 10., 1, 100));

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"omtf_nn_q12_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 12, 200, 0, 200));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 12, 10., 1, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 12, 10., 1, 100));

      ostr.str("");
      ostr<<"omtf_nn_q4_pTresh_"<<nn_pThreshold;
      omtfNNEfficiencyAnalysers.emplace_back(new PtGenVsPtCand(subDir, ostr.str().c_str(), 0.82, 1.24, 4, 200, 0, 200));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsPhi(subDir, ostr.str().c_str(), 0.82, 1.24, 4, 10., 1, 100));
      omtfNNEfficiencyAnalysers.emplace_back(new EfficiencyVsEta(subDir, ostr.str().c_str(), 4, 10., 1, 100));
    }
  }

  else if(analysisType == "rate") {
    TFileDirectory subDirRate = fs->mkdir("rate");

    TFileDirectory subDir = subDirRate.mkdir("omtfHighQ");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 12, 200, 0, 100));

    subDir = subDirRate.mkdir("omtfLowQ");
    omtfRateAnalysers.emplace_back(new  RateAnalyser(subDir, "", 1, 200, 0, 100));

    for(auto& nn_pThreshold : nn_pThresholds) {
      std::ostringstream ostr;
      ostr<<"omtf_q1_nn_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 1, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"omtf_q4_nn_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 4, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;

      ostr.str("");
      ostr<<"omtf_q12_nn_pTresh_"<<nn_pThreshold;
      subDir = subDirRate.mkdir(ostr.str().c_str());
      omtfNNRateAnalysers.emplace_back(new RateAnalyser(subDir, "", 12, 200, 0, 100));
      edm::LogImportant("l1tMuBayesEventPrint") <<" adding omtfNNEfficiencyAnalysers, nn_pThreshold "<<nn_pThreshold<<std::endl;
    }
  }
}


L1MuonAnalyzerOmtf::~L1MuonAnalyzerOmtf() {
  // TODO Auto-generated destructor stub
}

void L1MuonAnalyzerOmtf::beginJob() {

}

void L1MuonAnalyzerOmtf::analyze(const edm::Event& event, const edm::EventSetup& es) {
  if(analysisType == "efficiency") {
    analyzeEfficiency(event, es);
  }
  else if(analysisType == "rate") {
    analyzeRate(event, es);
  }
}

void L1MuonAnalyzerOmtf::analyzeEfficiency(const edm::Event& event, const edm::EventSetup& es) {
  LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::analyzeEfficiency"<<std::endl;

  edm::Handle<edm::SimTrackContainer> simTraksHandle;
  event.getByToken(simTrackToken, simTraksHandle);

  edm::Handle<l1t::RegionalMuonCandBxCollection> l1omtfHandle;
  event.getByToken(omtfToken, l1omtfHandle);

  for (unsigned int iSimTrack = 0; iSimTrack != simTraksHandle->size(); iSimTrack++ ) {
    edm::Ptr< SimTrack > simTrackPtr(simTraksHandle, iSimTrack);


    if (abs(simTrackPtr->type()) == 13 || abs(simTrackPtr->type()) == 1000015 ) {  //|| tpPtr->pt() > 20 //todo 1000015 is stau
      //only muons
    }
    else
      continue;

    /*    if (tpPtr->pt() < TP_minPt) TODO
      continue;
    if (fabs(tpPtr->eta()) > TP_maxEta)
      continue;*/

    if( (fabs(simTrackPtr->momentum().eta()) >= etaCutFrom ) && (fabs(simTrackPtr->momentum().eta()) <= etaCutTo) ) { //TODO it duplicates the above condition,
    }
    else
      continue;

    /*
    if(tpPtr->vertex().rho() > TP_maxRho)
      continue;
     */

    if(simTrackPtr->eventId().bunchCrossing() != 0)
      continue;


    LogTrace("l1tMuBayesEventPrint") <<"L1MuonAnalyzerOmtf::analyze, sim track type "<<simTrackPtr->type()<<" simTrack pt "<<simTrackPtr->momentum().pt()<<std::endl;

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

      if( matched(simTrackPtr, *omtfCand)) {
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
    }

    L1MuonCand l1MuonCand;  //empty cand
    if(bestOmtfCand) {
      L1MuonCand l1MuonCand1(*bestOmtfCand, bestCandFiredLayersCnt);
      l1MuonCand = l1MuonCand1;
      l1MuonCand.ptGev = hwPtToPtGeV(bestOmtfCand->hwPt() ); //TODO
    }
    else {
      LogTrace("l1tMuBayesEventPrint") <<" no matching candidate!!!!!!!!!!!!!!!!!" <<std::endl;
    }

    for(auto& efficiencyAnalyser : omtfEfficiencyAnalysers) {
      efficiencyAnalyser->fill(simTrackPtr->momentum().pt(), simTrackPtr->momentum().eta(), simTrackPtr->momentum().phi(), l1MuonCand);
    }

    for(unsigned int i = 0; i < omtfNNEfficiencyAnalysers.size(); i++) {
      if(bestOmtfCand) {
        l1MuonCand.ptGev = fabs(hwPtToPtGeV(bestOmtfCand->trackAddress().at(10 + i/6) ) ); //TODO check if abs is in the proper place, TODO i/2 because there are 2 Analysers, change if toehre addeds
      }
      omtfNNEfficiencyAnalysers[i]->fill(simTrackPtr->momentum().pt(), simTrackPtr->momentum().eta(), simTrackPtr->momentum().phi(), l1MuonCand);
    }

  }
}

double L1MuonAnalyzerOmtf::getDeltaR(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand) {

  return 0;
}

double hwPhiToGlobalPhi(int phi) {
  double phiGmtUnit = 2. * M_PI / 576.;
  return phi * phiGmtUnit;
}

double foldPhi(double phi) {
  if(phi > M_PI)
    return (phi - 2 * M_PI );
  else if(phi < -M_PI)
    return (phi +  2 * M_PI );

  return phi;
}

bool L1MuonAnalyzerOmtf::matched(const edm::Ptr< SimTrack >& simTrackPtr, const l1t::RegionalMuonCand& omtfCand) {
  if( (simTrackPtr->momentum().eta() * omtfCand.hwEta() ) > 0 ) { //eta has the same sign
    double globalPhi = l1t::MicroGMTConfiguration::calcGlobalPhi( omtfCand.hwPhi(), omtfCand.trackFinderType(), omtfCand.processor() );
    globalPhi = hwPhiToGlobalPhi(globalPhi );

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
    L1MuonCand l1MuonCand1(*bestOmtfCand, bestCandFiredLayersCnt);
    l1MuonCand = l1MuonCand1;
    l1MuonCand.ptGev = hwPtToPtGeV(bestOmtfCand->hwPt() ); //TODO
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

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MuonAnalyzerOmtf);

}
