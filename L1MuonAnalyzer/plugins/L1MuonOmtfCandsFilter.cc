/*
 * L1MuonOmtfCandsFilter.cpp
 *
 *  Created on: Jan 8, 2021
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/plugins/L1MuonOmtfCandsFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <bitset>

namespace L1MuAn {
L1MuonOmtfCandsFilter::L1MuonOmtfCandsFilter(const edm::ParameterSet& edmCfg) {
  omtfToken = consumes<l1t::RegionalMuonCandBxCollection >(edmCfg.getParameter<edm::InputTag>("L1OMTFInputTag"));

  qualityCut = edmCfg.getParameter<int>("qualityCut");
  hwPtCut = edmCfg.getParameter<int>("hwPtCut");

  edm::LogImportant("l1tOmtfEventPrint") <<" L1MuonOmtfCandsFilter: line "<<__LINE__<<" qualityCut "<<qualityCut<<" hwPtCut "<<hwPtCut<<std::endl;

  edm::Service<TFileService> fs;

  candPerEvent = fs->make<TH1D>("candPerEvent", "candPerEvent", 21, -0.5, 20.5);
}

L1MuonOmtfCandsFilter::~L1MuonOmtfCandsFilter() {
  // TODO Auto-generated destructor stub
}

bool L1MuonOmtfCandsFilter::filter(edm::Event& event, const edm::EventSetup&) {
edm::Handle<l1t::RegionalMuonCandBxCollection> l1omtfHandle;
 event.getByToken(omtfToken, l1omtfHandle);

 LogTrace("l1tOmtfEventPrint")<<std::endl<<"\nL1MuonAnalyzerOmtf::analyze. event "<<event.id().event() <<" number of candidates "<<l1omtfHandle.product()->size(0)<<" ---------------------------"<<std::endl;

 const l1t::RegionalMuonCandBxCollection* mtfCands = l1omtfHandle.product();

 candPerEvent->Fill(l1omtfHandle.product()->size(0));

 bool accept = false;
 for(unsigned int  i1 = 0; i1 < mtfCands->size(0); ++i1) {
   if(mtfCands->at( 0, i1 ).hwQual() >= qualityCut && mtfCands->at( 0, i1 ).hwPt() >= hwPtCut) {
     accept = true;
   }

   LogTrace("l1tOmtfEventPrint") <<"L1MuonOmtfCandsFilter::filter accept "<<accept<<" candidate pt "<<std::setw(3)<<mtfCands->at( 0, i1 ).hwPt()<<" qual "<<std::setw(2)<<mtfCands->at( 0, i1 ).hwQual()
       <<" proc "<<std::setw(2)<<mtfCands->at( 0, i1 ).processor()<<" eta "<<std::setw(4)<<mtfCands->at( 0, i1 ).hwEta()<<" gloablEta "<<std::setw(8)<<mtfCands->at( 0, i1 ).hwEta() * 0.010875
       <<" hwPhi "<<std::setw(3)<<mtfCands->at( 0, i1 ).hwPhi()
       //<<" globalPhi "<<std::setw(8)<<hwGmtPhiToGlobalPhi(l1t::MicroGMTConfiguration::calcGlobalPhi( mtfCands->at( 0, i1 ).hwPhi(), mtfCands->at( 0, i1 ).trackFinderType(), mtfCands->at( 0, i1 ).processor() ) )
       <<" fireadLayers "<<std::bitset<18>(mtfCands->at( 0, i1 ).trackAddress().at(0) )
       <<std::endl;
 }
 if(accept)
   acceptedEvents++;
 return accept;
}

void L1MuonOmtfCandsFilter::endJob() {
  edm::LogImportant("l1tOmtfEventPrint") <<" L1MuonOmtfCandsFilter: line "<<__LINE__<<" acceptedEvents "<<acceptedEvents<<std::endl;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MuonOmtfCandsFilter);

} //namespace L1MuAn
