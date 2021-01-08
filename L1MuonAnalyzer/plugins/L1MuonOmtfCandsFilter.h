/*
 * L1MuonOmtfCandsFilter.h
 *
 *  Created on: Jan 8, 2021
 *      Author: kbunkow
 */

#ifndef PLUGINS_L1MUONOMTFCANDSFILTER_H_
#define PLUGINS_L1MUONOMTFCANDSFILTER_H_

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

namespace L1MuAn {
class L1MuonOmtfCandsFilter: public edm::EDFilter {
public:
  L1MuonOmtfCandsFilter(const edm::ParameterSet& edmCfg);
  virtual ~L1MuonOmtfCandsFilter();

  //void beginJob() override;

  void endJob() override;

  //void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;

  bool filter(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection > omtfToken;

  int qualityCut = 12;
  int hwPtCut = 41;

  int acceptedEvents = 0;
};

} //namespace L1MuAn

#endif /* PLUGINS_L1MUONOMTFCANDSFILTER_H_ */
