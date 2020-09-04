/*
 * RateAnalyser.cc
 *
 *  Created on: Apr 1, 2020
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/interface/RateAnalyser.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace L1MuAn {

RateAnalyser::RateAnalyser(TFileDirectory subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo): qualityCut(qualityCut) {
  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<"candPt_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand pt, "<<" quality >= "<<qualityCut<<";cand pt [GeV]; #events";

  candPt = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut1GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 1 GeV "<<" quality >= "<<qualityCut<<";eta; #events";

  candEtaPtCut1 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut10GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 10 GeV "<<" quality >= "<<qualityCut<<";eta; #events";

  candEtaPtCut10 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut20GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 20 GeV "<<" quality >= "<<qualityCut<<";eta; #events";

  candEtaPtCut20 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

}

RateAnalyser::~RateAnalyser() {
  // TODO Auto-generated destructor stub
}

///center of eta bin
double hwEtaToEta(int hwEta) {
  double etaUnit = 0.010875; //TODO from the interface note, should be defined somewhere globally

  return (hwEta * etaUnit);
}

void RateAnalyser::fill(L1MuonCand& l1MuonCand) {
  auto candPtGev = l1MuonCand.ptGev;
  if(candPtGev >= candPt->GetXaxis()->GetXmax())
    candPtGev = candPt->GetXaxis()->GetXmax() - 0.01;

  if(l1MuonCand.hwQual >= qualityCut) {
    candPt->Fill(candPtGev);

    auto eta = hwEtaToEta(l1MuonCand.hwEta);
    if(candPtGev >= 1.)
      candEtaPtCut1->Fill(eta);

    if(candPtGev >= 10.)
      candEtaPtCut10->Fill(eta);

    if(candPtGev >= 20.)
      candEtaPtCut20->Fill(eta);
  }

}

void RateAnalyser::write() {
  candPt->Write();

  candEtaPtCut1->Write();
  candEtaPtCut10->Write();
  candEtaPtCut20->Write();
}

/*
PtGenVsPtCand promptMuonsPtGenVsPtCand;
PtGenVsPtCand nonPromptMuonsPtGenVsPtCand;

PtGenVsPtCand muonsFromPionsPtGenVsPtCand;
PtGenVsPtCand muonsFromKaonsPtGenVsPtCand;*/

CandsMatchingAnalyser::CandsMatchingHists::CandsMatchingHists(TFileDirectory& parrentSubDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo):
    subDir(parrentSubDir.mkdir(name)),
    qualityCut(qualityCut),
    rateAn(subDir, name, qualityCut, nBins, binsFrom, binsTo),
    ptGenVsPtCand(subDir, "", 0.82, 1.24, qualityCut, 100, 0, 100)
{
  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<"_simVertexRhoVsPtGen_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" simVertexRho Vs genPt, "<<" quality >= "<<qualityCut<<"; genPt [GeV]; Rho [cm]";

  simVertexRhoVsPtGen = subDir.make<TH2I>(histName.str().c_str(), histTitle.str().c_str(), 100, binsFrom, binsTo, 30, 0, 300);
}

void CandsMatchingAnalyser::CandsMatchingHists::fill(L1MuonCand& l1MuonCand, const TrackingParticle* matchedTrackingParticle,  double simVertexRho) {
  rateAn.fill(l1MuonCand);

  ptGenVsPtCand.fill(matchedTrackingParticle->pt(), matchedTrackingParticle->eta(), matchedTrackingParticle->phi(), l1MuonCand);

  if(l1MuonCand.hwQual >= qualityCut) {
    simVertexRhoVsPtGen->Fill(matchedTrackingParticle->pt(), simVertexRho);
  }

}

void CandsMatchingAnalyser::CandsMatchingHists::write() {
  rateAn.write();
  ptGenVsPtCand.write();
  simVertexRhoVsPtGen->Write();
}

CandsMatchingAnalyser::CandsMatchingAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo):
    promptMuons(subDir, "promptMuons", qualityCut, nBins, binsFrom, binsTo),
    nonPromptMuons(subDir, "nonPromptMuons", qualityCut, nBins, binsFrom, binsTo),
    muonsFromPions(subDir, "muonsFromPions", qualityCut, nBins, binsFrom, binsTo),
    muonsFromKaons(subDir, "muonsFromKaons", qualityCut, nBins, binsFrom, binsTo),
    notMatchedRateAn(subDir.mkdir("notMatched"), "notMatched", qualityCut, nBins, binsFrom, binsTo),
    likelihoodDistribution(subDir, name, qualityCut, 0, 20, 120)
{

  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<name<<"simVertexRhoVsPtGen_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" simVertexRho Vs genPt, "<<" quality >= "<<qualityCut<<"; genPt [GeV]; Rho [cm]";

  simVertexRhoVsPtGen = subDir.make<TH2I>(histName.str().c_str(), histTitle.str().c_str(), 100, binsFrom, binsTo, 30, 0, 300);

}

void CandsMatchingAnalyser::fill(L1MuonCand& l1MuonCand, const TrackingParticle* matchedTrackingParticle) {
  if(matchedTrackingParticle) {
    int parentTrackPdgId = 0;


    if(matchedTrackingParticle->parentVertex().isNonnull() ) {
      LogTrace("l1tMuBayesEventPrint")<<" CandsMatchingAnalyser parentVertex Rho "<<matchedTrackingParticle->parentVertex()->position().Rho()<<std::endl;
      for(auto& parentTrack : matchedTrackingParticle->parentVertex()->sourceTracks() ) {
        parentTrackPdgId = parentTrack->pdgId();
        LogTrace("l1tMuBayesEventPrint")<<" CandsMatchingAnalyser parentTrackPdgId "<<parentTrackPdgId<<std::endl;
      }
    }

    if (parentTrackPdgId == 0) {
      promptMuons.fill(l1MuonCand, matchedTrackingParticle, 0);
    }
    else {
      double simVertexRho = matchedTrackingParticle->parentVertex()->position().Rho();

      simVertexRhoVsPtGen->Fill(matchedTrackingParticle->pt(), simVertexRho);

      if(simVertexRho > 10) {
        nonPromptMuons.fill(l1MuonCand, matchedTrackingParticle, simVertexRho);
      }
      else {
        promptMuons.fill(l1MuonCand, matchedTrackingParticle, simVertexRho);
      }

      if(abs(parentTrackPdgId) == 211) {
        muonsFromPions.fill(l1MuonCand, matchedTrackingParticle, simVertexRho);
      }
      else if(abs(parentTrackPdgId) == 321) {
        muonsFromKaons.fill(l1MuonCand, matchedTrackingParticle, simVertexRho);
      }
    }
  }
  else {
    notMatchedRateAn.fill(l1MuonCand);
  }

  likelihoodDistribution.fill(0, 0, 0, l1MuonCand);
}

/*void CandsMatchingAnalyser::fill(L1MuonCand& l1MuonCand, const SimTrack* matchedSimTrack, const edm::SimVertexContainer* simVertices) {
  if(matchedSimTrack) {
    int vtxInd = matchedSimTrack->vertIndex();
    if (vtxInd < 0) {
      promptMuonsRateAn.fill(l1MuonCand);
      simVertexRhoVsPtGen->Fill(l1MuonCand.ptGev, 0.);
    }
    else {
      const SimVertex& simVertex = simVertices->at(vtxInd);

      //simVertex.parentIndex();
      //matchedSimTrack->trackId();

      double simVertexRho = simVertex.position().Rho();

      simVertexRhoVsPtGen->Fill(l1MuonCand.ptGev, simVertexRho);

      if(simVertexRho > 10) {
        nonPromptMuonsRateAn.fill(l1MuonCand);
      }
      else {
        promptMuonsRateAn.fill(l1MuonCand);
      }
    }
  }
  else {
    notMatchedRateAn.fill(l1MuonCand);
  }

  likelihoodDistribution.fill(0, 0, 0, l1MuonCand);
}*/


void CandsMatchingAnalyser::write() {
  promptMuons.write();
  nonPromptMuons.write();
  muonsFromPions.write();
  muonsFromKaons.write();
  notMatchedRateAn.write();

  simVertexRhoVsPtGen->Write();
}
} /* namespace L1MuAn */
