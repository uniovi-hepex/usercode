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
  histTitle<<name<<" cand pt, "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candPt = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut1GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 1 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut1 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut10GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 10 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

  candEtaPtCut10 = subDir.make<TH1D>(histName.str().c_str(), histTitle.str().c_str(), 100, -2.1, 2.1);

  histName.str("");
  histTitle.str("");

  histName<<"candEta_PtCut20GeV_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" cand Eta, ptCut 20 GeV "<<" quality >= "<<qualityCut<<";pt [GeV]; #events";

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

CandsMatchingAnalyser::CandsMatchingAnalyser(TFileDirectory& subDir, std::string name, int qualityCut, int nBins, double binsFrom, double binsTo):
    promptMuons(subDir.mkdir("promptMuons"), "promptMuons", qualityCut, nBins, binsFrom, binsTo),
    nonPromptMuons(subDir.mkdir("nonPromptMuons"), "nonPromptMuons", qualityCut, nBins, binsFrom, binsTo),
    muonsFromPions(subDir.mkdir("muonsFromPions"), "muonsFromPions", qualityCut, nBins, binsFrom, binsTo),
    muonsFromKaons(subDir.mkdir("muonsFromKaons"), "muonsFromKaons", qualityCut, nBins, binsFrom, binsTo),
    notMatched(subDir.mkdir("notMatched"), "notMatched", qualityCut, nBins, binsFrom, binsTo),
    likelihoodDistribution(subDir, name, qualityCut, 0, 20, 120)
{

  std::ostringstream histName;
  std::ostringstream histTitle;

  histName<<name<<"simVertexRhoVsCandPt_"<<"_qualityCut_"<<qualityCut;
  histTitle<<name<<" simVertexRho Vs CandPt, "<<" quality >= "<<qualityCut<<";pt [GeV]; Rho [cm]";

  simVertexRhoVsCandPt = subDir.make<TH2I>(histName.str().c_str(), histTitle.str().c_str(), nBins, binsFrom, binsTo, 30, 0, 300);

}

void CandsMatchingAnalyser::fill(L1MuonCand& l1MuonCand, const SimTrack* matchedSimTrack, const edm::SimVertexContainer* simVertices) {
  if(matchedSimTrack) {
    int vtxInd = matchedSimTrack->vertIndex();
    if (vtxInd < 0) {
      promptMuons.fill(l1MuonCand);
      simVertexRhoVsCandPt->Fill(l1MuonCand.ptGev, 0.);
    }
    else {
      const SimVertex& simVertex = simVertices->at(vtxInd);

      //simVertex.parentIndex();
      //matchedSimTrack->trackId();

      double simVertexRho = simVertex.position().Rho();

      simVertexRhoVsCandPt->Fill(l1MuonCand.ptGev, simVertexRho);

      if(simVertexRho > 10) {
        nonPromptMuons.fill(l1MuonCand);
      }
      else {
        promptMuons.fill(l1MuonCand);
      }
    }
  }
  else {
    notMatched.fill(l1MuonCand);
  }

  likelihoodDistribution.fill(0, 0, 0, l1MuonCand);
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
      promptMuons.fill(l1MuonCand);
      simVertexRhoVsCandPt->Fill(l1MuonCand.ptGev, 0.);
    }
    else {
      double simVertexRho = matchedTrackingParticle->parentVertex()->position().Rho();

      simVertexRhoVsCandPt->Fill(l1MuonCand.ptGev, simVertexRho);

      if(simVertexRho > 10) {
        nonPromptMuons.fill(l1MuonCand);
      }
      else {
        promptMuons.fill(l1MuonCand);
      }

      if(abs(parentTrackPdgId) == 211)
        muonsFromPions.fill(l1MuonCand);
      else if(abs(parentTrackPdgId) == 321)
        muonsFromKaons.fill(l1MuonCand);
    }
  }
  else {
    notMatched.fill(l1MuonCand);
  }

  likelihoodDistribution.fill(0, 0, 0, l1MuonCand);
}

void CandsMatchingAnalyser::write() {
  promptMuons.write();
  nonPromptMuons.write();
  muonsFromPions.write();
  muonsFromKaons.write();
  notMatched.write();

  simVertexRhoVsCandPt->Write();
}
} /* namespace L1MuAn */
