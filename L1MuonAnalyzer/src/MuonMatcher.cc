/*
 * MuonMatcher.cc
 *
 *  Created on: Jun 16, 2020
 *      Author: kbunkow
 */

#include "usercode/L1MuonAnalyzer/interface/MuonMatcher.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

#include <math.h>

#include "TFile.h"


namespace L1MuAn {

double hwGmtPhiToGlobalPhi(int phi) {
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

MuonMatcher::MuonMatcher(const edm::ParameterSet& edmCfg) {
  edm::Service<TFileService> fs;
  TFileDirectory subDir =  fs->mkdir("MuonMatcher");

  deltaPhiPropCand   = subDir.make<TH2F>("deltaPhiPropCand",  "delta Phi propagated track - muonCand Vs Pt", 100, 0, 100, 100, -0.5 -0.005, 0.5 -0.005); //delta phi between propagated track and muon candidate,
  deltaPhiVertexProp = subDir.make<TH2F>("deltaPhiVertexProp", "delta Phi vertex - propagated track Vs Pt", 100, 0, 100, 100, -1, 1);; //delta phi between phi at vertex and propagated track phi)

  if(edmCfg.exists("muonMatcherFile") ) {
    std::string muonMatcherFileName =  edmCfg.getParameter<edm::FileInPath>("muonMatcherFile").fullPath();
    TFile inFile(muonMatcherFileName.c_str());
    edm::LogImportant("l1tMuBayesEventPrint") <<" MuonMatcher: using muonMatcherFileName "<<muonMatcherFileName<<std::endl;

    deltaPhiPropCandMean = (TH1D*)inFile.Get("deltaPhiPropCandMean");
    deltaPhiPropCandStdDev = (TH1D*)inFile.Get("deltaPhiPropCandStdDev");

    fillMean =  false;
  }
  else {
    deltaPhiPropCandMean = new TH1D("deltaPhiPropCandMean", "mean #DELTA#phi propagated track - muonCand Vs Pt", 100, 0, 100); //TODO the bins must be the same as in the deltaPhiPropCand
    deltaPhiPropCandStdDev = new TH1D("deltaPhiPropCandStdDev", "StdDev #delta#phi propagated track - muonCand Vs Pt", 100, 0, 100); //TODO the bins must be the same as in the deltaPhiPropCand
    fillMean =  true;

  }

}

void MuonMatcher::setup(const edm::EventSetup& eventSetup) {
  eventSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
  eventSetup.get<IdealMagneticFieldRecord>().get(magField);

  eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator);
}

MuonMatcher::~MuonMatcher() {


}

void MuonMatcher::saveHists() {
  if(fillMean) {
    std::string rootFileName = "muonMatcherHists.root";
    //gStyle->SetOptStat(111111);
    TFile outfile(rootFileName.c_str(), "RECREATE");
    edm::LogImportant("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<" out fileName "<<rootFileName<<" outfile->GetName() "<<outfile.GetName()<<endl;
    outfile.cd();
    deltaPhiPropCand->Write();
    deltaPhiVertexProp->Write();

    for(int iBin = 0; iBin <= deltaPhiPropCand->GetXaxis()->GetNbins() +1; iBin++) {
      auto projection = deltaPhiPropCand->ProjectionY( (std::string("deltaPhiPropCand") + std::to_string(iBin)).c_str(), iBin, iBin);

      double mean = 0;
      double stdDev = 0;
      double entries = projection->Integral(0, projection->GetNbinsX() +1); //to iniclude the under and overflow bins
      if(entries) {
        mean = deltaPhiPropCandMean->GetBinContent(iBin) / entries ;

        stdDev = deltaPhiPropCandStdDev->GetBinContent(iBin) / entries;
        stdDev = sqrt(stdDev - mean * mean);
      }
      deltaPhiPropCandMean->SetBinContent(iBin, mean);
      deltaPhiPropCandStdDev->SetBinContent(iBin, stdDev);

      edm::LogImportant("l1tMuBayesEventPrint") <<" MuonMatcher::saveHists() "<<std::setw(5)<<iBin<<" mean "<<std::setw(8)<<mean<<" GetMean "<<projection->GetMean()<<
          " stdDev "<<stdDev<<" GetRMS "<<projection->GetRMS()<<" Entries "<<entries<<std::endl;

    }

    deltaPhiPropCandMean->Write();
    deltaPhiPropCandStdDev->Write();

    TH1* deltaPhiPropCandStdDevSmooth = (TH1*)deltaPhiPropCandStdDev->Clone("deltaPhiPropCandStdDevSmooth");
    deltaPhiPropCandStdDevSmooth->GetXaxis()->SetRangeUser(3, 100);
    deltaPhiPropCandStdDevSmooth->Smooth(1, "R");
    deltaPhiPropCandStdDevSmooth->Write();
  }
}


TrajectoryStateOnSurface MuonMatcher::atStation2(FreeTrajectoryState ftsStart, float eta) const {
  ReferenceCountingPointer<Surface> rpc;
  if (eta < -1.24)       rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,-790.),  TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) );
  else if (eta < 1.24)   rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 500, 500, -900, 900 ) ) );
  else                   rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,790.),   TkRotation<float>(), SimpleDiskBounds( 300., 810., -10., 10. ) ) );
//  if (eta < -1.04)       rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,-790.), TkRotation<float>(), SimpleDiskBounds( 300., 710., -10., 10. ) ) );
//  else if (eta < -0.72)  rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 520, 520, -700, 700 ) ) );
//  else if (eta < 0.72)   rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 500, 500, -700, 700 ) ) );
//  else if (eta < 1.04)   rpc = ReferenceCountingPointer<Surface>(new  BoundCylinder( GlobalPoint(0.,0.,0.), TkRotation<float>(), SimpleCylinderBounds( 520, 520, -700, 700 ) ) );
//  else                      rpc = ReferenceCountingPointer<Surface>(new  BoundDisk( GlobalPoint(0.,0.,790.), TkRotation<float>(), SimpleDiskBounds( 300., 710., -10., 10. ) ) );
  //theEs.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator);
  TrajectoryStateOnSurface trackAtRPC =  propagator->propagate( ftsStart, *rpc);
  return trackAtRPC;
}

FreeTrajectoryState MuonMatcher::simTrackToFts(const SimTrack& simTrackPtr, const SimVertex& simVertex) {
  int charge = simTrackPtr.type() > 0 ? -1 : 1; //works for muons

  CLHEP::Hep3Vector p3T(simTrackPtr.momentum().x(), simTrackPtr.momentum().y(), simTrackPtr.momentum().z());
  //if (p3T.mag()< 2.) continue;

  CLHEP::Hep3Vector  r3T = CLHEP::Hep3Vector(simVertex.position().x(),
                            simVertex.position().y(),
                            simVertex.position().z());


  GlobalVector p3GV(p3T.x(), p3T.y(), p3T.z());
  GlobalPoint r3GP(r3T.x(), r3T.y(), r3T.z());

  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, &*magField);

  //CartesianTrajectoryError tCov(cov);

  //return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;

  return FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState MuonMatcher::simTrackToFts(const TrackingParticle& trackingParticle) {
  int charge = trackingParticle.pdgId() > 0 ? -1 : 1; //works for muons

  CLHEP::Hep3Vector p3T(trackingParticle.momentum().x(), trackingParticle.momentum().y(), trackingParticle.momentum().z());
  //if (p3T.mag()< 2.) continue;

  CLHEP::Hep3Vector  r3T = CLHEP::Hep3Vector(trackingParticle.vx(), trackingParticle.vy(), trackingParticle.vz());


  GlobalVector p3GV(p3T.x(), p3T.y(), p3T.z());
  GlobalPoint r3GP(r3T.x(), r3T.y(), r3T.z());

  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, &*magField);

  //CartesianTrajectoryError tCov(cov);

  //return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;

  return FreeTrajectoryState(tPars) ;
}

TrajectoryStateOnSurface MuonMatcher::propagate(const SimTrack& simTrack, const edm::SimVertexContainer* simVertices) {
  SimVertex simVertex;
  int vtxInd = simTrack.vertIndex();
  if (vtxInd < 0){
    std::cout<<"Track with no vertex, defaulting to (0,0,0)"<<std::endl;
  }
  else {
    simVertex = simVertices->at(vtxInd);
    if(((int)simVertex.vertexId()) != vtxInd) {
      std::cout<<"simVertex.vertexId() != vtxInd !!!!!!!!!!!!!!!!!"<<std::endl;
      edm::LogImportant("l1tMuBayesEventPrint") <<"simVertex.vertexId() != vtxInd. simVertex.vertexId() "<<simVertex.vertexId()<<" vtxInd "<<vtxInd<<" !!!!!!!!!!!!!!!!!";
    }
  }

  FreeTrajectoryState ftsTrack = simTrackToFts(simTrack, simVertex);

  TrajectoryStateOnSurface tsof = atStation2(ftsTrack, simTrack.momentum().eta() ); //propagation

  return tsof;
}


TrajectoryStateOnSurface MuonMatcher::propagate(const TrackingParticle& trackingParticle) {
  FreeTrajectoryState ftsTrack = simTrackToFts(trackingParticle);

  TrajectoryStateOnSurface tsof = atStation2(ftsTrack, trackingParticle.momentum().eta() ); //propagation

  return tsof;
}


float normal_pdf(float x, float m, float s) {
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}

MatchingResult MuonMatcher::match(const l1t::RegionalMuonCand* muonCand, const SimTrack& simTrack, TrajectoryStateOnSurface& tsof) {
  MatchingResult result;
  result.simTrack = &simTrack;
  double candGloablEta  = muonCand->hwEta() * 0.010875;
  if( abs(simTrack.momentum().eta() - candGloablEta ) < 0.3 ) {
    double candGlobalPhi = l1t::MicroGMTConfiguration::calcGlobalPhi( muonCand->hwPhi(), muonCand->trackFinderType(), muonCand->processor() );
    candGlobalPhi = hwGmtPhiToGlobalPhi(candGlobalPhi );

    if(candGlobalPhi > M_PI)
      candGlobalPhi = candGlobalPhi -(2.*M_PI);

    result.deltaPhi = foldPhi(tsof.globalPosition().phi() - candGlobalPhi);
    result.deltaEta = tsof.globalPosition().eta() - candGloablEta;

    double mean = 0;
    double sigma = 1;
    if(!fillMean) {
      auto ptBin = deltaPhiPropCandMean->FindBin(simTrack.momentum().pt());

      mean = deltaPhiPropCandMean->GetBinContent(ptBin);
      sigma = deltaPhiPropCandStdDev->GetBinContent(ptBin);
    }

    result.matchingLikelihood = normal_pdf(result.deltaPhi, mean, sigma); //TODO temporary solution

    result.muonCand = muonCand;

    if( abs(result.deltaPhi) < (4. * sigma)) //TODO 4 sigma, because the distribution has non-gaussian tails
      result.result = MatchingResult::ResultType::matched;

    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: simTrack type "<<simTrack.type()<<" pt "<<std::setw(8)<<simTrack.momentum().pt()
        <<" eta "<<std::setw(8)<<simTrack.momentum().eta()<<" phi "<<std::setw(8)<<simTrack.momentum().phi()
        <<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
        <<" muonCand pt "<<std::setw(8)<<muonCand->hwPt()<<" candGloablEta "<<std::setw(8)<<candGloablEta<<" candGlobalPhi "<<std::setw(8)<<candGlobalPhi<<" hwQual "<<muonCand->hwQual()
        <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
        << std::endl;

/*    if(abs(result.deltaPhi) > 0.4)
      edm::LogImportant("l1tMuBayesEventPrint") <<"MuonMatcher::match: simTrack type "<<simTrack.type()<<" pt "<<std::setw(8)<<simTrack.momentum().pt()
          <<" eta "<<std::setw(8)<<simTrack.momentum().eta()<<" phi "<<std::setw(8)<<simTrack.momentum().phi()
          <<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
          <<" muonCand pt "<<std::setw(8)<<muonCand->hwPt()<<" candGloablEta "<<std::setw(8)<<candGloablEta<<" candGlobalPhi "<<std::setw(8)<<candGlobalPhi<<" hwQual "<<muonCand->hwQual()
          <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
          << std::endl;*/

  }

  return result;
}

MatchingResult MuonMatcher::match(const l1t::RegionalMuonCand* muonCand, const TrackingParticle& trackingParticle, TrajectoryStateOnSurface& tsof) {
  MatchingResult result;
  result.trackingParticle = &trackingParticle;
  double candGloablEta  = muonCand->hwEta() * 0.010875;
  if( abs(trackingParticle.momentum().eta() - candGloablEta ) < 0.3 ) {
    double candGlobalPhi = l1t::MicroGMTConfiguration::calcGlobalPhi( muonCand->hwPhi(), muonCand->trackFinderType(), muonCand->processor() );
    candGlobalPhi = hwGmtPhiToGlobalPhi(candGlobalPhi );

    if(candGlobalPhi > M_PI)
      candGlobalPhi = candGlobalPhi -(2.*M_PI);

    result.deltaPhi = foldPhi(tsof.globalPosition().phi() - candGlobalPhi);
    result.deltaEta = tsof.globalPosition().eta() - candGloablEta;

    double mean = 0;
    double sigma = 1;
    if(!fillMean) {
      auto ptBin = deltaPhiPropCandMean->FindBin(trackingParticle.pt());

      mean = deltaPhiPropCandMean->GetBinContent(ptBin);
      sigma = deltaPhiPropCandStdDev->GetBinContent(ptBin);
    }

    result.matchingLikelihood = normal_pdf(result.deltaPhi, mean, sigma); //TODO temporary solution

    result.muonCand = muonCand;

    if( abs(result.deltaPhi) < (4. * sigma)) //TODO 4 sigma, because the distribution has non-gaussian tails
      result.result = MatchingResult::ResultType::matched;

    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: simTrack type "<<trackingParticle.pdgId()<<" pt "<<std::setw(8)<<trackingParticle.pt()
        <<" eta "<<std::setw(8)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(8)<<trackingParticle.momentum().phi()
        <<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
        <<" muonCand pt "<<std::setw(8)<<muonCand->hwPt()<<" candGloablEta "<<std::setw(8)<<candGloablEta<<" candGlobalPhi "<<std::setw(8)<<candGlobalPhi<<" hwQual "<<muonCand->hwQual()
        <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
        << std::endl;

/*    if(abs(result.deltaPhi) > 0.4)
      edm::LogImportant("l1tMuBayesEventPrint") <<"MuonMatcher::match: simTrack type "<<simTrack.type()<<" pt "<<std::setw(8)<<simTrack.momentum().pt()
          <<" eta "<<std::setw(8)<<simTrack.momentum().eta()<<" phi "<<std::setw(8)<<simTrack.momentum().phi()
          <<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
          <<" muonCand pt "<<std::setw(8)<<muonCand->hwPt()<<" candGloablEta "<<std::setw(8)<<candGloablEta<<" candGlobalPhi "<<std::setw(8)<<candGlobalPhi<<" hwQual "<<muonCand->hwQual()
          <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
          << std::endl;*/

  }

  return result;
}

std::vector<MatchingResult> MuonMatcher::match(std::vector<const l1t::RegionalMuonCand*>& muonCands, const edm::SimTrackContainer* simTracks, const edm::SimVertexContainer* simVertices,
    std::function<bool(const SimTrack& )> const& simTrackFilter)
{
  std::vector<MatchingResult> matchingResults;

  for (auto& simTrack : *simTracks ) {

    if(!simTrackFilter(simTrack))
      continue;

    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match, simTrack type "<<std::setw(3)<<simTrack.type()<<" pt "<<std::setw(9)<<simTrack.momentum().pt()<<" eta "<<std::setw(9)<<simTrack.momentum().eta()<<" phi "<<std::setw(9)<<simTrack.momentum().phi()<<std::endl;

    bool matched = false;

    TrajectoryStateOnSurface tsof = propagate(simTrack, simVertices);
    if(!tsof.isValid()) {
      LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" propagation failed"<<std::endl;
      MatchingResult result;
      result.result = MatchingResult::ResultType::propagationFailed;
      continue; //no sense to do matching
    }


    double ptGen = simTrack.momentum().pt();
    if(ptGen >= deltaPhiVertexProp->GetXaxis()->GetXmax())
      ptGen = deltaPhiVertexProp->GetXaxis()->GetXmax() - 0.01;


    deltaPhiVertexProp->Fill(ptGen, simTrack.momentum().phi() - tsof.globalPosition().phi());


    for(auto& muonCand : muonCands) {
      //int refLayer = (int)omtfCand->trackAddress().at(1);
      //int layerHits = (int)omtfCand->trackAddress().at(0);
      //std::bitset<18> layerHitBits(layerHits);

      if(muonCand->hwQual() <= 1) //dropping very low quality candidates, as they are fakes usually
        continue;

      MatchingResult result = match(muonCand, simTrack, tsof);
      if(result.result == MatchingResult::ResultType::matched) {
        matchingResults.push_back(result);
        matched = true;
      }
    }

    if(!matched) { //we are adding also if it was not matching to any candidate
      MatchingResult result;
      result.simTrack = &simTrack;
      matchingResults.push_back(result);
      LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" no matching candidate found"<<std::endl;
    }
  }

  //Cleaning the matching
  std::sort(matchingResults.begin(), matchingResults.end(),
      [](const MatchingResult& a, const MatchingResult& b)->bool { return a.matchingLikelihood > b.matchingLikelihood; } );

  for(unsigned int i1 = 0; i1 < matchingResults.size(); i1++) {
    if(matchingResults[i1].result == MatchingResult::ResultType::matched) {
      for(unsigned int i2 = i1 + 1; i2 < matchingResults.size(); i2++) {
        if( (matchingResults[i1].simTrack == matchingResults[i2].simTrack) || (matchingResults[i1].muonCand == matchingResults[i2].muonCand) ) {
          //if matchingResults[i1].muonCand == false, then it is also OK here
          matchingResults[i2].result = MatchingResult::ResultType::duplicate;
        }
      }
    }
  }

  std::vector<MatchingResult> cleanedMatchingResults;
  for(auto& matchingResult : matchingResults) {
    if(matchingResult.result == MatchingResult::ResultType::matched  || matchingResult.muonCand == nullptr) //adding also the simTracks that are not matched at all, before it is assured that they are not duplicates
      cleanedMatchingResults.push_back(matchingResult);

    if(matchingResult.result == MatchingResult::ResultType::matched) {
      double ptGen = matchingResult.simTrack->momentum().pt();
      if(ptGen >= deltaPhiPropCand->GetXaxis()->GetXmax())
        ptGen = deltaPhiPropCand->GetXaxis()->GetXmax() - 0.01;

      deltaPhiPropCand->Fill(ptGen, matchingResult.deltaPhi);

      if(fillMean) {
        deltaPhiPropCandMean->Fill(matchingResult.simTrack->momentum().pt(), matchingResult.deltaPhi); //filling oveflow is ok here
        deltaPhiPropCandStdDev->Fill(matchingResult.simTrack->momentum().pt(), matchingResult.deltaPhi * matchingResult.deltaPhi);
      }
    }
  }

  //adding the muonCand-s that were not matched, i.e. in order to analyze them later
  for(auto& muonCand : muonCands) {
    bool isMatched = false;
    for(auto& matchingResult : cleanedMatchingResults) {
      if(matchingResult.muonCand == muonCand) {
        isMatched =  true;
        break;
      }
    }

    if(!isMatched) {
      MatchingResult result;
      result.muonCand = muonCand;
      cleanedMatchingResults.push_back(result);
    }
  }


  LogTrace("l1tMuBayesEventPrint")<<"MuonMatcher::match cleanedMatchingResults:"<<std::endl;
  for(auto& result : cleanedMatchingResults) {
    if(result.simTrack)
      LogTrace("l1tMuBayesEventPrint") <<"simTrack type "<<result.simTrack->type()<<" pt "<<std::setw(8)<<result.simTrack->momentum().pt()
        <<" eta "<<std::setw(8)<<result.simTrack->momentum().eta()<<" phi "<<std::setw(8)<<result.simTrack->momentum().phi();
    else
      LogTrace("l1tMuBayesEventPrint")<<"no sim track ";

        //<<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
    if(result.muonCand)
      LogTrace("l1tMuBayesEventPrint")<<" muonCand pt "<<std::setw(8)<<result.muonCand->hwPt()<<" hwQual "<<result.muonCand->hwQual()<<" hwEta "<<result.muonCand->hwEta()
        <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
        << std::endl;
    else
      LogTrace("l1tMuBayesEventPrint")<<" no muonCand "<<" result "<<(short)result.result<< std::endl;
  }
  LogTrace("l1tMuBayesEventPrint")<<" "<<std::endl;

  return cleanedMatchingResults;
}



std::vector<MatchingResult> MuonMatcher::match(std::vector<const l1t::RegionalMuonCand*>& muonCands, const TrackingParticleCollection* trackingParticles,
    std::function<bool(const TrackingParticle& )> const& simTrackFilter)
{
  std::vector<MatchingResult> matchingResults;
  LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match trackingParticles->size() "<<trackingParticles->size()<<std::endl;

  for (auto& trackingParticle : *trackingParticles ) {
    //LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match:"<<__LINE__<<" trackingParticle type "<<std::setw(3)<<trackingParticle.pdgId()<<" pt "<<std::setw(9)<<trackingParticle.pt()<<" eta "<<std::setw(9)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(9)<<trackingParticle.momentum().phi()<<std::endl;

    if(simTrackFilter(trackingParticle) == false)
      continue;

    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match, trackingParticle type "<<std::setw(3)<<trackingParticle.pdgId()<<" pt "<<std::setw(9)<<trackingParticle.pt()<<" eta "<<std::setw(9)<<trackingParticle.momentum().eta()<<" phi "<<std::setw(9)<<trackingParticle.momentum().phi()<<std::endl;

    bool matched = false;

    TrajectoryStateOnSurface tsof = propagate(trackingParticle);
    if(!tsof.isValid()) {
      LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" propagation failed"<<std::endl;
      MatchingResult result;
      result.result = MatchingResult::ResultType::propagationFailed;
      continue; //no sense to do matching
    }
    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;

    double ptGen = trackingParticle.pt();
    if(ptGen >= deltaPhiVertexProp->GetXaxis()->GetXmax())
      ptGen = deltaPhiVertexProp->GetXaxis()->GetXmax() - 0.01;


    deltaPhiVertexProp->Fill(ptGen, trackingParticle.momentum().phi() - tsof.globalPosition().phi());
    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;

    for(auto& muonCand : muonCands) {
      //int refLayer = (int)omtfCand->trackAddress().at(1);
      //int layerHits = (int)omtfCand->trackAddress().at(0);
      //std::bitset<18> layerHitBits(layerHits);

      if(muonCand->hwQual() <= 1) //dropping very low quality candidates, as they are fakes usually
        continue;
      LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
      MatchingResult result = match(muonCand, trackingParticle, tsof);
      LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
      if(result.result == MatchingResult::ResultType::matched) {
        matchingResults.push_back(result);
        matched = true;
      }
    }

    if(!matched) { //we are adding also if it was not matching to any candidate
      MatchingResult result;
      result.trackingParticle = &trackingParticle;
      matchingResults.push_back(result);
      LogTrace("l1tMuBayesEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" no matching candidate found"<<std::endl;
    }
  }
  LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
  //Cleaning the matching
  std::sort(matchingResults.begin(), matchingResults.end(),
      [](const MatchingResult& a, const MatchingResult& b)->bool { return a.matchingLikelihood > b.matchingLikelihood; } );
  LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
  for(unsigned int i1 = 0; i1 < matchingResults.size(); i1++) {
    if(matchingResults[i1].result == MatchingResult::ResultType::matched) {
      for(unsigned int i2 = i1 + 1; i2 < matchingResults.size(); i2++) {
        if( (matchingResults[i1].trackingParticle == matchingResults[i2].trackingParticle) || (matchingResults[i1].muonCand == matchingResults[i2].muonCand) ) {
          //if matchingResults[i1].muonCand == false, then it is also OK here
          matchingResults[i2].result = MatchingResult::ResultType::duplicate;
        }
      }
    }
  }
  LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
  std::vector<MatchingResult> cleanedMatchingResults;
  for(auto& matchingResult : matchingResults) {
    if(matchingResult.result == MatchingResult::ResultType::matched  || matchingResult.muonCand == nullptr) //adding also the simTracks that are not matched at all, before it is assured that they are not duplicates
      cleanedMatchingResults.push_back(matchingResult);
    LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
    if(matchingResult.result == MatchingResult::ResultType::matched) {
      double ptGen = matchingResult.trackingParticle->pt();
      if(ptGen >= deltaPhiPropCand->GetXaxis()->GetXmax())
        ptGen = deltaPhiPropCand->GetXaxis()->GetXmax() - 0.01;

      deltaPhiPropCand->Fill(ptGen, matchingResult.deltaPhi);
      LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
      if(fillMean) {
        deltaPhiPropCandMean->Fill(matchingResult.trackingParticle->pt(), matchingResult.deltaPhi); //filling oveflow is ok here
        deltaPhiPropCandStdDev->Fill(matchingResult.trackingParticle->pt(), matchingResult.deltaPhi * matchingResult.deltaPhi);
      }
      LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
    }
  }
  LogTrace("l1tMuBayesEventPrint") <<"MuonMatcher::match: "<<__LINE__<<std::endl;
  //adding the muonCand-s that were not matched, i.e. in order to analyze them later
  for(auto& muonCand : muonCands) {
    bool isMatched = false;
    for(auto& matchingResult : cleanedMatchingResults) {
      if(matchingResult.muonCand == muonCand) {
        isMatched =  true;
        break;
      }
    }

    if(!isMatched) {
      MatchingResult result;
      result.muonCand = muonCand;
      cleanedMatchingResults.push_back(result);
    }
  }


  LogTrace("l1tMuBayesEventPrint")<<":"<<__LINE__<<" MuonMatcher::match cleanedMatchingResults:"<<std::endl;
  for(auto& result : cleanedMatchingResults) {
    if(result.trackingParticle)
      LogTrace("l1tMuBayesEventPrint")<<":"<<__LINE__ <<" simTrack type "<<result.trackingParticle->pdgId()<<" pt "<<std::setw(8)<<result.trackingParticle->pt()
        <<" eta "<<std::setw(8)<<result.trackingParticle->momentum().eta()<<" phi "<<std::setw(8)<<result.trackingParticle->momentum().phi();
    else
      LogTrace("l1tMuBayesEventPrint")<<"no sim track ";

        //<<" propagation eta "<<std::setw(8)<<tsof.globalPosition().eta()<<" phi "<<tsof.globalPosition().phi()
    if(result.muonCand)
      LogTrace("l1tMuBayesEventPrint")<<" muonCand pt "<<std::setw(8)<<result.muonCand->hwPt()<<" hwQual "<<result.muonCand->hwQual()<<" hwEta "<<result.muonCand->hwEta()
        <<" deltaEta "<<std::setw(8)<<result.deltaEta<<" deltaPhi "<<std::setw(8)<<result.deltaPhi<<" Likelihood "<<std::setw(8)<<result.matchingLikelihood <<" result "<<(short)result.result
        << std::endl;
    else
      LogTrace("l1tMuBayesEventPrint")<<" no muonCand "<<" result "<<(short)result.result<< std::endl;
  }
  LogTrace("l1tMuBayesEventPrint")<<" "<<std::endl;

  return cleanedMatchingResults;
}

}
