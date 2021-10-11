/*
 * L1MuonCand.h
 *
 *  Created on: Mar 18, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_L1MUONCAND_H_
#define INTERFACE_L1MUONCAND_H_


#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

namespace L1MuAn {

/*
 * wrapper for using in L1MuonAnalyzerOmtf
 */
class L1MuonCand {
public:

  L1MuonCand();

  L1MuonCand(const l1t::RegionalMuonCand& muCand):
    ptGev( muCand.hwPt() > 0 ? (muCand.hwPt()-1)/2. : 0), //TODO watch out!!!!
    uptGev(muCand.hwPtUnconstrained() > 0 ? (muCand.hwPtUnconstrained()-1)/2. : 0), //TODO watch out!!!!

    hwPt(muCand.hwPt()),
    hwUPt(muCand.hwPtUnconstrained()),

    hwEta(muCand.hwEta()),
    hwPhi(muCand.hwPhi()),
    hwQual(muCand.hwQual()),
    hwSign(muCand.hwSign() ),
    hwBeta(0),
    firedLayerBits(muCand.trackAddress().at(0)),
    refLayer(muCand.trackAddress().at(1)),
    likelihood(muCand.trackAddress().at(2)),
    likelihoodUpt(muCand.trackAddress().at(3))
    //firedLayerCnt(firedLayerCnt)
    {
    }

  virtual ~L1MuonCand();

  // integer "hardware" values
  double ptGev = 0;
  double uptGev = 0;

  int hwPt = 0;
  int hwUPt = 0;

  int hwEta = 0;
  int hwPhi = 0;
  int hwQual = 0;

  int hwSign = 0;

  int hwBeta = 0;
  //float beta = 0;
  //float betaLikelihood = 0;

  //boost::dynamic_bitset<> firedLayerBits;
  unsigned int firedLayerBits = 0;
  //unsigned int firedLayerCnt = 0;
  //int pdfSum = 0;
  int refLayer = 0;
  float likelihood = 0;
  float likelihoodUpt = 0;
};

} /* namespace L1MuAn */

#endif /* INTERFACE_L1MUONCAND_H_ */
