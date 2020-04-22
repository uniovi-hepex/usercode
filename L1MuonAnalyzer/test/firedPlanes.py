from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from libPyROOT import TDirectory
import sys
    

#version = "v2_t" + sys.argv[0]
version = "v2_t37" 

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_rate_v0006.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v31.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t30/results/omtfAnalysis2.root' )
rateHistFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_' + version + '/results/omtfAnalysis2.root' )
effHistFile =  TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_' + version + '/results/omtfAnalysis2.root' )


#histFile.ls()

print (rateHistFile)

lhcFillingRatio = 2760./3564.;
lhcFreq = 40144896; #11264 * 3564

analyzerOmtfDirRate = rateHistFile.Get("L1MuonAnalyzerOmtf")
candPerEvent = analyzerOmtfDirRate.Get("candPerEvent")
print ("candPerEvent " + str(type(candPerEvent) ))

eventCnt = candPerEvent.Integral() 
scale = 1./eventCnt * lhcFreq * lhcFillingRatio;

print ("eventCnt " + str(eventCnt) );
print ("scale " + str(scale) );


#firedPlanesEventCntOmtfRate = analyzerOmtfDirRate.Get("firedPlanesEventCntOmtf")
firedPlanesEventCntOmtfRate = analyzerOmtfDirRate.Get("firedPlanesEventCntNN")

firedPlanesEventCntOmtfEff = effHistFile.Get("L1MuonAnalyzerOmtf").Get("firedPlanesEventCntOmtf")
effNorm = firedPlanesEventCntOmtfEff.Integral() 

firedPlanesStat = []

fullRate = 0
for firedPlanes in range(0, firedPlanesEventCntOmtfRate.GetNbinsX(), 1) : 
    rate = firedPlanesEventCntOmtfRate.GetBinContent(firedPlanes +1) * scale
    eff = firedPlanesEventCntOmtfEff.GetBinContent(firedPlanes +1) / effNorm
    
    fullRate += rate
    eff_rate = 10000000
    if rate > 0 :
        eff_rate = eff/rate

    #if eff > 0 :
    #    eff_rate = rate/eff #eff/rate

    if rate > 0 or eff > 0:
        #if rate > 100:
            firedPlanesStat.append( (firedPlanes, rate, eff, eff_rate) )
            print("%8i %018i %f" % (firedPlanes, firedPlanes,  rate) ) 

firedPlanesStat.sort(key = lambda x: x[3], reverse = False)  

totalRateDrop = 0
totalEff = 0
for firedPlaneStat in firedPlanesStat :
    #print (format(firedPlaneStat[0], '018b'), firedPlaneStat)
    totalRateDrop += firedPlaneStat[1]
    totalEff  += firedPlaneStat[2]
    #if firedPlaneStat[1] > 0: #rate > 0
    print("%8i %s rate: %8.1f eff: %.5f ratio %f totalEff %f totalRateDrop %f fullRate %f " % (firedPlaneStat[0], format(firedPlaneStat[0], '018b'),  firedPlaneStat[1], firedPlaneStat[2], firedPlaneStat[3], totalEff, totalRateDrop, fullRate - totalRateDrop) ) 
    #print("%s" % (format(firedPlaneStat[0], '018b')) ) 

#execfile('ratePlots.py')
