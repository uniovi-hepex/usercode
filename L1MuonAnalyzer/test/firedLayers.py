from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from libPyROOT import TDirectory
import sys
    

#version = "v2_t" + sys.argv[0]
version = "v2_t52" 

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_rate_v0006.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v31.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t30/results/omtfAnalysis2.root' )
rateHistFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_SingleNeutrino_PU250_'+ version + '/results/omtfAnalysis2.root' )

#version = "v2_t44" 
effHistFile =  TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_' + version + '/results/omtfAnalysis2.root' )


#histFile.ls()

print ("rateHistFile " + rateHistFile.GetName())
print ("effHistFile " + effHistFile.GetName())

lhcFillingRatio = 2760./3564.;
lhcFreq = 40144896; #11264 * 3564

analyzerOmtfDirRate = rateHistFile.Get("L1MuonAnalyzerOmtf")
candPerEvent = analyzerOmtfDirRate.Get("candPerEvent")
print ("candPerEvent " + str(type(candPerEvent) ))

eventCnt = candPerEvent.Integral() 
scale = 1./eventCnt * lhcFreq * lhcFillingRatio;

print ("eventCnt " + str(eventCnt) );
print ("scale " + str(scale) );


firedLayersEventCntOmtfRate = analyzerOmtfDirRate.Get("firedLayersEventCntOmtf")
#firedLayersEventCntOmtfRate = analyzerOmtfDirRate.Get("firedLayersEventCntNN") #firedLayersEventCntNN firedPlanesEventCntNN

firedLayersEventCntOmtfEff = effHistFile.Get("L1MuonAnalyzerOmtf").Get("firedLayersEventCntOmtf")
#firedLayersEventCntOmtfEff = effHistFile.Get("L1MuonAnalyzerOmtf").Get("firedLayersEventCntNN") #firedLayersEventCntOmtf firedPlanesEventCntOmtf

print("rate hist " + firedLayersEventCntOmtfRate.GetName() )
print("eff  hist " + firedLayersEventCntOmtfEff.GetName() )



effNorm = firedLayersEventCntOmtfEff.Integral() 

firedLayersStat = []

fullRate = 0
for firedLayers in range(0, firedLayersEventCntOmtfRate.GetNbinsX(), 1) : 
    rate = firedLayersEventCntOmtfRate.GetBinContent(firedLayers +1) * scale
    eff = firedLayersEventCntOmtfEff.GetBinContent(firedLayers +1) / effNorm
    
    fullRate += rate
    eff_rate = 10000000
    if rate > 0 :
        eff_rate = eff/rate

    #if eff > 0 :
    #    eff_rate = rate/eff #eff/rate

    if rate > 0 or eff > 0:
        #if rate > 100:
        firedLayersStat.append( (firedLayers, rate, eff, eff_rate) )
            #print("%8i %018i %f" % (firedLayers, firedLayers,  rate) ) 
        #print("%8i %s rate: %8.1f eff: %.5f ratio %f" % (firedLayers, format(firedLayers, '018b'), rate, eff, eff_rate) ) 

print("\nselected\n")
firedLayersStat.sort(key = lambda x: x[3], reverse = False)  

totalRateDrop = 0
totalEff = 0
for firedLayerStat in firedLayersStat :
    #print (format(firedLayerStat[0], '018b'), firedLayerStat)
    if (firedLayerStat[1] > -1) :
    #if (firedLayerStat[1] > 150) or (firedLayerStat[2] < 0.0001 and firedLayerStat[1] > 100): #rate > 100
        totalRateDrop += firedLayerStat[1]
        totalEff  += firedLayerStat[2]
        print("%8i %s rate: %8.1f eff: %.5f ratio %f totalEff %f totalRateDrop %f fullRate %f " % (firedLayerStat[0], format(firedLayerStat[0], '018b'),  firedLayerStat[1], firedLayerStat[2], firedLayerStat[3], totalEff, totalRateDrop, fullRate - totalRateDrop) ) 
        #print("%s" % (format(firedLayerStat[0], '018b')) ) 

#execfile('ratePlots.py')
