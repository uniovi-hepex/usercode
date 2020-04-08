from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from libPyROOT import TDirectory


def makeEfficiency(passed, total, title, lineColor):
    if TEfficiency.CheckConsistency(passed, total) :
        efficiency = TEfficiency(passed, total)
        #title = std::regex_replace(title, std::regex("\\muCandGenEtaMuons"), "tagging efficiency");
        efficiency.SetTitle( title );
        efficiency.SetStatisticOption(6  ); #TEfficiency.EStatOption.kBUniform
        efficiency.SetPosteriorMode();
        efficiency.SetLineColor(lineColor);
        return efficiency;
    else :
        print("makeEfficiency TEfficiency::CheckConsistency(*ptGenPtTTMuonNom, *ptGenPtTTMuonDenom) failed" )
        exit(1);
    

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_rate_v0006.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v28_10Files.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t28/results/omtfAnalysis2.root' )

#histFile.ls()


lhcFillingRatio = 2760./3564.;
lhcFreq = 40144896; #11264 * 3564

analyzerOmtfDir = histFile.Get("L1MuonAnalyzerOmtf")
candPerEvent = analyzerOmtfDir.Get("candPerEvent")
print ("candPerEvent " + str(type(candPerEvent) ))

eventCnt = candPerEvent.Integral() 
scale = 1./eventCnt * lhcFreq * lhcFillingRatio;

print ("eventCnt " + str(eventCnt) );
print ("scale " + str(scale) );


rateDir = histFile.Get("L1MuonAnalyzerOmtf/rate")
rateDir.ls()

canvases = []
rateCumuls = []
efficienciesHist1 = []
efficienciesHist2 = []

def makeRatePlots(algoDir, lineColor) :
    print (algoDir.GetName())
    algoDir.ls()
    
    c1 = TCanvas('canvas_' + algoDir.GetName(), algoDir.GetName(), 200, 10, 700, 900)
    c1.Divide(2, 2)
    c1.cd(1).SetGridx()
    c1.cd(1).SetGridy()
    print ('created canvas ' + c1.GetName())

    ##########################
    
    for obj in algoDir.GetListOfKeys():
        if "candPt" in obj.GetName() :
            candPt = obj.ReadObj()
            candPt.SetLineColor(lineColor)
            candPt.Draw("")


            candPt.Sumw2(False);
            candPt_rateCumul = candPt.GetCumulative(False, "_rate");
            candPt_rateCumul.SetBinContent(1, 0);
            #candPt_rateCumul.Sumw2(False);
            candPt_rateCumul.Scale(scale) #TODO maybe it should be before Sumw2
            candPt_rateCumul.Sumw2(False);
    
            c1.cd(2).SetGridx()
            c1.cd(2).SetGridy()   
            candPt_rateCumul.SetLineColor(lineColor)
            candPt_rateCumul.Draw("")
            
            rateCumuls.append(candPt_rateCumul)
       
    canvases.append(c1) 
       

for iAlgo, obj in enumerate(rateDir.GetListOfKeys() ) :
    algoDir = obj.ReadObj()
    if isinstance(algoDir, TDirectory):
        #makeRatePlots(5)
        lineColor = 2
        if iAlgo ==  0:
            lineColor = 1
        if iAlgo == 1:
            lineColor = 4
        makeRatePlots(algoDir, lineColor)
        
for iAlgo, canvas in enumerate(canvases ) :
    if iAlgo >= 1 :
        canvas.cd(2)
        rateCumuls[0].DrawCopy("same hist")
        canvas.cd(2).Modified()
        canvas.Update()
         
#         canvas.cd(4)
#         efficienciesHist2[0].DrawCopy("same hist")
#         canvas.cd(4).Modified()
#         canvas.Update()

#pad1 = TPad( 'pad1', 'The pad with the function',  0.03, 0.62, 0.50, 0.92, 21 )
#pad2 = TPad( 'pad2', 'The pad with the histogram', 0.51, 0.62, 0.98, 0.92, 21 )
#pad3 = TPad( 'pad3', 'The pad with the histogram', 0.03, 0.02, 0.97, 0.57, 21 )
#pad1.Draw()
#pad2.Draw()
#pad3.Draw()

#%jsroot on
#from ROOT import gROOT 
gROOT.GetListOfCanvases().Draw()

#execfile('efficiencyPlots.py')
