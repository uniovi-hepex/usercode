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
    

version = "v2_t" + sys.argv[0]

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_rate_v0006.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v31.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t30/results/omtfAnalysis2.root' )
histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_' + version + '/results/omtfAnalysis2.root' )


#histFile.ls()

print (histFile)

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
        
        
ratesOnThreshHist = TH1D("ratesOnThreshHist", "ratesOnThreshHist", rateCumuls.__len__(), 0, rateCumuls.__len__()) 
relativeRatesOnThreshHist = TH1D("relativeRatesOnThreshHist", "relativeRatesOnThreshHist", rateCumuls.__len__(), 0, rateCumuls.__len__()) 
        
for iAlgo, canvas in enumerate(canvases ) :
    if iAlgo >= 1 :
        canvas.cd(2)
        rateCumuls[0].DrawCopy("same hist")
        canvas.cd(2).Modified()
        canvas.Update()
    
    ptCutGev = 21.5
    if iAlgo < 2 :
        ptCutGev = 20
    
    ptCutBin = rateCumuls[iAlgo].GetXaxis().FindBin(ptCutGev)        
    rateOnThresh = rateCumuls[iAlgo].GetBinContent(ptCutBin)   
    ratesOnThreshHist.Fill( iAlgo, rateOnThresh)
    ratesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() )   
    

    relativeRatesOnThresh = rateOnThresh / ratesOnThreshHist.GetBinContent(1)

    relativeRatesOnThreshHist.Fill( iAlgo, relativeRatesOnThresh)
    relativeRatesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() + " ptCut " + str(ptCutGev) + "GeV")    

canvasComapre = TCanvas('canvasComapre' , "compare " + version, 200, 10, 1400, 700)    
canvasComapre.Divide(2, 1)

canvasComapre.cd(1)
   
canvasComapre.cd(1).SetLeftMargin(0.4)
canvasComapre.cd(1).SetGridx()
canvasComapre.cd(1).SetGridy()
ratesOnThreshHist.GetYaxis().SetRangeUser(8000, 15000)
ratesOnThreshHist.GetYaxis().SetLabelSize(0.02)
ratesOnThreshHist.SetFillColor(0)
ratesOnThreshHist.SetFillStyle(3001)
ratesOnThreshHist.Draw("HBAR")

canvasComapre.cd(2)
   
canvasComapre.cd(2).SetLeftMargin(0.4)
canvasComapre.cd(2).SetGridx()
canvasComapre.cd(2).SetGridy()
relativeRatesOnThreshHist.GetYaxis().SetRangeUser(0.5, 1)
relativeRatesOnThreshHist.GetYaxis().SetLabelSize(0.02)
relativeRatesOnThreshHist.SetFillColor(0)
relativeRatesOnThreshHist.SetFillStyle(3001)
relativeRatesOnThreshHist.Draw("HBAR")

#%jsroot on
#from ROOT import gROOT 
gROOT.GetListOfCanvases().Draw()

#execfile('ratePlots.py')
