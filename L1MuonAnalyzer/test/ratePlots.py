from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle

from libPyROOT import TDirectory
import os
import sys

gStyle.SetOptStat(0)

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
    

#version = "PU200_v2_t" + sys.argv[1] #PU200_mtd5_v2_t
version = sys.argv[1] 
inputResults = 'SingleNeutrino_' + version 

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_rate_v0006.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v31.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t30/results/omtfAnalysis2.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_' + version + '/results/omtfAnalysis2.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )

if version < "MuFlatPt_PU200_v3_t73" :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )
else :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/L1Trigger/L1TMuonOverlapPhase1/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )


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

rateCumul_withTEffs = []
paintedGraphs = []
#gStyle.SetOptStat(111111111)

if not os.path.exists(inputResults):
    os.mkdir(inputResults)
    
outFile = TFile(inputResults + "/ratePlots.root", "RECREATE")




def makeRatePlotWithEff(candPt_rateCumul_copy, lineColor, canvasRate1, canvasRate2) :
    allEventsHist = candPt_rateCumul_copy.Clone(candPt_rateCumul_copy.GetName() + "_allEventsHist");
    for iBin in range(0, allEventsHist.GetNbinsX() ) :
        allEventsHist.SetBinContent(iBin, eventCnt)
    
    candPt_rateCumul_copy.Sumw2(False);
    allEventsHist.Sumw2(False);
    
    title = ("; ttTrack p_{T} [GeV]; rate [kHz]");
    rateCumul_withTEff = makeEfficiency(candPt_rateCumul_copy, allEventsHist, title, lineColor)
    
    canvasRate1.cd()
    rateCumul_withTEff.Draw("APZ")
    canvasRate1.Update()
    rateCumul_withTEff.GetPaintedGraph().GetXaxis().SetRangeUser(0, 100)
    
    paintedGraph = rateCumul_withTEff.GetPaintedGraph().Clone(rateCumul_withTEff.GetName() + "_copy" ) 
    scalekHz = 0.001
    for i in range(0, paintedGraph.GetN()) :
        paintedGraph.GetY()[i] *= lhcFreq * lhcFillingRatio * scalekHz
        paintedGraph.GetEYhigh()[i] *= lhcFreq * lhcFillingRatio * scalekHz
        paintedGraph.GetEYlow()[i] *= lhcFreq * lhcFillingRatio * scalekHz
          
    canvasRate2.cd()     
    canvasRate2.SetGridx()
    canvasRate2.SetGridy()   
    canvasRate2.SetLogy()
    paintedGraph.Draw("APZ")
    canvasRate2.Update();
    paintedGraph.GetXaxis().SetRangeUser(0, 100);
    #paintedGraph.GetYaxis().SetRangeUser(10 * scalekHz, 50000000 * scalekHz);
    paintedGraph.GetYaxis().SetRangeUser(1, 1000);
    canvasRate2.Update();     
          
    print ("createad rateCumul_withTEff "  + rateCumul_withTEff.GetName() )           
    print ("createad paintedGraph       "  + paintedGraph.GetName() )       
    return rateCumul_withTEff , paintedGraph


def makeRatePlots(algoDir, lineColor) :
    print (algoDir.GetName())
    algoDir.ls()
    
    c1 = TCanvas('canvas_' + algoDir.GetName(), algoDir.GetName().replace("_", " "), 200, 10, 700, 900)
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
            candPt_rateCumul = candPt.GetCumulative(False, "_");
            
#             trying to ge agreement woth Carlos plot, but it is not so easy, the last bin contains both pt 100 GeV and overflows, 
#             corr = candPt.GetBinContent(200)
#             print ("corr " + str(corr))
#             for iBin in range(0, candPt_rateCumul.GetNbinsX(), 1) : 
#                 candPt_rateCumul.AddBinContent(iBin, -corr)
            
            candPt_rateCumul.SetName(algoDir.GetName() + "_" + candPt_rateCumul.GetName().replace("candPt_", "rate") ) #+ "_" + version
            candPt_rateCumul.SetTitle(algoDir.GetName().replace("_", " ") + ", " + version)
            
            
            c1.cd(3).SetGridx()
            c1.cd(3).SetGridy()   
            c1.cd(3).SetLogy()
            
            rateCumul_withTEff, paintedGraph = makeRatePlotWithEff(candPt_rateCumul.Clone(candPt_rateCumul.GetName() + "_copy"), lineColor, c1.cd(3), c1.cd(4))
            rateCumul_withTEffs.append(rateCumul_withTEff)
            paintedGraphs.append(paintedGraph)
            
            candPt_rateCumul.Scale(0.001)
            candPt_rateCumul.GetYaxis().SetTitle("rate [kHz]")
            
            print("candPt: " + obj.GetName() + " candPt_rateCumul " + candPt_rateCumul.GetName() + " " + candPt_rateCumul.GetTitle() )
            candPt_rateCumul.SetBinContent(1, 0);
            #candPt_rateCumul.Sumw2(False);
            candPt_rateCumul.Scale(scale) #TODO maybe it should be before Sumw2
            candPt_rateCumul.Sumw2(False);
    
            c1.cd(2).SetGridx()
            c1.cd(2).SetGridy()   
            candPt_rateCumul.SetLineColor(lineColor)
            candPt_rateCumul.Draw("")
            
            rateCumuls.append(candPt_rateCumul)
            print ("created rate plot " + candPt_rateCumul.GetName() + " name: " + candPt_rateCumul.GetTitle() )
       
    canvases.append(c1) 
# makeRatePlots #######################################################################       

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
        
        
ratesOnThreshHist = TH1D("ratesOnThreshHist", inputResults.replace("_", " "), rateCumuls.__len__(), 0, rateCumuls.__len__()) 
relativeRatesOnThreshHist = TH1D("relativeRatesOnThreshHist", "relativeRatesOnThreshHist", rateCumuls.__len__(), 0, rateCumuls.__len__()) 
ratesOnThreshHist.GetYaxis().SetTitle("rate [kHz]")
relativeRatesOnThreshHist.GetYaxis().SetTitle("rate [kHz]")

referenceRate = 13.679002  #omtf q12, PU200_v2_t35

for iAlgo, canvas in enumerate(canvases ) :
    if iAlgo >= 1 :
        canvas.cd(2)
        rateCumuls[0].DrawCopy("same hist")
        canvas.cd(2).Modified()
        canvas.Update()
    
    ptCutGev = 21.5
        
    if rateCumuls[iAlgo].GetName().find("nn_omtf") >= 0:
        ptCutGev = 22
    elif version.find("t58") >= 0 or version.find("t65") >= 0 or version.find("t74") >= 0 or version.find("t78") >= 0 or version.find("t80") >= 0 or version.find("t82") >= 0 or version.find("t98") >= 0:
        ptCutGev = 18
    else :
        lineColor = 1
        ptCutGev = 20    
    
    ptCutBin = rateCumuls[iAlgo].GetXaxis().FindBin(ptCutGev)        
    rateOnThresh = rateCumuls[iAlgo].GetBinContent(ptCutBin)   
    ratesOnThreshHist.Fill( iAlgo, rateOnThresh)
    ratesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() + " ptCut " + str(ptCutGev) + "GeV" )   
    

    relativeRatesOnThresh = rateOnThresh / referenceRate

    print("%s rate %f realitve rate %f " % (ratesOnThreshHist.GetXaxis().GetBinLabel(iAlgo +1), rateOnThresh, relativeRatesOnThresh) )

    relativeRatesOnThreshHist.Fill( iAlgo, relativeRatesOnThresh)
    relativeRatesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() + " ptCut " + str(ptCutGev) + "GeV")    
        
    outFile.cd()
    rateCumuls[iAlgo].Write()
    rateCumul_withTEffs[iAlgo].Write()
    paintedGraphs[iAlgo].Write()
    canvas.Write()

canvasComapre = TCanvas('canvasComapre' , "compare_ " + inputResults, 200, 10, 1400, 700)    
canvasComapre.Divide(2, 1)

canvasComapre.cd(1)
   
canvasComapre.cd(1).SetLeftMargin(0.4)
canvasComapre.cd(1).SetGridx()
canvasComapre.cd(1).SetGridy()
ratesOnThreshHist.GetYaxis().SetRangeUser(6, 20)
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

outFile.cd()
ratesOnThreshHist.Write()
canvasComapre.Write()
#%jsroot on
#from ROOT import gROOT 
gROOT.GetListOfCanvases().Draw()

#outFile.Close()

raw_input("Press ENTER to exit")

#execfile('ratePlots.py')
