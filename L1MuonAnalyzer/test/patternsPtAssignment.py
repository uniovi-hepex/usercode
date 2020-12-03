from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
from ROOT import kBlack, kBlue, kRed, kGreen, kMagenta, kCyan

import sys
import os
from array import array

gStyle.SetOptStat(0)

def rebin(hist, bins) :
    newHist = TH1D( (hist.GetName() + "_rebined"), hist.GetTitle(), bins.__len__()-1, bins)
    newHist.Sumw2(False)    
    for iBin in range(0, hist.GetNbinsX() +1) : 
        newHist.Fill(hist.GetBinCenter(iBin), hist.GetBinContent(iBin) )
        #print hist.GetBinCenter(iBin), " ", hist.GetBinContent(iBin)
    
    newHist.Sumw2(False)  
    return  newHist


def makeEfficiency(passed, total, title, lineColor):
    if TEfficiency.CheckConsistency(passed, total) :
        #print "makeEfficiency passed ", passed.GetName() 
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


xBins = []
edge = 0
for i in range(0, 200, 1) :
    #print ("edge", edge)
    xBins.append(edge) 
    
    if edge < 30 :
        edge = edge + 1
    elif edge < 50 :
        edge = edge + 2    
    elif edge < 100 :
        edge = edge + 5
    elif edge < 140 :
        edge = edge + 10     
    elif edge < 180: 
        edge = edge + 20
    elif edge == 180:
        break
    
xBins.append(199)     
xBins.append(300)    
xBins.append(400)  
xBins.append(500) 
xBins.append(600)
xBins.append(700)
xBins.append(800)
xBins.append(900)    
xBins = array('d', xBins)

print('xBins', xBins)

#version = "PU200_v2_t" + sys.argv[1]
version = sys.argv[1]
#inputResults = 'MuFlatPt_' + version #+ "_test" 
inputResults = version #+ "_test" 

omtf_type = 2022 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v57_1efficiency.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v28_10Files.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t33/results/omtfAnalysis2.root' )

if "SingleMu_" in version :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/eff_SingleMu/omtfAnalysis2_eff_' + version + '.root' )
    omtf_type = 2022
elif version < "MuFlatPt_PU200_v3_t70" :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )
else :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/L1Trigger/L1TMuonOverlapPhase1/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )

if "0x0006" in version  or "t35" in version :
    omtf_type = 2018


print("omtf_type: ", omtf_type)

print (histFile)

if not os.path.exists(inputResults):
    os.mkdir(inputResults)
    
outFile = TFile(inputResults + "/efficiencyPlots.root", "RECREATE")

#histFile.ls()

efficiencyDir = histFile.Get("L1MuonAnalyzerOmtf/efficiency")
efficiencyDir.ls()
canvases = []
efficiencies = []
efficienciesHist1 = []
efficienciesHist2 = []

efficienciesOnThresh = []

efficienciesVsEta = []

thresholdHists = []

def makeEfficiencyPlots(ptCutGev, platCutGev, lineColor) :
    ptCutBin = ptGenVsPtCand.GetYaxis().FindBin(ptCutGev)
    accpetedVsPtGen = ptGenVsPtCand.ProjectionX(ptGenVsPtCand.GetName() + "_accpetedVsPtGen_" + str(ptCutGev) + "_GeV", ptCutBin, -1)
    allVsPtGen = ptGenVsPtCand.ProjectionX(ptGenVsPtCand.GetName() + "_allVsPtGen" + str(ptCutGev) + "_GeV", -1, -1)
    
    effName = canvasTitle.replace("_", " ") + " " + version
    eff = makeEfficiency(rebin(accpetedVsPtGen, xBins), rebin(allVsPtGen, xBins), effName + ";ptGen [GeV];efficiency", 4)
    eff.SetName(eff.GetName().replace("ptGenVsPtCand", "efficiency").replace("allVsPtGen", "ptCut_").replace("_rebined_clone", "")   ) 
    #print "eff.GetName()", eff.GetName()
    #eff.Draw("")
    #accpetedVsPtGen.Draw("same")
    efficiencies.append(eff)
    
    ##########################
    effHistName = accpetedVsPtGen.GetName().replace("_ptGenVsPtCand", "").replace("accpetedVsPtGen", "effOnPtCut")
    effHist1 = accpetedVsPtGen.Clone(effHistName + "_1")
    effHist1.Divide(allVsPtGen)
    effHist1.SetLineColor(lineColor)
    effHist1.Draw("hist")
    effHist1.GetYaxis().SetTitle("efficiency")
    effHist1.GetYaxis().SetRangeUser(0, 1.05)
    effHist1.SetTitle(effName)
    efficienciesHist1.append(effHist1)
    #print ("effHist1 " + effHist1.GetName())
     
    #################### calulating efficiency on the plataou
    platCutBin = allVsPtGen.GetXaxis().FindBin(platCutGev)
    allIntegrated = allVsPtGen.Integral(platCutBin, -1);
    
    accpetedIntegrated = accpetedVsPtGen.Integral(platCutBin, -1)
    
    effONPtCut = 0
    if allIntegrated > 0 :
        effONPtCut = accpetedIntegrated / allIntegrated
   
   
    thresholdHist = effHist1.Clone(effHist1.GetName() + "_threshold")
    for iBin in range(0, thresholdHist.GetNbinsX() +1, +1) :
        thresholdHist.SetBinContent(iBin, 1)
        
    thresholdHist.SetBinContent(0, 0)    
    for iBin in range(2, effHist1.GetNbinsX(), +1) :
        effOnSelectedPtCut = (effHist1.GetBinContent(iBin) + effHist1.GetBinContent(iBin-1) ) / 2.
        selectedPtCut = effHist1.GetBinLowEdge(iBin)
        #print("    iBin %i selectedPtCut %i effOnSelectedPtCut %f" % (iBin, selectedPtCut, effOnSelectedPtCut))
        if effOnSelectedPtCut > (0.85 *  effONPtCut) :
            break
        
        thresholdHist.SetBinContent(iBin, 0)
    
    #print("%s - %.3f" % (ptGenVsPtCand.GetName(), effONPtCut ) ) 
    #print("ptCutBin %f ptCutGev %f platCutBin %f platCutGev %f  accpetedIntegrated %i allIntegrated %i effONPtCut - %.3f" % (ptCutBin, ptCutGev, platCutBin, platCutGev, accpetedIntegrated, allIntegrated, effONPtCut ) )
    print("ptCutGev %f platCutGev %f  accpetedIntegrated %i allIntegrated %i effONPtCut - %.3f - effOnSelectedPtCut %.3f - selectedPtCut %f " % (ptCutGev, platCutGev, accpetedIntegrated, allIntegrated, effONPtCut, effOnSelectedPtCut, selectedPtCut ) )
    
    if selectedPtCut > 40 :
        effHist1.GetXaxis().SetRangeUser(0, selectedPtCut + 100)
    else :
        effHist1.GetXaxis().SetRangeUser(0, selectedPtCut + 20)
    
    thresholdHist.SetLineColor(kRed)
    thresholdHist.Draw("hist same")
    thresholdHists.append(thresholdHist)
    efficienciesOnThresh.append(effONPtCut)


orgPtBins = [   2  ,
                3.5,
                4  ,
                4.5,
                5  ,
                6  ,
                7  ,
                8  ,
                10 ,
                12 ,
                14 ,
                16 ,
                18 ,
                20 ,
                22 ,
                24 ,
                26 ,
                28 ,
                30 ,
                35 ,
                40 ,
                45 ,
                50 ,
                60 ,
                70 ,
                80 ,
                100,
                200 ]

#print (ptGenVsPtCand)
#ptGenVsPtCand.GetName()
canvasTitle = version
canvasTitle = canvasTitle.replace("_", " ")
c1 = TCanvas('canvas_1_' + canvasTitle, canvasTitle, 200, 10, 900, 900)
c1.Divide(4, 4)

c2 = TCanvas('canvas_12_' + canvasTitle, canvasTitle, 200, 10, 900, 900)
c2.Divide(4, 4)

#print ('created canvas ' + c1.GetName())

canvases.append(c1)
        
for iAlgo, obj in enumerate(efficiencyDir.GetListOfKeys() ) :
    ptGenVsPtCand = obj.ReadObj()
    if isinstance(ptGenVsPtCand, TH2D) and "q12" in ptGenVsPtCand.GetName(): 
        iPad = 1
        for orgPtBin in orgPtBins :
            if iPad <= 16 :
                c1.cd(iPad).SetGridx()
                c1.cd(iPad).SetGridy()
            else :
                c2.cd(iPad-16).SetGridx()
                c2.cd(iPad-16).SetGridy()        
            #makeEfficiencyPlots(5)
            lineColor = kBlue
            ptCut = orgPtBin
    
            if ptCut > 10 :
                offset = 10
            elif ptCut > 20 :
                offset = 20
            elif ptCut > 40 :
                offset = 30
            else :
                offset = 7
            makeEfficiencyPlots(ptCut, ptCut + offset, lineColor)
            
            iPad += 1

c1.Modified()
c1.Update()    
    
c2.Modified()
c2.Update()    

#%jsroot on
#from ROOT import gROOT 
gROOT.GetListOfCanvases().Draw()

raw_input("Press ENTER to exit")
#execfile('efficiencyPlots.py')
