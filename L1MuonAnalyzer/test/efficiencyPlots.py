from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
import sys
import os

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
    

#version = "PU200_v2_t" + sys.argv[1]
version = sys.argv[1]
inputResults = 'MuFlatPt_' + version #+ "_test" 

#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1_10Files_withMatching.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v21_1.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis2_v45efficiency.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/omtf/omtfAnalysis_newerSAmple_v28_10Files.root' )
#histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t33/results/omtfAnalysis2.root' )
histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )


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

def makeEfficiencyPlots(ptCutGev, platCutGev, lineColor) :
    #print (ptGenVsPtCand)
    #ptGenVsPtCand.GetName()
    canvasTitle = ptGenVsPtCand.GetName()[ : ptGenVsPtCand.GetName().find("ptGenVsPtCand")] + "ptCut " + str(ptCutGev) + "GeV" #+ "_" + version
    canvasTitle = canvasTitle.replace("_", " ")
    c1 = TCanvas('canvas_' + ptGenVsPtCand.GetName() + "_" + str(ptCutGev), canvasTitle, 200, 10, 700, 900)
    c1.Divide(2, 2)
    c1.cd(1).SetGridx()
    c1.cd(1).SetGridy()
    #print ('created canvas ' + c1.GetName())
    ptGenVsPtCand.Draw('colz')
    # c1.Draw()
    ##########################
    c1.cd(2).SetGridx()
    c1.cd(2).SetGridy()
    
    ptCutBin = ptGenVsPtCand.GetYaxis().FindBin(ptCutGev)
    accpetedVsPtGen = ptGenVsPtCand.ProjectionX(ptGenVsPtCand.GetName() + "_accpetedVsPtGen_" + str(ptCutGev) + "_GeV", ptCutBin, -1)
    allVsPtGen = ptGenVsPtCand.ProjectionX(ptGenVsPtCand.GetName() + "_allVsPtGen" + str(ptCutGev) + "_GeV", -1, -1)
    
    effName = canvasTitle.replace("_", " ") + " " + version
    eff = makeEfficiency(accpetedVsPtGen, allVsPtGen, effName + ";ptGen [GeV];efficiency", 4)
    eff.Draw("")
    efficiencies.append(eff)
    
    ##########################
    c1.cd(3).SetGridx()
    c1.cd(3).SetGridy()
    effHistName = accpetedVsPtGen.GetName().replace("_ptGenVsPtCand", "").replace("accpetedVsPtGen", "effOnPtCut")
    effHist1 = accpetedVsPtGen.Clone(effHistName + "_1")
    effHist1.Divide(allVsPtGen)
    effHist1.SetLineColor(lineColor)
    effHist1.Draw("hist")
    effHist1.GetYaxis().SetTitle("efficiency")
    effHist1.GetYaxis().SetRangeUser(0, 1.05)
    effHist1.SetTitle(effName)
    efficienciesHist1.append(effHist1)
    print ("effHist1 " + effHist1.GetName())

    ##########################
    rebin = 2
    c1.cd(4).SetGridx()
    c1.cd(4).SetGridy()
    effHist2 = accpetedVsPtGen.Rebin(rebin, effHistName + "_2")
    allVsPtGen2 = allVsPtGen.Rebin(rebin, allVsPtGen.GetName() + "_2")
    effHist2.Divide(allVsPtGen2)
    effHist2.SetLineColor(lineColor)
    effHist2.GetYaxis().SetRangeUser(0, 1.05)
    effHist2.GetYaxis().SetTitle("efficiency")
    effHist2.Draw("hist")
    effHist2.SetTitle(effName)
    efficienciesHist2.append(effHist2)        
    canvases.append(c1)
    
    #################### calulating efficiency on the plataou
    platCutBin = allVsPtGen.GetXaxis().FindBin(platCutGev)
    allIntegrated = allVsPtGen.Integral(platCutBin, -1);
    
    accpetedIntegrated = accpetedVsPtGen.Integral(platCutBin, -1);
    #print (ptGenVsPtCand.GetName() +  " " + str(accpetedIntegrated / allIntegrated) ) 
    print("%s - %.3f" % (ptGenVsPtCand.GetName(), (accpetedIntegrated / allIntegrated) ) ) # 4.000000
    efficienciesOnThresh.append(accpetedIntegrated / allIntegrated)

 
        
for iAlgo, obj in enumerate(efficiencyDir.GetListOfKeys() ) :
    ptGenVsPtCand = obj.ReadObj()
    if isinstance(ptGenVsPtCand, TH2D):
        #makeEfficiencyPlots(5)
        lineColor = 2
        ptCut = 21.5
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 21.5 #+3
        else :
            lineColor = 1
            ptCut = 20 #+5
            
        makeEfficiencyPlots(ptCut, 25, lineColor)
        #makeEfficiencyPlots(0, lineColor)
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 21.5 + 3
        else :
            lineColor = 1
            ptCut = 20 + 5
            
        makeEfficiencyPlots(ptCut, 30, lineColor)
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 21.5 + 20
        else :
            lineColor = 1
            ptCut = 20 + 20
            
        makeEfficiencyPlots(ptCut, 45, lineColor)

    
    
efficienciesOnThreshHist = TH1D("efficienciesOnThresh", inputResults, efficienciesOnThresh.__len__(), 0, efficienciesOnThresh.__len__())   
        
for iAlgo, canvas in enumerate(canvases ) :
    efficienciesOnThreshHist.Fill( iAlgo, efficienciesOnThresh[iAlgo])
    efficienciesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() )
    
        
    outFile.cd()
    efficiencies[iAlgo].Write()
    efficienciesHist1[iAlgo].Write()
    efficienciesHist2[iAlgo].Write()
    canvas.Write()
    
    if iAlgo >= 1 :
        canvas.cd(3)
        efficienciesHist1[0].DrawCopy("same hist")
        canvas.cd(3).Modified()
        canvas.Update()
        
        canvas.cd(4)
        efficienciesHist2[0].DrawCopy("same hist")
        canvas.cd(4).Modified()
        canvas.Update()
            
canvasComapre = TCanvas('canvasComapre' , "compare " + inputResults, 200, 10, 700, 700)    
canvasComapre.cd()
   
canvasComapre.SetLeftMargin(0.42)
canvasComapre.SetGridx()
canvasComapre.SetGridy()
efficienciesOnThreshHist.GetXaxis().SetLabelSize(0.02)
efficienciesOnThreshHist.GetYaxis().SetRangeUser(0.92, 0.97)
efficienciesOnThreshHist.GetYaxis().SetLabelSize(0.02)
efficienciesOnThreshHist.GetYaxis().SetTitle("efficiency")
efficienciesOnThreshHist.SetFillColor(0)
efficienciesOnThreshHist.SetFillStyle(3001)
efficienciesOnThreshHist.Draw("HBAR")

outFile.cd()
efficienciesOnThreshHist.Write()
canvasComapre.Write()

#pad1 = TPad( 'pad1', 'The pad with the function',  0.03, 0.62, 0.50, 0.92, 21 )
#pad2 = TPad( 'pad2', 'The pad with the histogram', 0.51, 0.62, 0.98, 0.92, 21 )
#pad3 = TPad( 'pad3', 'The pad with the histogram', 0.03, 0.02, 0.97, 0.57, 21 )
#pad1.Draw()
#pad2.Draw()
#pad3.Draw()

#%jsroot on
#from ROOT import gROOT 
gROOT.GetListOfCanvases().Draw()

raw_input("Press ENTER to exit")
#execfile('efficiencyPlots.py')
