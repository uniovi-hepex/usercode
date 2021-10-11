from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
import sys
import os
from array import array

gStyle.SetOptStat(0)

def rebin(hist, bins) :
    newHist = TH1D( (hist.GetName() + "_rebined"), hist.GetTitle(), bins.__len__()-1, bins)
    newHist.Sumw2(False)    
    for iBin in range(0, hist.GetNbinsX() +1) : 
        newHist.Fill(hist.GetBinCenter(iBin), hist.GetBinContent(iBin) )
        #print(hist.GetBinCenter(iBin), " ", hist.GetBinContent(iBin) )
    
    newHist.Sumw2(False)  
    return  newHist


def makeEfficiency(passed, total, title, lineColor):
    if TEfficiency.CheckConsistency(passed, total) :
        print("makeEfficiency passed ", passed.GetName()) 
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

elif "t41" in version or "t66" in version :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )
else :
    histFile = TFile( '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/L1Trigger/L1TMuonOverlapPhase1/test/crab/crab_omtf_nn_MC_analysis_' + inputResults + '/results/omtfAnalysis2.root' )


#histFile = TFile( '/home/kbunkow/CMSSW/CMSSW_12_1_0_pre3/src/L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/omtfAnalysis2_eff_SingleMu_tt10_displ_test.root' ) # displaced
histFile = TFile( '/home/kbunkow/CMSSW/CMSSW_12_1_0_pre3/src/L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/omtfAnalysis2_eff_SingleMu_tt10_test_50files.root' )


#if "_t1" in version :
#    omtf_type = 2018
        
if "0x0006" in version  or "t35" in version or "t85" in version :
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

def makeEfficiencyPlots(ptCutGev, platCutGev, lineColor) :
    #print (ptGenVsPtCand)
    #ptGenVsPtCand.GetName()
    canvasTitle = ptGenVsPtCand.GetName()[ : ptGenVsPtCand.GetName().find("ptGenVsPtCand")] + "ptCut " + str(ptCutGev) + "GeV" #+ "_" + version
    canvasTitle = canvasTitle.replace("_", " ")
    c1 = TCanvas('canvas_' + ptGenVsPtCand.GetName() + "_" + str(ptCutGev), canvasTitle, 200, 10, 900, 900)
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
    #eff = makeEfficiency(accpetedVsPtGen, allVsPtGen, effName + ";ptGen [GeV];efficiency", 4)
    eff = makeEfficiency(rebin(accpetedVsPtGen, xBins), rebin(allVsPtGen, xBins), effName + ";ptGen [GeV];efficiency", 4)
    #eff = rebin(accpetedVsPtGen, xBins)
    eff.SetName(eff.GetName().replace("ptGenVsPtCand", "efficiency").replace("allVsPtGen", "ptCut_").replace("_rebined_clone", "")   ) 
    #omtf_q12_ptGenVsPtCand_eta_0.82_1.24_qualityCut_12_allVsPtGen18_GeV_rebined_clone
    print("eff.GetName()", eff.GetName() )
    eff.Draw("")
    #accpetedVsPtGen.Draw("same")
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
    effHist1.GetXaxis().SetRangeUser(0, 100)
    effHist1.GetYaxis().SetRangeUser(0, 1.05)
    effHist1.SetTitle(effName)
    efficienciesHist1.append(effHist1)
    print ("effHist1 " + effHist1.GetName())

    ##########################
    rebinFactor = 2
    c1.cd(4).SetGridx()
    c1.cd(4).SetGridy()

    
    effHist2 = accpetedVsPtGen.Rebin(rebinFactor, effHistName + "_2")
    allVsPtGen2 = allVsPtGen.Rebin(rebinFactor, allVsPtGen.GetName() + "_2")
    effHist2.Divide(allVsPtGen2)
    effHist2.SetLineColor(lineColor)
    effHist2.GetXaxis().SetRangeUser(2, 100)
    effHist2.GetYaxis().SetRangeUser(0, 1.05)
    effHist2.GetYaxis().SetTitle("efficiency")
    effHist2.SetTitle(effName)
    efficienciesHist2.append(effHist2)   
    #c1.cd(4).SetLogx()
    #effHist2.Draw("hist")

    # efficiency vs eta inseed of the above turn on curve
    #omtf_q12_ptGenVsPtCand_eta_0.82_1.24_qualityCut_12
    #omtf_q12_allCandsEta__qualityCut_12_ptGenCut_10
#omtf_q12_allCandsEta__qualityCut_12_ptGenCut_25
#omtf_q12_aceptedCandsEta__qualityCut_12_ptGenCut_25_ptL1Cut_20
    
    ####################eff vs eta
    allCandsEtaName = ptGenVsPtCand.GetName().replace("ptGenVsPtCand_eta_0.82_1.24", "allCandsEta_") + "_ptGenCut_25"
    #print ("allCandsEtaName " + allCandsEtaName)
    allCandsEta = efficiencyDir.Get(allCandsEtaName)
    
    aceptedCandsEtaName = ptGenVsPtCand.GetName().replace("ptGenVsPtCand_eta_0.82_1.24", "aceptedCandsEta_") + "_ptGenCut_25_ptL1Cut_20"
        
    aceptedCandsEta = efficiencyDir.Get(aceptedCandsEtaName)
    if not aceptedCandsEta :
        aceptedCandsEtaName = ptGenVsPtCand.GetName().replace("ptGenVsPtCand_eta_0.82_1.24", "aceptedCandsEta_") + "_ptGenCut_25_ptL1Cut_18"
        aceptedCandsEta = efficiencyDir.Get(aceptedCandsEtaName)
    
    if aceptedCandsEta :
        effVsEta = makeEfficiency(aceptedCandsEta, allCandsEta, aceptedCandsEta.GetTitle().replace("allCands", "efficiency"), lineColor)
        effVsEta.SetName(effVsEta.GetName().replace("_clone", "").replace("allCandsEta", "efficiencyVsEta"))
        effVsEta.Draw("")
        c1.cd(4).Update()
        effVsEta.GetPaintedGraph().GetYaxis().SetRangeUser(0.8, 1.05)
        efficienciesVsEta.append(effVsEta)
    
    canvases.append(c1)
    
    #################### calulating efficiency on the plataou
    platCutBin = allVsPtGen.GetXaxis().FindBin(platCutGev)
    allIntegrated = allVsPtGen.Integral(platCutBin, -1) +1; # TODD check why sometimes is 0
    
    accpetedIntegrated = accpetedVsPtGen.Integral(platCutBin, -1);
    #print (ptGenVsPtCand.GetName() +  " " + str(accpetedIntegrated / allIntegrated) ) 
    print("%s - %.3f" % (ptGenVsPtCand.GetName(), (accpetedIntegrated / allIntegrated) ) ) # 4.000000
    efficienciesOnThresh.append(accpetedIntegrated / allIntegrated)

     
        
for iAlgo, obj in enumerate(efficiencyDir.GetListOfKeys() ) :
    ptGenVsPtCand = obj.ReadObj()
    if isinstance(ptGenVsPtCand, TH2D) and "ptGenVsPtCand" in ptGenVsPtCand.GetName() :
        #makeEfficiencyPlots(5)
        lineColor = 2
        ptCut = 22
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 22 #+3
        elif omtf_type ==  2018:
            lineColor = 1
            ptCut = 20 #20 #+5
        elif omtf_type ==  2022 :   
            lineColor = 1
            ptCut = 20 #20 #+5
            
        makeEfficiencyPlots(ptCut, 25, lineColor)
        #makeEfficiencyPlots(0, lineColor)
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 22 + 2
        elif omtf_type ==  2018:
            lineColor = 1
            ptCut = 22 #20 #+5    
        elif omtf_type ==  2022 :
            lineColor = 1
            ptCut = 22 #22
            
        makeEfficiencyPlots(ptCut, 30, lineColor)
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 26 + 2
        elif omtf_type ==  2018:
            lineColor = 1
            ptCut = 26 #20 #+5    
        elif omtf_type ==  2022 :
            lineColor = 1
            ptCut = 26 # 26
        
        makeEfficiencyPlots(ptCut, 34, lineColor)
            
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 22 + 20
        else :
            lineColor = 1
            ptCut = 40
            
        makeEfficiencyPlots(ptCut, 45, lineColor)
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 10
        elif omtf_type ==  2018:
            lineColor = 1
            ptCut = 10 #20 #+5    
        elif omtf_type ==  2022 :
            lineColor = 1
            ptCut = 10 # 26
            
        makeEfficiencyPlots(ptCut, 12, lineColor)
        
        
        if ptGenVsPtCand.GetName().find("nn_omtf") >= 0 :
            ptCut = 5
        elif omtf_type ==  2018:
            lineColor = 1
            ptCut = 5 #20 #+5    
        elif omtf_type ==  2022 :
            lineColor = 1
            ptCut = 5# 26
            
        makeEfficiencyPlots(ptCut, 8, lineColor)

        ptGenVsPtCand.GetXaxis().SetRangeUser(0, 100)
        ptGenVsPtCand.GetYaxis().SetRangeUser(0, 100)

    
    
efficienciesOnThreshHist = TH1D("efficienciesOnThresh", inputResults, efficienciesOnThresh.__len__(), 0, efficienciesOnThresh.__len__())   
        
for iAlgo, canvas in enumerate(canvases ) :
    efficienciesOnThreshHist.Fill( iAlgo, efficienciesOnThresh[iAlgo])
    efficienciesOnThreshHist.GetXaxis().SetBinLabel(iAlgo +1, canvases[iAlgo].GetTitle() )
    
        
    outFile.cd()
    efficiencies[iAlgo].Write()
    efficienciesHist1[iAlgo].Write()
    efficienciesHist2[iAlgo].Write()
    #if outFile.FindObject(efficienciesVsEta[iAlgo].GetName() ) == None :
    if efficienciesVsEta.__len__() > 0:
        if (iAlgo == 0) or (efficienciesVsEta[iAlgo].GetName() != efficienciesVsEta[iAlgo -1].GetName() ) :
            efficienciesVsEta[iAlgo].Write()
        
    canvas.Write()
    
    if iAlgo >= 1 :
        canvas.cd(3)
        efficienciesHist1[0].DrawCopy("same hist")
        canvas.cd(3).Modified()
        canvas.Update()
        
        canvas.cd(4)
        #efficienciesHist2[0].DrawCopy("same hist")
        canvas.cd(4).Modified()
        canvas.Update()
            
canvasComapre = TCanvas('canvasComapre' , "compare " + inputResults, 200, 10, 700, 700)    
canvasComapre.cd()
   
canvasComapre.SetLeftMargin(0.42)
canvasComapre.SetGridx()
canvasComapre.SetGridy()
efficienciesOnThreshHist.GetXaxis().SetLabelSize(0.02)
efficienciesOnThreshHist.GetYaxis().SetRangeUser(0.92, 0.98)
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

input("Press ENTER to exit")
#execfile('efficiencyPlots.py')
