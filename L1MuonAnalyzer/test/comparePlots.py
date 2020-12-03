from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TLegend
from ROOT import kBlack, kBlue, kRed, kGreen, kMagenta, kCyan
from array import array

from libPyROOT import TDirectory
import os
import sys
from __builtin__ import True
    


def makeUniqueFileName(path, name):
    fn = os.path.join(path, name)
    if not os.path.exists(fn):
        return fn

    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s%d%s' % (name, i, ext))

    for i in xrange(2, sys.maxint):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

    return None
###########################################


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0);


# leg -> SetHeader("here is a beautiful header")

plotsDir = '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/usercode/L1MuonAnalyzer/test/'

first = True

logScalePads = []
logScalePadNum = 0
effHistCopys = []

def drawEff(effFile, type, quality, ptCut, lineColor, legend, pTresh = "0.5") :
    global first
    doEff = True
    
    if not doEff :
        if type == "nn_omtf" :
            histNametype + "_q" + quality + "_pTresh_" + pTresh + "_eta_0.82_1.24_qualityCut_" + quality + "_effOnPtCut_" + ptCut + "_GeV_1"
        else :
            histName = type + "_q" + quality + "_eta_0.82_1.24_qualityCut_" + quality + "_effOnPtCut_" + ptCut + "_GeV_1"
            #omtf_q12_eta_0.82_1.24_qualityCut_12_effOnPtCut_20_GeV_1
    else  :    
        if type == "nn_omtf" :
            histNametype + "_q" + quality + "_pTresh_" + pTresh + "_eta_0.82_1.24_qualityCut_" + quality + "_effOnPtCut_" + ptCut + "_GeV_1"
        else :
            histName = type + "_q" + quality + "_efficiency_eta_0.82_1.24_qualityCut_" + quality + "_ptCut_" + ptCut + "_GeV"    
            #omtf_q12_efficiency_eta_0.82_1.24_qualityCut_12_ptCut_18_GeV
        
    print (effFile)    
    print (histName) 

    effHist = effFile.Get(histName)
    if effHist is None :
        print ("no histogram found: ", histName)
           
    
    effHist.SetLineColor(lineColor)
    print ("first " + str(first) )
    if first :
        if not doEff :
            effHist.Draw("hist")
            effHist.GetXaxis().SetRangeUser(0, 200)
            effHist.GetYaxis().SetRangeUser(0, 1.05)
        else :    
            effHist.Draw("AE")

        print ("line 56")
    else:
        if not doEff :
            effHist.Draw("hist same")
        else :
            effHist.Draw("P same")       
        print ("line 63")
        
    legend.AddEntry(effHist)  # , "OMTF", "lep");
    
    
#     if first :
#         pad = TPad('pad_' + str(logScalePads.__len__()), 'pad', 0.4,  0.27,  0.99,  0.75)
#         
#         pad.Draw()
#         pad.cd()
#         #pad.SetLogy()
#         pad.SetGridx()
#         pad.SetGridy()
#         pad.SetRightMargin(0.01)
#         pad.SetTopMargin(0.01)
#         pad.SetLeftMargin(0.1)
#         
#         #pad.SetLogx()
#         #pad.SetLogy()
#         if not doEff :
#         effHistCopy = effHist.DrawCopy("hist")
#         effHistCopy.GetXaxis().SetRangeUser(2, 20)
#         effHistCopy.GetYaxis().SetRangeUser(0.00, 0.1)
#         effHistCopy.GetYaxis().SetTitleOffset(1.5)
#         effHistCopys.append(effHistCopy)
#         logScalePads.append(pad)
#     else:
#         logScalePads[logScalePadNum].cd()
#         print("pad name " + logScalePads[logScalePadNum].GetName() )
#         effHistCopy = effHist.DrawCopy("hist same")    
#         effHistCopys.append(effHistCopy)
#         print ("line 84")
         
###################################################

firstEtaPlot = True
def drawEffVsEta(effFile, type, quality, lineColor) :
    global firstEtaPlot 
    
    #mtf_q12_efficiencyVsEta__qualityCut_12_ptGenCut_25
    
    if type == "nn_omtf" :
        effHistName = type + "_q" + quality + "_pTresh_" + pTresh + "_efficiencyVsEta__qualityCut_" + quality + "_ptGenCut_25"
    else :
        effHistName = type + "_q" + quality + "_efficiencyVsEta__qualityCut_" + quality + "_ptGenCut_25"
        #omtf_q12_eta_0.82_1.24_qualityCut_12_effOnPtCut_20_GeV_1
       
    effHist = effFile.Get(effHistName)
    print ("line 109 ", effHistName, effHist)   
    if effHist :    
        effHist.SetLineColor(lineColor)
        print ("first " + str(first) )
        if firstEtaPlot :
            #effHist.GetPaintedGraph().GetYaxis().SetRangeUser(0.8, 1.05)
            effHist.Draw("")
            firstEtaPlot =  False
        else:
            effHist.Draw("same")   


rateFiles = []
fillPat = 3002
def drawRate(rateFileDir, type, quality, lineColor, pTresh = "0.5") :
    global first
    global fillPat
    
    rateFile = TFile(plotsDir + rateFileDir + 'ratePlots.root' )
    rateFiles.append(rateFile)
    
    withTEff = "" # _copy_allEventsHist_clone_copy"
    if type == "nn_omtf" :
        effHist = rateFile.Get(type + "_q" + quality + "_pTresh_" + pTresh + "_rate_qualityCut_" + quality + "_" + withTEff)
    else :
        effHist = rateFile.Get(type + "_q" + quality + "_rate_qualityCut_" + quality + "_" + withTEff)
        
    print (rateFile) 
      
    
    effHist.SetLineColor(lineColor)
    #effHist.SetFillColor(lineColor);
    #effHist.SetFillStyle(fillPat)
    fillPat += 1
    #effHist.SetFillColorAlpha(lineColor, 0.5)
   
    global first
    if first :
        effHist.GetXaxis().SetRangeUser(0, 70)
        effHist.GetYaxis().SetRangeUser(1, 100)
        if withTEff == "" :
            effHist.Draw("hist")
        else :
            effHist.Draw("APZ")
        
    else:
        if withTEff == "" :
            effHist.Draw("hist same")   
        else :
            effHist.Draw("PZ")
         
    legend.AddEntry(effHist)  # , "OMTF", "lep");
    first = False
##################################################

######################################
c1 = TCanvas('canvas_efficiency_1', 'canvas_efficiency', 200, 10, 950, 500)
c1.Divide(2, 1)
c1.cd(1)
c1.cd(1).SetGridx()
c1.cd(1).SetGridy()
c1.cd(1).cd()

c1.cd(2).SetGridx()
c1.cd(2).SetGridy()

c2 = TCanvas('canvas_efficiency_2', 'canvas_efficiency_2', 200, 510, 950, 500)
c2.Divide(2, 1)
c2.cd(1).SetGridx()
c2.cd(1).SetGridy()
c2.cd(1).cd()

c2.cd(2).SetGridx()
c2.cd(2).SetGridy()

c3 = TCanvas('canvas_efficiency_3', 'canvas_efficiency_3', 200, 510, 950, 500)
c3.Divide(1, 1)
c3.cd(1).SetGridx()
c3.cd(1).SetGridy()
c3.cd(1).cd()

c3.cd(2).SetGridx()
c3.cd(2).SetGridy()

#legendEff1 = TLegend(0.06, 0.8, 0.57, 0.997)
legendEff1 = TLegend(0.21, 0.12, 0.75, 0.25)

#legend.SetHeader(header.c_str())
# leg -> SetBorderSize(0);
legendEff1.SetFillStyle(0)
legendEff1.SetBorderSize(0)
legendEff1.SetTextSize(0.03)
legendEff1.SetMargin(0.2)

legendEff2 = legendEff1.Clone()
legendEff3 = legendEff1.Clone()

eff_c1 = c1.cd(1)
eff_c2 = c2.cd(1)
eff_c3 = c2.cd(2)

# eff_c1.SetTopMargin(0.2)
# eff_c2.SetTopMargin(0.2)
# eff_c3.SetTopMargin(0.2)

effFiles = []

def drawEffs(fileDir, type, quality, lineColor, pTresh = "0.5" ) :
    global first
    global logScalePadNum
    print ("first " + str(first) )
    effFile = TFile(plotsDir + fileDir + 'efficiencyPlots.root' )
    effFiles.append(effFile)
    if type == "omtf" :
        eff_c1.cd()
        print (eff_c1.GetName() )
        logScalePadNum = 0
        drawEff(effFile, type, quality, "20", lineColor, legendEff1)
        
        logScalePadNum = 1 
        eff_c2.cd()
        drawEff(effFile, type, quality, "22", lineColor, legendEff2)

        logScalePadNum = 2
        eff_c3.cd()
        drawEff(effFile, type, "12", "26", lineColor, legendEff3)
        
    if type == "nn_omtf" :
        logScalePadNum = 0
        print (eff_c1.GetName() )
        eff_c1.cd()
        drawEff(effFile, type, quality, "22", lineColor, legendEff1, pTresh)
        logScalePadNum = 1
        eff_c2.cd()
        drawEff(effFile, type, quality, "24", lineColor, legendEff2, pTresh)
        logScalePadNum = 2
        eff_c3.cd()
        drawEff(effFile, type, quality, "42", lineColor, legendEff3, pTresh)
    
    if type == "omtf_patsKB" :
        logScalePadNum = 0
        eff_c1.cd()
        drawEff(effFile, "omtf", "12", "18", lineColor, legendEff1)
        logScalePadNum = 1
        eff_c2.cd()
        drawEff(effFile, "omtf", "12", "20", lineColor, legendEff2)
        logScalePadNum = 2
        eff_c3.cd()
        drawEff(effFile, "omtf", "12", "24", lineColor, legendEff3)
        
    c3.cd(1)
    
    if type == "omtf_patsKB" :
        drawEffVsEta(effFile, "omtf", quality, lineColor)
    else :    
        drawEffVsEta(effFile, type, quality, lineColor)
        
    first = False

#OMTF 2018+
#c1.cd(1)
#drawEffs('MuFlatPt_PU200_v2_t44/', "omtf", "12", kBlack)
#drawEffs('MuFlatPt_PU200_v2_t35/', "omtf", "12", kBlack)


#drawEffs('MuFlatPt_PU200_v2_t51/', "omtf", "12", kGreen+1)


#drawEffs('MuFlatPt_PU200_v2_t55/', "omtf_patsKB", "12", kRed)
#drawEffs('MuFlatPt_PU200_v2_t56/', "omtf_patsKB", "12", kBlue)

#drawEffs('MuFlatPt_PU200_v2_t65/', "omtf_patsKB", "12", kGreen+1) #omtf_patsKB

#drawEffs('MuFlatPt_PU200_v3_t78/', "omtf_patsKB", "12", kBlue) #omtf_patsKB
#drawEffs('MuFlatPt_PU200_v3_t79/', "omtf_patsKB", "12", kRed)
#drawEffs('MuFlatPt_PU200_v3_t80/', "omtf_patsKB", "12", kRed)
#drawEffs('MuFlatPt_PU200_v3_t82/', "omtf_patsKB", "12", kRed)


#drawEffs('MuFlatPt_PU200_v2_t46/', "omtf", "12", kGreen+1)

#OMTF 2018
#drawEffs('MuFlatPt_PU200_v2_t52/', "omtf", "12", kGreen+1)
#effFile = TFile(plotsDir + 'MuFlatPt_PU200_v2_t56/' + 'efficiencyPlots.root' )


#drawEffs('MuFlatPt_PU200_v2_t41/', "nn_omtf", "12", kBlue, "0.4")
#drawEffs('MuFlatPt_PU200_v2_t66/', "nn_omtf", "12", kRed, "0.4")

#drawEffs('MuFlatPt_PU300_v2_t41/', "nn_omtf", "12", kGreen+1, "0.5")


#drawEff(effFile_t42_pu200, "nn_omtf", "12", "21.5", kRed)

#drawEff(eff_t43_pu300, "nn_omtf", "12", "21.5", kRed)

#drawEffs('MuFlatPt_PU200_v3_t70/', "omtf_patsKB", "12", kRed)
#drawEffs('MuFlatPt_PU200_v3_t71/', "omtf_patsKB", "12", kGreen)
#drawEffs('MuFlatPt_PU200_v3_t73/', "omtf", "12", kRed)

#drawEffs('SingleMu_t85_10f/', "omtf", "12", kGreen) # Arturs' pattern GPs_parametrised_v1_classProb3.xml

#drawEffs('SingleMu_t80/', "omtf_patsKB", "12", kCyan)

#drawEffs('SingleMu_0x0006_t79/', "omtf", "12", kBlack)
 #drawEffs('SingleMu_t74/', "omtf_patsKB", "12", kGreen)
#drawEffs('SingleMu_t76/', "omtf_patsKB", "12", kRed)
#drawEffs('SingleMu_t77/', "omtf_patsKB", "12", kBlue)

#drawEffs('SingleMu_t78/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t78_1/', "omtf_patsKB", "12", kRed)

#drawEffs('SingleMu_t81_-16/', "omtf_patsKB", "12", kGreen)
#drawEffs('SingleMu_t81_-8/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t82/', "omtf_patsKB", "12", kCyan)
#drawEffs('SingleMu_t82_05GeVBinning/', "omtf_patsKB", "12", kCyan)
#drawEffs('SingleMu_t83_test/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t83/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t90/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t91/', "omtf_patsKB", "12", kMagenta)

#drawEffs('SingleMu_t92/', "omtf_patsKB", "12", kGreen) 
#drawEffs('SingleMu_t93/', "omtf_patsKB", "12", kBlue)
#drawEffs('SingleMu_t94/', "omtf_patsKB", "12", kMagenta) #this one is incomplete
#drawEffs('SingleMu_t97/', "omtf_patsKB", "12", kGreen)
#drawEffs('SingleMu_t98/', "omtf_patsKB", "12", kGreen)

#drawEffs('SingleMu_t80/', "omtf_patsKB", "12", kRed)
#drawEffs('SingleMu_t80_test/', "omtf_patsKB", "12", kBlue) #finalize8 !!!!!!!!!!!!!!!!!!
#drawEffs('SingleMu_t103/', "omtf_patsKB", "12", kMagenta) #finalize8 !!!!!!!!!!!!!!!!!!!!
drawEffs('SingleMu_t100/', "omtf", "12", kCyan)

#  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#80 ZprimeToMuMu, 82 i 84 - all with finzalize 9 and penalty -16
#82 and 84 - the same config, but in the 84 fixed muon matching

drawEffs('ZprimeToMuMu_PU140_0x0006_v3_t84/', "omtf", "12", kBlack)

#drawEffs('ZprimeToMuMu_PU140_v3_t82/', "omtf_patsKB", "12", kBlue) #before fixing the matching for the hight pt, but this fix was very small, the setup is the same as in the t84
#drawEffs('ZprimeToMuMu_PU140_v3_t80/', "omtf_patsKB", "12", kCyan) #by mistake "no matching hit penatly" i finalise 9 was -16
#drawEffs('ZprimeToMuMu_PU140_v3_t84/', "omtf_patsKB", "12", kGreen) #the setup is the same as in the t82, but fix fixed matching for the hight pt
#drawEffs('ZprimeToMuMu_PU140_v3_t98/', "omtf_patsKB", "12", kBlue)

drawEffs('ZprimeToMuMu_PU140_v3_t100/', "omtf", "12", kRed)
#drawEffs('ZprimeToMuMu_PU140_v3_t101/', "omtf", "12", kBlue)
drawEffs('ZprimeToMuMu_PU140_v3_t104/', "omtf", "12", kBlue)
drawEffs('ZprimeToMuMu_PU140_v3_t105/', "omtf", "12", kMagenta)

drawEffs('ZprimeToMuMu_NoPU_v3_t100/', "omtf", "12", kGreen)

eff_c1.cd()
legendEff1.Draw()

eff_c2.cd()
legendEff2.Draw()

eff_c3.cd()
legendEff3.Draw()

# c1.cd(1).Modified()
# c1.cd(1).Update()
# 
# c1.cd(2).Modified()
# c1.cd(2).Update()

#################################
first = True
 
#c2 = TCanvas('canvas_rate_1', 'canvas_rate', 800, 100, 500, 500)

rate_c1 = c1.cd(2)

rate_c1.SetGridx()
rate_c1.SetGridy()
rate_c1.SetLogy()
rate_c1.cd()
rate_c1.SetLeftMargin(0.15)

 
legendRate = TLegend(0.3, 0.65, 0.7, 0.8)

legend = legendRate

 #legend.SetHeader(header.c_str())
# leg -> SetBorderSize(0);
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.03)
legend.SetMargin(0.2)
 
#OMTF 2018+ 
#drawRate('SingleNeutrino_PU200_mtd5_v2_t44/', "omtf", "12", kBlack)
#drawRate('SingleNeutrino_PU200_v2_t44/', "omtf", "12", kBlack)
#drawRate('SingleNeutrino_PU200_v2_t35/', "omtf", 12, kRed)

#drawRate('SingleNeutrino_PU200_v2_t46/', "omtf", "12", kGreen+1) 
#drawRate('SingleNeutrino_PU200_v2_t51/', "omtf", "12", kGreen+1)

drawRate('SingleNeutrino_PU200_v2_t68/', "omtf", "12", kBlack) 

#drawRate('SingleNeutrino_PU200_v2_t55/', "omtf", "12", kRed)

#drawRate('SingleNeutrino_PU200_v2_t56/', "omtf", "12", kBlue)
 
#drawRate('SingleNeutrino_PU200_v2_t41/', "nn_omtf", "12", "0.4", kRed)

#drawRate('SingleNeutrino_PU200_v2_t65/', "omtf", "12", kBlue) 

#drawRate('SingleNeutrino_PU200_mtd5_v2_t44/', "omtf", "12", kBlack)
#drawRate('SingleNeutrino_PU200_v2_t35/', "nn_omtf", 12, kRed)
 
#drawRate('SingleNeutrino_PU200_mtd5_v2_t47/', "omtf", "12", kBlue)

#drawRate('SingleNeutrino_PU200_mtd5_v2_t51/', "omtf", "12", kGreen+1)

#drawRate('SingleNeutrino_PU200_mtd5_v2_t55/', "omtf", "12", kRed)
 
#drawRate('SingleNeutrino_PU200_mtd5_v2_t41/', "nn_omtf", "12", kBlue, "0.4")
#drawRate('SingleNeutrino_PU200_mtd5_v2_t41/', "nn_omtf", "12", kRed, "0.5")

#drawRate('SingleNeutrino_PU200_v2_t41/', "nn_omtf", "12", kBlue, "0.4")
#drawRate('SingleNeutrino_PU200_v2_t44/', "nn_omtf", "12", kRed, "0.5")

#drawRate('SingleNeutrino_PU200_v2_t67/', "nn_omtf", "12", kRed, "0.4")

#drawRate('SingleNeutrino_PU200_v3_t74/', "omtf", "12", kGreen+1)
#drawRate('SingleNeutrino_PU200_v3_t78/', "omtf", "12", kBlue)
#drawRate('SingleNeutrino_PU200_v3_t80/', "omtf", "12", kCyan)
#drawRate('SingleNeutrino_PU200_v3_t82/', "omtf", "12", kRed)
drawRate('SingleNeutrino_PU200_v3_t100/', "omtf", "12", kRed)
drawRate('SingleNeutrino_PU200_v3_t104/', "omtf", "12", kBlue) #good, GoldenPatternResult::finalise9() pdfSum -= 16 (first job failed, for the good one there is no commit

#OMTF 2018+
#drawRate('SingleNeutrino_PU250_v2_t45/', "omtf", "12", kGreen+1)

#drawRate('SingleNeutrino_PU250_v2_t42/', "nn_omtf", "12", kRed)

#drawRate('SingleNeutrino_PU250_v2_t41/', "nn_omtf", "12", kRed)



legend.Draw()

c1.Modified()
c1.Update()
 
c2.Modified()
c2.Update()

c3.Modified()
c3.Update()

if True :
    fileName = makeUniqueFileName("/eos/user/k/kbunkow/public/omtf_nn_plots/", "eff_rate_curves_.png")
    print ("saving as " + fileName)
    c1.SaveAs(fileName)
    c2.SaveAs(fileName.replace("eff_rate_", "eff_eff_"))
    c3.SaveAs(fileName.replace("eff_rate_", "effVsEta_"))

# c3 = TCanvas('canvas_efficiency_rate_1', 'canvas_efficiency_rate_1', 200, 10, 950, 500)
# c3.Divide(2, 1)
# c3.cd(1)
# c3.cd(1).SetGridx()
# c3.cd(1).SetGridy()
# c3.cd(1)
# 
# c3.Modified()
# c3.Update()

raw_input("Press ENTER to exit")

#execfile('ratePlots.py')
