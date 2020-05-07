from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TLegend
from ROOT import kBlack, kBlue, kRed, kGreen, kMagenta

from libPyROOT import TDirectory
import os
import sys
    


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

plotsDir = '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/usercode/L1MuonAnalyzer/test/'

first = True

def drawEff(effFile, type, quality, ptCut, lineColor, legend, pTresh = "0.5") :
    global first
    if type == "nn_omtf" :
        effHist = effFile.Get(type + "_q" + quality + "_pTresh_" + pTresh + "_eta_0.82_1.24_qualityCut_" + quality + "_effOnPtCut_" + ptCut + "_GeV_1")
    else :
        effHist = effFile.Get(type + "_q" + quality + "_eta_0.82_1.24_qualityCut_" + quality + "_effOnPtCut_" + ptCut + "_GeV_1")
        #omtf_q12_eta_0.82_1.24_qualityCut_12_effOnPtCut_20_GeV_1
        
    print (effFile)    
    
    effHist.SetLineColor(lineColor)
    print ("first " + str(first) )
    if first :
        effHist.GetXaxis().SetRangeUser(0, 100)
        effHist.GetYaxis().SetRangeUser(0, 1.05)
        effHist.Draw("hist")
        print ("line 56")
    else:
         effHist.Draw("hist same")   
         print ("line 59")
    legend.AddEntry(effHist)  # , "OMTF", "lep");

###################################################


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

#legendEff1 = TLegend(0.06, 0.8, 0.57, 0.997)
legendEff1 = TLegend(0.21, 0.12, 0.75, 0.31)

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
    print ("first " + str(first) )
    effFile = TFile(plotsDir + fileDir + 'efficiencyPlots.root' )
    effFiles.append(effFile)
    if type == "omtf" :
        eff_c1.cd()
        print (eff_c1.GetName() )
        drawEff(effFile, type, "12", "20", lineColor, legendEff1)
        eff_c2.cd()
        drawEff(effFile, type, "12", "25", lineColor, legendEff2)
        eff_c3.cd()
        drawEff(effFile, type, "12", "40", lineColor, legendEff3)
    if type == "nn_omtf" :
        print (eff_c1.GetName() )
        eff_c1.cd()
        drawEff(effFile, type, "12", "21.5", lineColor, legendEff1, pTresh)
        eff_c2.cd()
        drawEff(effFile, type, "12", "24.5", lineColor, legendEff2, pTresh)
        eff_c3.cd()
        drawEff(effFile, type, "12", "41.5", lineColor, legendEff3, pTresh)
    
    first = False

#OMTF 2018+
#c1.cd(1)
drawEffs('MuFlatPt_PU200_v2_t35/', "omtf", "12", kBlack)

#drawEffs('MuFlatPt_PU200_v2_t46/', "omtf", "12", kGreen)

#OMTF 2018
drawEffs('MuFlatPt_PU200_v2_t47/', "omtf", "12", kGreen)

drawEffs('MuFlatPt_PU200_v2_t41/', "nn_omtf", "12", kBlue, "0.4")
drawEffs('MuFlatPt_PU200_v2_t41/', "nn_omtf", "12", kRed, "0.5")

#drawEffs('MuFlatPt_PU300_v2_t41/', "nn_omtf", "12", kBlue, "0.4")


#drawEff(effFile_t42_pu200, "nn_omtf", "12", "21.5", kRed)

#drawEff(eff_t43_pu300, "nn_omtf", "12", "21.5", kRed)

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
#drawRate('SingleNeutrino_PU200_v2_t35/', "omtf", "12", kBlack)
#drawRate('SingleNeutrino_PU200_v2_t35/', "nn_omtf", 12, kRed)

#drawRate('SingleNeutrino_PU200_v2_t46/', "omtf", "12", kGreen) 
#drawRate('SingleNeutrino_PU200_v2_t47/', "omtf", "12", kGreen)
 
#drawRate('SingleNeutrino_PU200_v2_t41/', "nn_omtf", "12", "0.4", kRed)


drawRate('SingleNeutrino_PU200_mtd5_v2_t44/', "omtf", "12", kBlack)
#drawRate('SingleNeutrino_PU200_v2_t35/', "nn_omtf", 12, kRed)
 
drawRate('SingleNeutrino_PU200_mtd5_v2_t47/', "omtf", "12", kGreen)
 
drawRate('SingleNeutrino_PU200_mtd5_v2_t41/', "nn_omtf", "12", kBlue, "0.4")
drawRate('SingleNeutrino_PU200_mtd5_v2_t41/', "nn_omtf", "12", kRed, "0.5")

#OMTF 2018+
#drawRate('SingleNeutrino_PU250_v2_t45/', "omtf", "12", kGreen)

#drawRate('SingleNeutrino_PU250_v2_t42/', "nn_omtf", "12", kRed)

#drawRate('SingleNeutrino_PU250_v2_t41/', "nn_omtf", "12", kRed)



legend.Draw()

c1.Modified()
c1.Update()
 
c2.Modified()
c2.Update()

if True :
    fileName = makeUniqueFileName("/eos/user/k/kbunkow/public/omtf_nn_plots/", "eff_rate_curves_.png")
    print ("saving as " + fileName)
    c1.SaveAs(fileName)


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
