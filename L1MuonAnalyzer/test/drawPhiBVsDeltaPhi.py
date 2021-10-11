from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D
from ROOT import gROOT
from ROOT import gStyle
from ROOT import TLegend
from ROOT import kBlack, kBlue, kRed, kGreen, kMagenta, kCyan
from array import array

#from libPyROOT import TDirectory
from ROOT import TDirectory
import os
import sys
from matplotlib import legend

gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);

patternCnt = 56

ptBins = []
ptBins.append("q = 1  0 - 3.5 GeV    ")
ptBins.append("q = 1  3.5 - 4 GeV    ")
ptBins.append("q = 1  4 - 4.5 GeV    ")
ptBins.append("q = 1  4.5 - 5 GeV    ")
ptBins.append("q = 1  5 - 6 GeV      ")
ptBins.append("q = 1  6 - 7 GeV      ")
ptBins.append("q = 1  7 - 8 GeV      ")
ptBins.append("q = 1  8 - 10 GeV     ")
ptBins.append("q = 1  10 - 12 GeV    ")
ptBins.append("q = 1  12 - 14 GeV    ")
ptBins.append("q = 1  14 - 16 GeV    ")
ptBins.append("q = 1  16 - 18 GeV    ")
ptBins.append("q = 1  18 - 20 GeV    ")
ptBins.append("q = 1  20 - 22 GeV    ")
ptBins.append("q = 1  22 - 24 GeV    ")
ptBins.append("q = 1  24 - 26 GeV    ")
ptBins.append("q = 1  26 - 28 GeV    ")
ptBins.append("q = 1  28 - 30 GeV    ")
ptBins.append("q = 1  30 - 35 GeV    ")
ptBins.append("q = 1  35 - 40 GeV    ")
ptBins.append("q = 1  40 - 45 GeV    ")
ptBins.append("q = 1  45 - 50 GeV    ")
ptBins.append("q = 1  50 - 60 GeV    ")
ptBins.append("q = 1  60 - 70 GeV    ")
ptBins.append("q = 1  70 - 80 GeV    ")
ptBins.append("q = 1  80 - 100 GeV   ")
ptBins.append("q = 1  100 - 200 GeV  ")
ptBins.append("q = 1  200 - 10000 GeV")
ptBins.append("q =-1  0 - 3.5 GeV    ")
ptBins.append("q =-1  3.5 - 4 GeV    ")
ptBins.append("q =-1  4 - 4.5 GeV    ")
ptBins.append("q =-1  4.5 - 5 GeV    ")
ptBins.append("q =-1  5 - 6 GeV      ")
ptBins.append("q =-1  6 - 7 GeV      ")
ptBins.append("q =-1  7 - 8 GeV      ")
ptBins.append("q =-1  8 - 10 GeV     ")
ptBins.append("q =-1  10 - 12 GeV    ")
ptBins.append("q =-1  12 - 14 GeV    ")
ptBins.append("q =-1  14 - 16 GeV    ")
ptBins.append("q =-1  16 - 18 GeV    ")
ptBins.append("q =-1  18 - 20 GeV    ")
ptBins.append("q =-1  20 - 22 GeV    ")
ptBins.append("q =-1  22 - 24 GeV    ")
ptBins.append("q =-1  24 - 26 GeV    ")
ptBins.append("q =-1  26 - 28 GeV    ")
ptBins.append("q =-1  28 - 30 GeV    ")
ptBins.append("q =-1  30 - 35 GeV    ")
ptBins.append("q =-1  35 - 40 GeV    ")
ptBins.append("q =-1  40 - 45 GeV    ")
ptBins.append("q =-1  45 - 50 GeV    ")
ptBins.append("q =-1  50 - 60 GeV    ")
ptBins.append("q =-1  60 - 70 GeV    ")
ptBins.append("q =-1  70 - 80 GeV    ")
ptBins.append("q =-1  80 - 100 GeV   ")
ptBins.append("q =-1  100 - 200 GeV  ")
ptBins.append("q =-1  200 - 10000 GeV")

if len(ptBins) != patternCnt :
    print("len(ptBins) != patternCnt")
    exit(1)


inputFiles = []
canvases = []

legend1 = TLegend(0.8, 0.01, 0.99, 0.9)
#legend1.SetFillStyle(0)

patternFileDir = '/home/kbunkow/CMSSW/CMSSW_12_1_0_pre3/src/L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/'
displacePatternFile = TFile(patternFileDir + 'Patterns_dispalced_test_displHighPt_t10_200files.root' )
inputFiles.append(displacePatternFile)
displacePatternFile.ls()

layersStatDirDispl = displacePatternFile.Get("layerStats")
#layersStatDirDispl.ls()

refLayer_layer = 'refLayer_0_Layer_7'

displHist = layersStatDirDispl.Get("histLayerStat_PatNum_18_" + refLayer_layer)

c1 = TCanvas('canvas_1' , 'canvas_1', 200, 10, 900, 900)
c1.SetLeftMargin(0.13)
c1.SetRightMargin(0.22)
c1.cd().SetGridx()
c1.cd().SetGridy()
canvases.append(c1)

displHist.SetMarkerStyle(6)
displHist.SetMarkerColor(kBlack)

if "Layer_1" in refLayer_layer or "Layer_3" in refLayer_layer:
    displHist.GetYaxis().SetRangeUser(-300, 300)
else : 
    displHist.GetYaxis().SetRangeUser(-60, 60)

#displHist.GetYaxis().SetTitle("#phi_{B2} - #phi_{B2extrpol}")
displHist.GetYaxis().SetTitle("#Delta#phi_{2} - #Delta#phi_{2extrpol}")

displHist.Draw()


promptPatternFile = TFile(patternFileDir + 'Patterns_dispalced_test_allPt_t10.root' )
inputFiles.append(promptPatternFile)
promptPatternFile.ls()

layersStatDirPromt = promptPatternFile.Get("layerStats")

col = 2

step = 1
for iPat in range(0, patternCnt, step) :  
    promtHist = layersStatDirPromt.Get("histLayerStat_PatNum_" + str(iPat) + "_" + refLayer_layer)
    print(promtHist.GetName())
    
    promtHist.SetMarkerColor(col)
    promtHist.SetFillColor(col)
    promtHist.SetMarkerStyle(1)
    promtHist.Draw("same")
    
    legend1.AddEntry(promtHist, ptBins[iPat], "f") 
    col += step
    
    if iPat == patternCnt/2 -step:
        col = 2

displHist.Draw("same")
legend1.Draw()

c1.Modified()
c1.Update()

input("Press ENTER to exit")
 