from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1D, TEfficiency, TH2D, THStack, TLegend
from ROOT import gROOT
from ROOT import gStyle
from ROOT import kBlack, kBlue, kRed, kGreen, kMagenta
from libPyROOT import TDirectory
import os
import sys
from collections import namedtuple
#from __builtin__ import None

gStyle.SetOptStat(0)


#fileDir = '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/crab/'
fileDir = '/afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_11_x_x_l1tOfflinePhase2/CMSSW_11_1_3/src/L1Trigger/L1TMuonOverlapPhase1/test/crab/'


candCategories = ('promptMuons', 'nonPromptMuons', 'notMatched') #, 'muonsFromPions', 'muonsFromKaons'

canvases = []
rateCumuls = []
rateCumulStacks = []
histFiles = {}
legends = []

acceptedMuonsPts = []

def getRate(dir, ptCutGev,name):
    results = []
    for obj in dir.GetListOfKeys():
        if "candPt" in obj.GetName() :
            candPt = obj.ReadObj()
            #candPt.SetLineColor(lineColor)
            #candPt.Draw("")
            
            candPt.Sumw2(False);
            candPt_rateCumul = candPt.GetCumulative(False, "_" + name);
            
            ptCutBin = candPt_rateCumul.GetXaxis().FindBin(ptCutGev)        
            rateOnThresh = candPt_rateCumul.GetBinContent(ptCutBin) 
            candPt_rateCumul.GetXaxis().SetTitle("cand pT thresh [GeV]")
            results.append(rateOnThresh)
            results.append(candPt_rateCumul)
        if "_ptGenVsPtCand" in obj.GetName() :
            ptGenVsPtCand = obj.ReadObj()
            
            ptCutBin = ptGenVsPtCand.GetYaxis().FindBin(ptCutGev)  
            acceptedMuonsPt = ptGenVsPtCand.ProjectionX(ptGenVsPtCand.GetName() + "_" + name + "_ProjectionX", ptCutBin, -1)
            acceptedMuonsPt.SetTitle(ptGenVsPtCand.GetTitle() + " muon spectrum, L1 cut " + str(ptCutGev) + " " + name )
            results.append(acceptedMuonsPt)
            
    return results 
                
def readRateFile(version, algoName, ptCutGev):
    folder = 'crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_v2_' + version
    if version == "t74" or version == "t78" or version == "t80" or version == "t100"  or version == "t104":
        folder = 'crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_v3_' + version
        
    if folder in histFiles:
        histFile = histFiles.get(folder)
    else :
        histFile = TFile( fileDir + folder + '/results/omtfAnalysis2.root' )
        histFiles[folder] = histFile
    
    print ("\nreading file " + histFile.GetName())
    
    rateDir = histFile.Get("L1MuonAnalyzerOmtf/rate")
    #rateDir.ls()
    
    algoDir =  None
    #algoDir.ls()
    
    for obj in rateDir.GetListOfKeys():
        if algoName == obj.GetName() :
            algoDir = obj.ReadObj()
    
    if algoDir is None:
        print("coiuld not find algo " + algoName)
        return
    
    print(" algoDir %s ptCutGev %f " % (algoDir.GetName(), ptCutGev) )
    
    name = version + "_" + algoDir.GetName()
    title = version + " " + algoDir.GetName().replace("_", " ")
    
    c1 = TCanvas('canvas_' + name, title, 200, 10, 700 *2, 700)
    canvases.append(c1)
    c1.Divide(2, 1)
    c1.cd(1).SetGridx()
    c1.cd(1).SetGridy()
    c1.SetLogy()
    print ('created canvas ' + c1.GetName())
    
    legend = TLegend(0.5, 0.6, 0.8, 0.8)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.027)
    legend.SetHeader(title)
    legends.append(legend)
    
    rateResult = getRate(algoDir, ptCutGev, name + "_total")
    rateOnThresh = rateResult[0]
    
    candPt_rateCumul =  rateResult[1] 
    rateCumuls.append(candPt_rateCumul)
    candPt_rateCumul.SetLineColor(kBlack)
    candPt_rateCumul.GetYaxis().SetTitle("#events")
    candPt_rateCumul.Draw("hist")
    candPt_rateCumul.GetYaxis().SetRangeUser(10, 20000)
    print ("added candPt_rateCumul " + candPt_rateCumul.GetName() )
    
    legend.AddEntry(candPt_rateCumul, "total rate")
    
    print ("total rateOnThresh %f" % (rateOnThresh) )
    
    rateCumulStack = THStack("rateCumulStack_" + name, title )
    rateCumulStacks.append(rateCumulStack) 
    
     
    totalRate = 0
    color = 2
    
    thisAcceptedMuonsPts = []
    for candCategory in candCategories :
        #categoryDir = histFile.Get("L1MuonAnalyzerOmtf/rate/" +  algoName + "/" + candCategory)
        categoryDir = algoDir.Get(candCategory)
        
        print ("categoryDir " + categoryDir.GetName() )
        
        rateResult  = getRate(categoryDir, ptCutGev, name + "_" + candCategory)
        rateOnThresh = rateResult[0]
        
        rateVsCategory[candCategory].Fill(title + " " + str(ptCutGev) +" GeV", rateOnThresh) 
        rateVsCategory[candCategory].SetFillColor(color)
        
        totalRate = totalRate + rateOnThresh
        print("%s rate %f " % (categoryDir.GetName(), rateOnThresh) )
                
        candPt_rateCumul =  rateResult[1]   
        #candPt_rateCumul.SetFillColor(color)
        candPt_rateCumul.SetLineColor(color)

        #rateCumulStack.Add(candPt_rateCumul)
        rateCumuls.append(candPt_rateCumul)
        print ("added candPt_rateCumul " + candPt_rateCumul.GetName() )
        
        candPt_rateCumul.Draw("hist same")
        legend.AddEntry(candPt_rateCumul, candCategory)
        
        if rateResult.__len__() > 2 :
            acceptedMuonsPt = rateResult[2]
            acceptedMuonsPt.SetLineColor(color)
            acceptedMuonsPts.append(acceptedMuonsPt)
            thisAcceptedMuonsPts.append(acceptedMuonsPt)
        
        color += 1
        if color == 5:
            color = 6
            
    legend.Draw()
    #rateCumulStack.Draw()    
    
    c1.cd(2).SetGridx()
    c1.cd(2).SetGridy()
    
    first = True
    for acceptedMuonsPt in thisAcceptedMuonsPts :
        print "Drawing", acceptedMuonsPt.GetTitle()
        if first :
            acceptedMuonsPt.GetYaxis().SetRangeUser(0, 15)
            acceptedMuonsPt.Draw("hist")
        else :
            acceptedMuonsPt.Draw("hist same")
        first = False 
    
    
    c1.Modified()
    c1.Update(); 
    print ("totalRate %f" % (totalRate) )



VersionAlgo = namedtuple("VersionAlgo", "version algoName ptCutGev")

toCompareList = ( 
    VersionAlgo("t68", "omtf_q12", 10) , 
    #VersionAlgo("t65", "omtf_q8", 5), 
    #VersionAlgo("t65", "omtf_q12", 18), 
    
    #VersionAlgo("t74", "omtf_q8", 18), 
    #VersionAlgo("t74", "omtf_q12", 18), 


    #VersionAlgo("t78", "omtf_q12", 18), 
    
    #VersionAlgo("t80", "omtf_q12", 18),
     
    #VersionAlgo("t68", "omtf_q1", 20) , 
    #VersionAlgo("t67", "nn_omtf_q12_pTresh_0.4", 22),
    #VersionAlgo("t67", "nn_omtf_q12_pTresh_0.5", 22),
    
    
    VersionAlgo("t100", "omtf_q12", 9),
    #VersionAlgo("t104", "omtf_q12", 20),
    
    #VersionAlgo("t68", "nn_omtf_q12_pTresh_0.5", 22),
    #VersionAlgo("t80", "nn_omtf_q12_pTresh_0.4", 22),
    #VersionAlgo("t80", "nn_omtf_q12_pTresh_0.5", 22),
    
    #VersionAlgo("t66", "nn_omtf_q12_pTresh_0.5", 22), #no mathcin for nn 
    
    VersionAlgo("t67", "nn_omtf_q12_pTresh_0.4", 10), #this is good nn 
    #VersionAlgo("t67", "nn_omtf_q12_pTresh_0.5", 22), #this is good nn 
    )

# toCompareList = ( 
#     VersionAlgo("t68", "omtf_q1", 20) , 
#     VersionAlgo("t65", "omtf_q1", 18), 
#     #VersionAlgo("t68", "omtf_q1", 20) , 
#     VersionAlgo("t67", "nn_omtf_q1_pTresh_0.4", 22),
#     VersionAlgo("t67", "nn_omtf_q1_pTresh_0.5", 22),
#     )


#stackQ12 = THStack("hs","");
#histQ12 = []


# for candCategory in candCategories :
#     
#     histQ12.append(candPt_rateCumul)
#def makeHistograms(stack, hists):

rateVsCategory = {}
rateVsCategoryStack = THStack("rateVsCategoyStack", "rateVsCategoryStack")
rateVsCategoryStackLegend = TLegend(0.65, 0.8, 0.93, 0.95)
rateVsCategoryStackLegend.SetBorderSize(0)
rateVsCategoryStackLegend.SetTextSize(0.027)


    
for candCategory in candCategories :
    rateVsCategory[candCategory] = TH1D("rateVsCategory_" + candCategory, candCategory, toCompareList.__len__(), 0, toCompareList.__len__())
    rateVsCategory[candCategory].GetYaxis().SetTitle("#events")
    rateVsCategoryStack.Add(rateVsCategory[candCategory])
    rateVsCategoryStackLegend.AddEntry(rateVsCategory[candCategory], candCategory)
    
for toCompare in toCompareList:
    readRateFile(toCompare.version, toCompare.algoName, toCompare.ptCutGev)
    
rateVsCategoryCanvas = TCanvas('rateVsCategoryCanvas', "rateVsCategoryCanvas", 200, 10, 700, 900)   
rateVsCategoryCanvas.SetRightMargin(0.25)
rateVsCategoryCanvas.SetBottomMargin(0.16)
rateVsCategoryStack.Draw("hist text")
rateVsCategoryStack.GetYaxis().SetTitle("#events")
#rateVsCategoryStack.GetYaxis().SetRangeUser(0, 10000)
rateVsCategoryStackLegend.Draw()
#rateVsCategoryStack.SetLeftMargin(0.02)

#gROOT.GetListOfCanvases().Draw()

#outFile.Close()

for canvas in canvases :
    canvas.Modified()
    canvas.Update(); 
    print ("updated canvas " + canvas.GetName() )
    
for rateCumul in rateCumuls :     
    print (" rateCumul "  + rateCumul.GetName() )

raw_input("Press ENTER to exit")    
