// Drawing template (originally used for TRG-17-001)
// Author: O. Davignon (CERN)
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMarker.h>
#include <TLine.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <string>

TCanvas* CreateCanvas(TString CanvasName = "myPlot", bool LogY = false, bool Grid = true)
{
  TCanvas* c = new TCanvas(CanvasName.Data(),CanvasName.Data(),800,800);
  c->SetLeftMargin(0.11);
  if(Grid)
    {
      c->SetGrid();
    }
  if(LogY)
    {
      c->SetLogy();
    }
  return c;
}

void DrawPrelimLabel(TCanvas* c)
{
  c->cd();

  TLatex tex;
  tex.SetTextSize(0.03);
  tex.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS}");
  tex.Draw("same");

  return;
}

void DrawLumiLabel(TCanvas* c, TString Lumi = "35.9")
{
  c->cd();

  TLatex tex;
  tex.SetTextSize(0.035);
  TString toDisplay = Lumi + " fb^{-1} (13 TeV)";
  tex.DrawLatexNDC(0.66,0.91,toDisplay.Data());
  tex.Draw("same");

  return;
}

void DrawCmsSimulationLabel(TVirtualPad* c)
{
  c->cd();

  TLatex tex;
  tex.SetTextSize(0.03);

  tex.DrawLatexNDC(0.11,0.91,"#scale[1.5]{CMS} Phase-2 Simulation");//typically for Phase-2
  //TString toDisplay = "CMS Simulation, Phase-2, <PU> = " + pu + ", 14 TeV";
  //tex.DrawLatexNDC(0.27,0.91,toDisplay.Data());
  tex.Draw("same");

  return;
}

void DrawPuLabel(TVirtualPad* c, TString pu = "200 PU" ) //14 TeV,
{
  c->cd();

  TLatex tex;
  tex.SetTextSize(0.035);
  tex.SetTextAlign(31);
  TString toDisplay = pu;//typically for Phase-2
  // TString toDisplay = Lumi + " fb^{-1} (13 TeV)";//typically for Phase-1
  tex.DrawLatexNDC(0.90,0.91,toDisplay.Data());
  tex.Draw("same");

  return;
}

void DrawLabel(TVirtualPad* c, std::string label)
{
  c->cd();

  TLatex tex;
  tex.SetTextSize(0.03);
  TString toDisplay = label.c_str();
  tex.DrawLatexNDC(0.27,0.45,toDisplay.Data());
  tex.Draw("same");

  return;
}

void SaveCanvas(TVirtualPad* c, std::string plotsDir, std::string PlotName = "myPlotName")
{
  c->cd();
  c->SaveAs( (plotsDir + "/" + PlotName + ".pdf").c_str() );
  //c->SaveAs( (plotsDir + "/" + PlotName + ".root").c_str() );
  c->SaveAs( (plotsDir + "/" + PlotName + ".png").c_str() );

  return;
}
