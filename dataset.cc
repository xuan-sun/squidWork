// Standard C and C++ libraries
#include         <vector>
#include         <iostream>
#include         <algorithm>
#include         <functional>
#include         <fstream>
#include         <sstream>
#include         <stdlib.h>
#include         <math.h>
#include         <string.h>
#include         <time.h>
#include         <assert.h>

// Pretty much all the ROOT libraries I have ever used.
#include         <TROOT.h>
#include         <TSystem.h>
#include         <TMath.h>
#include         <TF1.h>
#include         <TGaxis.h>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <TH2F.h>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TObjArray.h>
#include         <TFractionFitter.h>
#include         <TLatex.h>
#include         <TMatrixD.h>
#include	 <TRandom3.h>

using            namespace std;

double fullTimeWindow = 1000;	// global time window to sample over, in seconds

// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);

// Useful functions for getting the math done.
void FillOneEvent(TH1D* h, TRandom3* factor);
double omega(double time);
double phase(double time);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  TCanvas *C = new TCanvas("canvas", "canvas");

  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  TH1D* hEvts = new TH1D("Events", "Events", fullTimeWindow / 10, 0, fullTimeWindow);

  int maxNbEvents = 10000;
  for(int i = 0; i < maxNbEvents; i++)
  {
    FillOneEvent(hEvts, engine);
  }

  PlotHist(C, 1, 1, hEvts, "Sampled Events", "");

  // Save our plot and print it out as a pdf.
  C -> Print("output_dataset.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillOneEvent(TH1D* h, TRandom3* factor)
{
  double testProb = 1;
  double pdfValue = 0;

  double testTime;

  while(testProb > pdfValue)
  {
    testTime = fullTimeWindow*(factor->Rndm());	// this is the accept-reject part

    pdfValue = 0.5*exp(-testTime/fullTimeWindow)*(sin(omega(testTime)*testTime + phase(testTime)) + 1);	// value of the function we are testing against

    if(pdfValue > 1)
    {
      cout << "WARNING: Monte Carlo pdf value greater than 1. Value is " << pdfValue << endl;
    }

    testProb = factor->Rndm();
  }

  h->Fill(testTime);	// After we get a good time signature that passes the accept-reject, save it.

}

double omega(double time)
{
  return 0.1;	// this returns the value of omega as a function of time
}

double phase(double time)
{
  return 0;	// returns value of phi, the phase, as a function of time
}

void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command)
{
  C -> cd(canvasIndex);

  fPlot->SetLineColor(styleIndex % 50);	// only 50 colors in set line color.
  fPlot->GetYaxis()->SetRangeUser(-40, 40);
  fPlot->GetYaxis()->SetTitle("Erecon error");
  fPlot->GetXaxis()->SetTitle("Evis");

  fPlot->Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command)
{
  C -> cd(canvasIndex);
  gPlot->SetLineColor(styleIndex);
  gPlot->GetYaxis()->SetRangeUser(-40, 40);
  gPlot->Draw(command);
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Time (s)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}

