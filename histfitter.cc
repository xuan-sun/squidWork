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

// global stuff
using            namespace std;

// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);

// Useful functions for getting the math done.
double DetectionRate(double time);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  // initial, and approximate, experimental values for use in fit.
  double N0 = 3.8e5;            // number of UCNs at t=0 in each 3000 cc cell
  double Tau_beta = 885;        // beta decay lifetime, in seconds
  double Tau_3 = 500;           // UCN-Helium3 absorption time, in seconds
//  double Tau_cell = 2000;     // UCN-wall absorption time, in seconds
//  double T_m = 1000;          // measurement time, in seconds
//  double T_f = 1000;          // cold neutron fill time, in seconds
//  double T_d = 400;           // dead time between cycles
  double P3 = 0.98;             // Helium3 initial polarization, fraction
  double Pn = 0.98;             // UCN initial polarization, fraction
  double Gamma_p = 1.0/20000;   // He3 and UCN depolarization rate, in inverse seconds
  double epsilon_3 = 0.93;      // detection efficiency for UCN-He3 absorption, fraction
  double epsilon_beta = 0.5;    // detection efficiency for beta decay, fraction
  double phi_B = 5;             // other background, in Hertz

  // read in a data file
  TFile fIn("SimData/dataset.root","READ");

  // create a canvas for visualization of histogram + fit function
  TCanvas *C = new TCanvas("canvas", "canvas");

  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  // save our histogram from the read-in file
  TH1D* hEvts = (TH1D*)fIn.Get("Events");

  // Create our fitting function. Set the seed fit values and give the parameters some names for bookkeeping
  double seedFrequencyValue = 0.1;
  TF1* fit = new TF1("rate", "[0]*exp(-[1]*x) + [2]*exp(-[1]*x)*(1 - [3]*cos([4]*x)) + [5]", 0, 10000);
  fit->SetParName(0, "Decay amplitude");
  fit->SetParameter(0, N0*(epsilon_beta/Tau_beta));
  fit->SetParName(1, "Time constant");
  fit->SetParameter(1, Gamma_p);
  fit->SetParName(2, "Oscillatory amplitude");
  fit->SetParameter(2, N0*(epsilon_3/Tau_3));
  fit->SetParName(3, "Total polarization");
  fit->SetParameter(3, P3*Pn);
  fit->SetParName(4, "Frequency");
  fit->SetParameter(4, seedFrequencyValue);
  fit->SetParName(5, "Background offset");
  fit->SetParameter(5, phi_B);

  // fit histogram and plot.
  hEvts->Fit("rate");
  PlotHist(C, 1, 1, hEvts, "Events", "");

  // Save our plot and print it out as a pdf.
  C -> Print("output_histfitter.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
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

