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
double fullTimeWindow = 1000;	// time window to sample over, in seconds

  double N0 = 3.8e5;            // number of UCNs at t=0 in each 3000 cc cell
  double Tau_beta = 885;        // beta decay lifetime, in seconds
  double Tau_3 = 500;           // UCN-Helium3 absorption time, in seconds
  double Tau_cell = 2000;       // UCN-wall absorption time, in seconds
//  double T_m = 1000;          // measurement time, in seconds
//  double T_f = 1000;          // cold neutron fill time, in seconds
//  double T_d = 400;           // dead time between cycles
  double P3 = 0.98;             // Helium3 initial polarization, fraction
  double Pn = 0.98;             // UCN initial polarization, fraction
  double Gamma_p = 1.0/20000;   // He3 and UCN depolarization rate, in inverse seconds
  double Gamma_T = (1.0/Tau_beta + 1.0/Tau_3 + 1.0/Tau_cell);   // Gamma T for rate function
  double epsilon_3 = 0.93;      // detection efficiency for UCN-He3 absorption, fraction
  double epsilon_beta = 0.5;    // detection efficiency for beta decay, fraction
  double phi_B = 5;             // other background, in Hertz



// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);

// Useful functions for getting the math done.
void FillOneEvent(TH1D* h, TRandom3* factor, double normalizer);
double DetectionRate(double time);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  TFile fOut("dataset.root","RECREATE");
  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas");

  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  int nBins = fullTimeWindow*100;	// this is 1/100 s time bins (i.e. time windows)

  TH1D* hEvts = new TH1D("Events", "Events", nBins, 0, fullTimeWindow);

  double max = 0;
  // getting the max value of the sample function in time steps 100x smaller than bin width
  for(int i = 0; i < nBins*100; i++)
  {
    double valueAti = DetectionRate(0.01*i);
    if(valueAti > max)
    {
      max = valueAti;
    }
  }

  int maxNbEvents = 1e9;
  for(int i = 0; i < maxNbEvents; i++)
  {
    FillOneEvent(hEvts, engine, max);
  }

  // OPTIONAL: adding the histfitter code here to ensure no carry-over issues.
/*  TF1* fit = new TF1("rate",
                Form("[0]*(%f*exp(-%f*x) + %f*exp(-%f*x)*(1 - %f*cos(2*TMath::Pi()*[1]*x + [2])) + %f)",
                        N0*(epsilon_beta/Tau_beta),
                        Gamma_T,
                        N0*(epsilon_3/Tau_3),
                        Gamma_T,
                        P3*Pn,
                        phi_B),
                0, 1000);
  fit->SetParName(0, "Global normalization");
  fit->SetParName(1, "Frequency");
  fit->SetParameter(1, 10);
  fit->SetParName(2, "Phase offset");
  fit->SetParameter(2, 0);
*/
  TF1* fit = new TF1("rate",
		Form("[0]*exp(-[5]*x) + [1]*exp(-[5]*x)*(1 - %f*cos(2*TMath::Pi()*[2]*x + [3])) + [4]" , P3*Pn),
		0, 1000);
  fit->SetParName(0, "Amplitude beta decay");
  fit->SetParName(1, "Amplitude He3 capture");
  fit->SetParName(2, "Frequency of He3");
  fit->SetParName(3, "Phase offset of He3");
  fit->SetParName(4, "Background fraction");
  fit->SetParName(5, "Gamma T (decay timescale)");
  fit->SetParameter(0, N0*(epsilon_beta/Tau_beta));
  fit->SetParameter(1, N0*(epsilon_3/Tau_3));
  fit->SetParameter(2, 10);
  fit->SetParameter(3, 0);
  fit->SetParameter(4, phi_B);
  fit->SetParameter(5, Gamma_T);


  hEvts->Fit("rate", "0");
  TF1* fitResult = hEvts->GetFunction("rate");
  cout << "The Chi-squared is " << fitResult->GetChisquare() << " \n"
       << "NDF is " << fitResult->GetNDF() << " \n"
       << "Chi-squared per degree of freedom is " << fitResult->GetChisquare() / fitResult->GetNDF() << endl;
  PlotHist(C, 1, 1, hEvts, "Events", "");
  fitResult->SetLineWidth(5);
  fitResult->SetLineColor(3);
  fitResult->Draw("SAMEAL");
  // The code above here is a snippet from histfitter that is here for minimizing human error.



//  PlotHist(C, 1, 1, hEvts, "Sampled Events", "");

  fOut.Write();

  // Save our plot and print it out as a pdf.
  C -> Print("output_lightoutput_dataset.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void FillOneEvent(TH1D* h, TRandom3* factor, double normalizer)
{
  double testProb = 1;
  double pdfValue = 0;

  double testTime;

  while(testProb > pdfValue)
  {
    testTime = fullTimeWindow*(factor->Rndm());	// this is the accept-reject part

    pdfValue = DetectionRate(testTime) / normalizer;

    if(pdfValue > 1)
    {
      cout << "WARNING: Monte Carlo pdf value greater than 1. Value is " << pdfValue << endl;
    }

    testProb = factor->Rndm();
  }

  h->Fill(testTime);	// After we get a good time signature that passes the accept-reject, save it.

}

double DetectionRate(double time)
{
/*
  double N0 = 3.8e5;		// number of UCNs at t=0 in each 3000 cc cell
  double Tau_beta = 885;	// beta decay lifetime, in seconds
  double Tau_3 = 500;		// UCN-Helium3 absorption time, in seconds
  double Tau_cell = 2000;	// UCN-wall absorption time, in seconds
//  double T_m = 1000;		// measurement time, in seconds
//  double T_f = 1000;		// cold neutron fill time, in seconds
//  double T_d = 400;		// dead time between cycles
  double P3 = 0.98;		// Helium3 initial polarization, fraction
  double Pn = 0.98;		// UCN initial polarization, fraction
  double Gamma_p = 1.0/20000;	// He3 and UCN depolarization rate, in inverse seconds
  double Gamma_T = (1.0/Tau_beta + 1.0/Tau_3 + 1.0/Tau_cell);	// Gamma T for rate function
  double epsilon_3 = 0.93;	// detection efficiency for UCN-He3 absorption, fraction
  double epsilon_beta = 0.5;	// detection efficiency for beta decay, fraction
  double phi_B = 5;		// other background, in Hertz
*/
  double rateAtTime_time = N0*(epsilon_beta/Tau_beta)*exp(-Gamma_T*time)
			 + N0*(epsilon_3/Tau_3)*exp(-Gamma_T*time)*(1 - P3*Pn*cos(2*M_PI*10*time))
			 + phi_B;

  return rateAtTime_time;
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

