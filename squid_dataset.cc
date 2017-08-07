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
#include	 <cmath>
#include	 <complex>
#include	 <iomanip>

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
#include	 <TRandom.h>
#include 	 <TVirtualFFT.h>

// global stuff
using            namespace std;

double N0 = 3.8e5;            // number of UCNs at t=0 in each 3000 cc cell
double Tau_beta = 885;        // beta decay lifetime, in seconds
double Tau_3 = 500;           // UCN-Helium3 absorption time, in seconds
double Tau_cell = 2000;       // UCN-wall absorption time, in seconds
double T_m = 1000;            // measurement time, in seconds
//  double T_f = 1000;          // cold neutron fill time, in seconds
//  double T_d = 400;           // dead time between cycles
double P3 = 0.98;             // Helium3 initial polarization, fraction
double Pn = 0.98;             // UCN initial polarization, fraction
double Gamma_p = 1.0/20000;   // He3 and UCN depolarization rate, in inverse seconds
double Gamma_T = (1.0/Tau_beta + 1.0/Tau_3 + 1.0/Tau_cell);   // Gamma T for rate function
double epsilon_3 = 0.93;      // detection efficiency for UCN-He3 absorption, fraction
double epsilon_beta = 0.5;    // detection efficiency for beta decay, fraction
double phi_B = 5;             // other background, in Hertz

const   complex<double> i(0.0,1.0);

// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString title, TString command);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  TFile fOut("squid_dataset.root","RECREATE");
  // creating canvas for plotting
  TCanvas *C = new TCanvas("canvas", "canvas", 800, 400);
  C->Divide(2,1);

  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  cout << "Begin generating the fake SQUID signal..." << endl;

  // begin copying of Chris code.
  double alpha = 1;	// alpha is a multiplier to T2 or T_m
			// measurement window alpha*T2 determines sensitivity (standard deviation)
  double S2N_amplitude = 1;     // strength of signal relative to noise of amplitude 1

  double T2m = T_m*1.000;

  double padding = 1e7;         // number of zeros we'll pad on to the end for FFT

  vector < pair < double, double > > timeAndFrequencyOfSQUID;
  double timeStep = 0.01;
  double frequency = 10;	// setting frequency value at t = 0 to 10 Hz
  for(double t = 0; t < alpha*T_m; t = t + timeStep)
  {
    frequency = 10 + t*10;	// frequency offset of 10 Hz, slope of 0.1Hz "slow drift"
    if(fmod(t, 42) < timeStep)	// at each 42s increment, sample and change our frequency in SQUID function
    {
      // with equal probability, change the frequency from previous time "chunk" by 3Hz
      if(engine->Rndm() > 0.5)
      {
        frequency = frequency + 3;
      }
      else
      {
        frequency = frequency - 3;
      }
    }

    timeAndFrequencyOfSQUID.push_back(make_pair(t, frequency));
  }

  vector <double> noise;
  vector <double> squid;
  for(unsigned int j = 0; j < timeAndFrequencyOfSQUID.size(); j++)
  {
    noise.push_back(engine->Gaus(0,1));
    squid.push_back(real(S2N_amplitude*exp(-2*M_PI*timeAndFrequencyOfSQUID[j].second*timeAndFrequencyOfSQUID[j].first*i)*exp(-timeAndFrequencyOfSQUID[j].first/T2m)));
  }

  vector <double> total;
  for(unsigned int j = 0; j < timeAndFrequencyOfSQUID.size(); j++)
  {
    total.push_back(noise[j] + squid[j]);
  }

  cout << "Completed generating the dataset..." << endl;

  // this loop below exists so the TGraph* totalSignal can be referenced correctly
  vector <double> time;
  for(unsigned int i = 0; i < timeAndFrequencyOfSQUID.size(); i++)
  {
    time.push_back(timeAndFrequencyOfSQUID[i].first);
  }

  TGraph* totalSignal = new TGraph(time.size(), &(time[0]), &(total[0]));

  cout << "Performing FFT..." << endl;

  // do the FFT and get the values we care about.
  total.resize(total.size() + padding);
  int fftSize = total.size();
  TVirtualFFT* fft_total = TVirtualFFT::FFT(1, &fftSize, "R2C P K");
  fft_total->SetPoints(&(total[0]));
  fft_total->Transform();

  cout << "Completed FFT..." << endl;

  cout << "Begin displaying results..." << endl;

  TH1* fft_total_hist = 0;
  fft_total_hist = TH1::TransformHisto(fft_total, fft_total_hist, "Re");

  // the fftSize-1 is for an overflow/underflow bin offset that biases the frequency
  TH1D* hFFT = new TH1D("Hist FFT", "Histogram of FFT SQUID + noise"
			, (fftSize - 1) / 2
			, fft_total_hist->GetXaxis()->GetXmin() / ((fftSize-1)*timeStep)
			, fft_total_hist->GetXaxis()->GetXmax() / (2.0*(fftSize-1)*timeStep));

  for(int j = 0; j < hFFT->GetNbinsX(); j++)
  {
    hFFT->SetBinContent(j, fft_total_hist->GetBinContent(j));
  }
  delete fft_total_hist;
  delete fft_total;
  total.clear();

  PlotGraph(C, 1, 1, totalSignal, "SQUID time domain signal", "");
  PlotHist(C, 1, 2, hFFT, "SQUID frequency domain signal", "");

  fOut.Write();

  // Save our plot and print it out as a pdf.
  C -> Print("output_squid_dataset.png");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString title, TString command)
{
  C -> cd(canvasIndex);

  fPlot->SetLineColor(styleIndex % 50);	// only 50 colors in set line color.
//  fPlot->GetYaxis()->SetRangeUser(-40, 40);
  fPlot->GetYaxis()->SetTitle("Erecon error");
  fPlot->GetXaxis()->SetTitle("Evis");
  fPlot->SetTitle(title);

  fPlot->Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  gPlot->SetLineColor(styleIndex);
  gPlot->SetTitle(title);
  gPlot->GetXaxis()->SetTitle("Time (s)");
  gPlot->GetYaxis()->SetTitle("Amplitude (arbitrary)");
//  gPlot->GetYaxis()->SetRangeUser(-40, 40);
  gPlot->Draw(command);
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Frequency (Hz)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Amplitude (arbitrary)");
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

