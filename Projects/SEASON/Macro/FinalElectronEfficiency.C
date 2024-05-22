#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>

#include "TH1.h"
#include "TMath.h"
#include "TSEASONData.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"

using namespace ROOT;
using namespace std;

Double_t fitfunction(Double_t *x, Double_t *par){
  Double_t E = x[0];
  Double_t temp = (E-par[1])/par[2];
  Double_t G = par[0]/(sqrt(2*TMath::Pi())*par[2])*exp(-0.5*temp*temp);
  Double_t r = 0.5-0.5*TMath::TanH(temp);
  return G + r*(par[3]*E + par[4]) + (1-r)*(par[5]*E + par[6]);
}

Double_t ElectronPeak(Double_t *x, Double_t *par){
  Double_t I = par[0], mean = par[1], sigma = par[2], b_rate = par[3];
  Double_t value = (x[0]-mean)/(TMath::Sqrt(2)*sigma);
  Double_t amp = I/(sigma*TMath::Sqrt(2.*TMath::Pi()));
  Double_t G = (1-b_rate) * amp * TMath::Exp(-value*value);
  Double_t backscattering = b_rate * I * TMath::Erfc(value) / (2*mean);
  return G + backscattering;
}

void TreatSim(TString filename, double eff[4]){
  
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("SimulatedTree",filename.Data());
  
  auto DSSD1Hit = [](const ROOT::RVec<UShort_t> DetNumber, const ROOT::RVec<Double_t> Energy) { return 1000*Energy[DetNumber==1]; };
  auto TunnelHit = [](const ROOT::RVec<UShort_t> DetNumber, const ROOT::RVec<Double_t> Energy) { return 1000*Energy[DetNumber>=2]; };
  
  auto d2 = d.Define("Energy","(ROOT::RVec<Double_t>)SEASON.VGetY_Energy()").Define("DetNumber","(ROOT::RVec<UShort_t>)SEASON.VGetYE_DetectorNbr()").Define("EDSSD1",DSSD1Hit,{"DetNumber","Energy"}).Define("ETunnel",TunnelHit,{"DetNumber","Energy"});
  auto h_DSSD1 = d2.Histo1D({"Energy_DSSD1","Electron Energy in DSSD1;DSSD1 Energy (keV);Counts/keV",40,80,120},"EDSSD1");
  auto h_Tunnel = d2.Histo1D({"Energy_DSSD2","Electron Energy in DSSD2;DSSD2 Energy (keV);Counts/keV",40,80,120},"ETunnel");
  
  TCanvas *c1 = (TCanvas*) gROOT->FindObject("c1");
  if(c1) delete c1;
  c1 = new TCanvas("c1","c1",700,500);
  TH1D* h1 = (TH1D*)h_DSSD1->DrawClone();
  
  TCanvas *c2 = (TCanvas*) gROOT->FindObject("c2");
  if(c2) delete c2;
  c2 = new TCanvas("c2","c2",700,500);
  TH1D* h2 = (TH1D*)h_Tunnel->DrawClone();
  
  c1->cd();
  TF1 *FitFunction = new TF1("FitFunction", fitfunction, 0, 1000, 7);
  FitFunction->SetParameters(50000,100,3,0,0,0,0);
  // TF1 *FitFunction = new TF1("FitFunction", ElectronPeak, 0, 1000, 4);
  // FitFunction->SetParameters(50000,100,3,0.5);
  TFitResultPtr fitRes = h1->Fit("FitFunction", "SMQE", "SAME", 80, 120);
  //cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
  double IDSSD1_100 = fitRes->Parameter(0);
  double IDSSD1_100_err = fitRes->ParError(0);
  
  c1->Update();
  
  c2->cd();
  FitFunction->SetParameters(50000,100,3,0,0,0,0);
  // FitFunction->SetParameters(50000,100,3,0.5);
  fitRes = h2->Fit("FitFunction", "SMQE", "SAME", 80, 120);
  
  double ITunnel_100 = fitRes->Parameter(0);
  double ITunnel_100_err = fitRes->ParError(0);
  
  c2->Update();
  
  eff[0] = IDSSD1_100/100000.;
  eff[1] = IDSSD1_100_err/100000.;
  eff[2] = ITunnel_100/100000.;
  eff[3] = ITunnel_100_err/100000.;
}


void ElectronEfficiency(){
  ofstream outFile("Results/FinalGeometryElectronsEfficiency.txt", ofstream::out | ofstream::app);
  
  outFile << "# Dist DSSD1_Eff DSSD2_Eff " << endl;
  int dist[] = {2,3,4,5,6,7,8,9,10,12,15,20,25,30,40,50};
  for(int i=0;i<16;i++){
    double eff[4];
    TreatSim(TString::Format("FinalElectronGeometries/Electron_%imm.root",dist[i]), eff);
    
    outFile << dist[i] << " " << eff[0]  << " " << eff[2] << endl;
    cout << dist[i] << " " << eff[0]  << " " << eff[2] << endl;
  }
}