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


void TreatSim(TString filename, double eff[4]){
  
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("SimulatedTree",filename.Data());
  
  auto DSSD1Hit = [](const ROOT::RVec<UShort_t> DetNumberX, const ROOT::RVec<UShort_t> DetNumberY, const ROOT::RVec<Double_t> EnergyX, const ROOT::RVec<Double_t> EnergyY) { 
    if(EnergyX.size()>=EnergyY.size()) return 1000*EnergyX[DetNumberX==1]; 
    else if (EnergyX.size()<EnergyY.size()) return 1000*EnergyY[DetNumberY==1];
  };
  auto TunnelHit = [](const ROOT::RVec<UShort_t> DetNumberX, const ROOT::RVec<UShort_t> DetNumberY, const ROOT::RVec<Double_t> EnergyX, const ROOT::RVec<Double_t> EnergyY) { 
    if(EnergyX.size()>=EnergyY.size()) return 1000*EnergyX[DetNumberX>=2]; 
    else if (EnergyX.size()<EnergyY.size()) return 1000*EnergyY[DetNumberY>=2];
  };
  
  auto d2 = d.Define("EnergyX","(ROOT::RVec<Double_t>)SEASON.VGetX_Energy()");
  auto d3 = d2.Define("EnergyY","(ROOT::RVec<Double_t>)SEASON.VGetY_Energy()");
  auto d4 = d3.Define("DetNumberX","(ROOT::RVec<UShort_t>)SEASON.VGetXE_DetectorNbr()");
  auto d5 = d4.Define("DetNumberY","(ROOT::RVec<UShort_t>)SEASON.VGetYE_DetectorNbr()");
  // auto d5 = d4.Define("StripX","(ROOT::RVec<UShort_t>)SEASON.VGetXE_StripNbr()");
  // auto d6 = d5.Define("StripY","(ROOT::RVec<UShort_t>)SEASON.VGetYE_StripNbr()");
  auto d6 = d5.Define("EDSSD1",DSSD1Hit,{"DetNumberX","DetNumberY","EnergyX","EnergyY"});
  auto d7 = d6.Define("ETunnel",TunnelHit,{"DetNumberX","DetNumberY","EnergyX","EnergyY"});
  
  
  auto h_DSSD1 = d7.Histo1D({"Energy_DSSD1","Electron Energy in DSSD1;DSSD1 Energy (keV);Counts/keV",1000,0,1000},"EDSSD1");
  auto h_Tunnel = d7.Histo1D({"Energy_DSSD2","Electron Energy in Tunnel;Tunnel Energy (keV);Counts/keV",1000,0,1000},"ETunnel");
  
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
  FitFunction->SetParameters(100000,50,2,0,0,0,0);
  TFitResultPtr fitRes = h1->Fit("FitFunction", "SMQE", "SAME", 30, 70);
  //cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
  double IDSSD1_50 = fitRes->Parameter(0);
  double IDSSD1_50_err = fitRes->ParError(0);
  
  // c1->Update();
  // 
  // FitFunction->SetParameters(100000,150,2,0,0,0,0);
  // fitRes = h1->Fit("FitFunction", "SMQE", "SAME", 130, 170);
  // 
  // double IDSSD1_150 = fitRes->Parameter(0);
  // 
  // c1->Update();
  // 
  // FitFunction->SetParameters(100000,250,2,0,0,0,0);
  // fitRes = h1->Fit("FitFunction", "SMQE", "SAME", 230, 270);
  // 
  // double IDSSD1_250 = fitRes->Parameter(0);
  // 
  // c1->Update();
  // 
  // FitFunction->SetParameters(100000,500,2,0,0,0,0);
  // fitRes = h1->Fit("FitFunction", "SMQE", "SAME", 480, 520);
  // 
  // double IDSSD1_500 = fitRes->Parameter(0);
  
  c1->Update();
  
  c2->cd();
  FitFunction->SetParameters(100000,50,2,0,0,0,0);
  fitRes = h2->Fit("FitFunction", "SMQE", "SAME", 30, 70);
  
  double ITunnel_50 = fitRes->Parameter(0);
  double ITunnel_50_err = fitRes->ParError(0);
  
  c2->Update();
  
  // FitFunction->SetParameters(100000,150,2,0,0,0,0);
  // fitRes = h2->Fit("FitFunction", "SMQE", "SAME", 130, 170);
  // 
  // double ITunnel_150 = fitRes->Parameter(0);
  // 
  // c2->Update();
  // 
  // FitFunction->SetParameters(100000,250,2,0,0,0,0);
  // fitRes = h2->Fit("FitFunction", "SMQE", "SAME", 230, 270);
  // 
  // double ITunnel_250 = fitRes->Parameter(0);
  // 
  // c2->Update();
  // 
  // FitFunction->SetParameters(100000,500,2,0,0,0,0);
  // fitRes = h2->Fit("FitFunction", "SMQE", "SAME", 480, 520);
  // 
  // double ITunnel_500 = fitRes->Parameter(0);
  // 
  // c2->Update();
  // 
  //cout << IDSSD1_50/400000. << " " << IDSSD1_150/400000. << " " << IDSSD1_250/400000. << " " << IDSSD1_500/400000. << endl;
  //cout << ITunnel_50/400000. << " " << ITunnel_150/400000. << " " << ITunnel_250/400000. << " " << ITunnel_500/400000. << endl;
  
  
  
  eff[0] = IDSSD1_50/100000.;
  eff[1] = ITunnel_50/100000.;
  eff[2] = IDSSD1_50_err/100000.;
  eff[3] = ITunnel_50_err/100000.;
  // 
  // cout << eff[0] << " " << eff[2] << endl;
  // cout << eff[1] << " " << eff[3] << endl;
  // cout << eff[0]+eff[1] << " " << TMath::Sqrt(eff[2]*eff[2]+eff[3]*eff[3]) << endl;
}


void ElectronEfficiency(){
  
  
  ofstream outFile("Results/MultipleElectronsMultipleGeometriesEfficiency.txt", ofstream::out | ofstream::app);
  
  outFile << "#Nemitted Dist DSSD1_Eff Err DSSD2_Eff Err Total_Eff Err " << endl;
  double eff[4];
  int dist[7] = {2,3,4,5,7,10,15};
  int NPartMax=6;
  for(int i=0;i<7;i++){
    TreatSim(TString::Format("MultipleEmissionMultipleGeometries/0Alpha1Electron_%imm.root",dist[i]), eff);
    outFile << 1 << " " << dist[i] << " " << eff[0] << " " << eff[2] << " " << eff[1] << " " << eff[3] << " " << eff[0]+eff[1] << TMath::Sqrt(eff[2]*eff[2]+eff[3]*eff[3]) << endl;
    cout << 1 << " " << dist[i] << " " << eff[0] << " " << eff[2] << " " << eff[1] << " " << eff[3] << " " << eff[0]+eff[1] << TMath::Sqrt(eff[2]*eff[2]+eff[3]*eff[3]) << endl;
  }
  for(int NPart=2;NPart<=NPartMax;NPart++){
    for(int i=0;i<7;i++){
      TreatSim(TString::Format("MultipleEmissionMultipleGeometries/1Alpha%iElectron_%imm.root",NPart-1,dist[i]), eff);
      outFile << NPart << " " << dist[i] << " " << eff[0] << " " << eff[2] << " " << eff[1] << " " << eff[3] << " " << eff[0]+eff[1] << TMath::Sqrt(eff[2]*eff[2]+eff[3]*eff[3]) << endl;
      cout << NPart << " " << dist[i] << " " << eff[0] << " " << eff[2] << " " << eff[1] << " " << eff[3] << " " << eff[0]+eff[1] << TMath::Sqrt(eff[2]*eff[2]+eff[3]*eff[3]) << endl;
    }
  }
}