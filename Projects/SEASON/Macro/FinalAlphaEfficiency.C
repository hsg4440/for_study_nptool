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


void TreatSim(TString filename, double eff[4]){
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("SimulatedTree",filename.Data());
  
  auto d2 = d.Define("DetNbrX","(ROOT::RVec<UShort_t>)SEASON.VGetXE_DetectorNbr()").Define("DetNbrY","(ROOT::RVec<UShort_t>)SEASON.VGetYE_DetectorNbr()").Define("StripNbrX","(ROOT::RVec<UShort_t>)SEASON.VGetXE_StripNbr()").Define("StripNbrY","(ROOT::RVec<UShort_t>)SEASON.VGetYE_StripNbr()");
  
  auto dEnergyRaw = d2.Define("EnergyX", "(ROOT::RVec<Double_t>)SEASON.VGetX_Energy()").Define("EnergyY", "(ROOT::RVec<Double_t>)SEASON.VGetY_Energy()").Define("E_DSSD1_Front", "EnergyX[DetNbrX==1 && EnergyX>7.5] * 1000").Define("E_DSSD1_Back", "EnergyY[DetNbrY==1 && EnergyY>7.5] * 1000").Define("E_Tunnel_Front", "EnergyX[DetNbrX>1 && EnergyX>7.9] * 1000").Define("E_Tunnel_Back", "EnergyY[DetNbrY>1 && EnergyY>7.9] * 1000");//.Filter("DetNbrX.size()==DetNbrY.size()");
  
  auto dEnergyFilter = d2.Define("EnergyX","(ROOT::RVec<Double_t>)SEASON.VGetX_Energy()").Define("EnergyY","(ROOT::RVec<Double_t>)SEASON.VGetY_Energy()").Filter("DetNbrX.size()==DetNbrY.size()").Define("Ediff","1000*abs(EnergyX-EnergyY)").Define("E_DSSD1_Front","EnergyX[DetNbrX==1 && Ediff<40]*1000").Define("E_DSSD1_Back","EnergyY[DetNbrY==1 && Ediff<30]*1000").Define("E_Tunnel_Front","EnergyX[DetNbrX>1 && Ediff<40]*1000").Define("E_Tunnel_Back","EnergyY[DetNbrY>1 && Ediff<30]*1000");//.Filter("DetNbrX.size()==DetNbrY.size()");
  
  
  auto EnergyRaw_DSSD1_Front = dEnergyRaw.Histo1D({"EnergyRaw_DSSD1_Front","DSSD1 Front Energy no filter;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_DSSD1_Front");
  auto EnergyRaw_DSSD1_Back = dEnergyRaw.Histo1D({"EnergyRaw_DSSD1_Back","DSSD1 Back Energy no filter;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_DSSD1_Back");
  auto EnergyRaw_Tunnel_Front = dEnergyRaw.Histo1D({"EnergyRaw_Tunnel_Front","Tunnel Front Energy no filter;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_Tunnel_Front");
  auto EnergyRaw_Tunnel_Back = dEnergyRaw.Histo1D({"EnergyRaw_Tunnel_Back","Tunnel Back Energy no filter;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_Tunnel_Back");
  
  auto EnergyFilter_DSSD1_Front = dEnergyFilter.Histo1D({"EnergyFilter_DSSD1_Front","DSSD1 Front Energy;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_DSSD1_Front");
  auto EnergyFilter_DSSD1_Back = dEnergyFilter.Histo1D({"EnergyFilter_DSSD1_Back","DSSD1 Back Energy;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_DSSD1_Back");
  auto EnergyFilter_Tunnel_Front = dEnergyFilter.Histo1D({"EnergyFilter_Tunnel_Front","Tunnel Front Energy;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_Tunnel_Front");
  auto EnergyFilter_Tunnel_Back = dEnergyFilter.Histo1D({"EnergyFilter_Tunnel_Back","Tunnel Back Energy;Energy (keV);Counts/10keV", 1200, -2000, 10000},"E_Tunnel_Back");
  
  auto Ediff = dEnergyFilter.Histo1D({"Ediff","Ediff;E (keV);Counts",10000,0,10000},"Ediff");
  
  TCanvas *c1 = (TCanvas*) gROOT->FindObject("c1");
  if(c1) delete c1;
  c1 = new TCanvas("c1","c1",700,500);
  c1->cd();
  TH1D* h1 = (TH1D*)EnergyRaw_DSSD1_Front->DrawClone("");
  EnergyFilter_DSSD1_Front->SetLineColor(2);
  TH1D* h1_2 = (TH1D*)EnergyFilter_DSSD1_Front->DrawClone("same");
  // cout << "DSSD1 Front Raw : " << EnergyRaw_DSSD1_Front->GetEntries() << endl;
  // cout << "DSSD1 Front Filter : " << EnergyFilter_DSSD1_Front->GetEntries() << endl;
  // c1->Update();
  
  TCanvas *c2 = (TCanvas*) gROOT->FindObject("c2");
  if(c2) delete c2;
  c2 = new TCanvas("c2","c2",700,500);
  c2->cd();
  TH1D* h2 = (TH1D*)EnergyRaw_DSSD1_Back->DrawClone("");
  EnergyFilter_DSSD1_Back->SetLineColor(2);
  TH1D* h2_2 = (TH1D*)EnergyFilter_DSSD1_Back->DrawClone("same");
  // cout << "DSSD1 Back Raw : " << EnergyRaw_DSSD1_Back->GetEntries() << endl;
  // cout << "DSSD1 Back Filter : " << EnergyFilter_DSSD1_Back->GetEntries() << endl;
  // c2->Update();
  
  TCanvas *c3 = (TCanvas*) gROOT->FindObject("c3");
  if(c3) delete c3;
  c3 = new TCanvas("c3","c3",700,500);
  c3->cd();
  TH1D* h3 = (TH1D*)EnergyRaw_Tunnel_Front->DrawClone("");
  EnergyFilter_Tunnel_Front->SetLineColor(2);
  TH1D* h3_2 = (TH1D*)EnergyFilter_Tunnel_Front->DrawClone("same");
  // cout << "Tunnel Front Raw : " << EnergyRaw_Tunnel_Front->GetEntries() << endl;
  // cout << "Tunnel Front Filter : " << EnergyFilter_Tunnel_Front->GetEntries() << endl;
  // c3->Update();
  
  TCanvas *c4 = (TCanvas*) gROOT->FindObject("c4");
  if(c4) delete c4;
  c4 = new TCanvas("c4","c4",700,500);
  c4->cd();
  TH1D* h4 = (TH1D*)EnergyRaw_Tunnel_Back->DrawClone("");
  EnergyFilter_Tunnel_Back->SetLineColor(2);
  TH1D* h4_2 = (TH1D*)EnergyFilter_Tunnel_Back->DrawClone("same");
  // cout << "Tunnel Back Raw : " << EnergyRaw_Tunnel_Back->GetEntries() << endl;
  // cout << "Tunnel Back Filter : " << EnergyFilter_Tunnel_Back->GetEntries() << endl;
  // c4->Update();
  
  TCanvas *C5 = (TCanvas*) gROOT->FindObject("C5");
  if(C5) delete C5;
  C5 = new TCanvas("C5","C5",700,500);
  C5->cd();
  TH1D* h5 = (TH1D*)Ediff->DrawClone("");
  
  
  eff[0] = EnergyRaw_DSSD1_Front->GetEntries()/100000.;
  eff[1] = EnergyRaw_DSSD1_Back->GetEntries()/100000.;
  eff[2] = EnergyRaw_Tunnel_Front->GetEntries()/100000.;
  eff[3] = EnergyRaw_Tunnel_Back->GetEntries()/100000.;
}

void FinalAlphaEfficiency(){    
  
  ofstream outFile("Results/FinalAlphaEfficiency.txt", ofstream::out | ofstream::app);
  
  outFile << "# Dist DSSD1_Front_Eff DSSD1_Back_Eff Tunnel_Front_Eff Tunnel_Back_Eff" << endl;
  int dist[] = {2,3,4,5,6,7,8,9,10,12,15,20,25,30,40,50};
  for(int i=0;i<16;i++){
    double eff[4];
    TreatSim(TString::Format("FinalAlphaGeometries/Alpha_%imm.root",dist[i]), eff);
    
    outFile << dist[i] << " " << eff[0]  << " " << eff[1]  << " " << eff[2]  << " " << eff[3] << endl;
    cout << dist[i] << " " << eff[0]  << " " << eff[1]  << " " << eff[2]  << " " << eff[3] << endl;
  }
}