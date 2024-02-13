#ifndef recuit_h
#define recuit_h 1

#include<iostream>
#include<vector>

// ROOT
#include "TChain.h"
#include "TGraph.h"
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TSpectrum2.h"
#include "TSpectrum.h"
#include "TCutG.h"
#include "TF1.h"

// NPTOOL
#include "NPReaction.h"
#include "TCATSPhysics.h"
#include "TZDDPhysics.h"
#include "TTACPhysics.h"
#include "TMust2Physics.h"

  
class MCRMethod
{
 public:
  MCRMethod();
  MCRMethod(double lambd, double ini, double final,  int n1);
  ~MCRMethod();

  void MetroHast(); //Calcul de la distance minimale et du chemin pour le problème du voyageur de commerce
  double read(); // Affichage
  void plot();
  void RandomStep();
  double MinimizationFunction();
  void Init();
  float ProjectOnCats(unsigned int i, double x1, double x2);
  float ProjectOnTarget(unsigned int i, double x1, double x2);
  float ProjectAndEx(unsigned int i,string side, double x1, double x2);
  void NappeBis();
  bool condition(unsigned int);
  void Unallocate();
  void Nappe();

 private:

  double Lambda = 0.98;
  std::vector<double> dist;
  std::vector<Double_t> pos[4];
  unsigned int it;
  unsigned int itmin;
  TString runmask[2];
  
  
  
  double Step;
  double TemperatureInitiale;
  double Temperature;
  double TemperatureFinale;
  double DistRatio[2];
  double CATSPosXY[4] = {8,-6,-2,-1}; // CATS1X, 1Y, 2X, 2Y;
  double CATSPosXYCan[4] = {8,-6,-2,-1}; // CATS1X, 1Y, 2X, 2Y;
  double CATSPosXYLim[8] = {-30,30,-30,30,-30,30,-30,30}; // CATS1X min max, 1Y etc, 2X, 2Y;
  double CATSPosZ[2] = {-1587.1,-1090.1};
  double MASKPosZ[2] = {-1732.1,-1235.1};

  private:
  TChain *Chain[2];
  TTreeReader *TreeReader[2];
  TString path = "./ssd/";
  TFile * f_cut_Cr = new TFile("./CUT_Cr.root");
  TCutG* cut_Cr =  (TCutG*) f_cut_Cr->FindObjectAny("CUTCr");
  TFile * f_cut_deuton = new TFile("./CUT_deuton.root");
  TCutG* cut_deuton =  (TCutG*) f_cut_deuton->FindObjectAny("CUT_deuton");
  TTreeReaderValue<TCATSPhysics> *CATSPhysics_[2];
  TTreeReaderValue<TTACPhysics> *TACPhysics_;
  TTreeReaderValue<TZDDPhysics> *ZDDPhysics_;
  TTreeReaderValue<TMust2Physics> *MUST2Physics_;
  TCATSPhysics *CATSPhysics;
  TTACPhysics TACPhysics;
  TZDDPhysics ZDDPhysics;
  TMust2Physics Must2Physics;
  TCutG *CUT[2]; 
  TFile *CFile[2]; 
  Double_t TopLeft[2][2];
  Double_t BotLeft[2][2];
  Double_t BotRight[2][2];

  Double_t TopLeftGeo[4] = {-12.5,12.7,-12.4,11.4};
  Double_t BotLeftGeo[4] = {-12.5,-12.3,-12.4,-13.6};
  Double_t BotRightGeo[4] = {12.5,-12.3,12.6,-13.6};

  std::map<int,std::map<string,TH1F*>>* TH1Map; 
  TGraph* Graph;

  double MaskPosX[2] = {0,0.1};
  double MaskPosY[2] = {0.2,-1.1};
  double MaskTargetX = -2.4;
  double MaskTargetY = -2.1;

  double tini;
  double tfinale;
  double pas;
  double lambda;
  double reponse;
  int n;
  
  unsigned short M2_TelescopeM;
  std::vector<double> M2_Ex_p;
  std::vector<double> M2_Ex_d;
  std::vector<double> M2_Ex_t;
  std::vector<double> M2_Ex_a;
  std::vector<double> M2_CsI_E_p;
  std::vector<double> M2_CsI_E_d;
  std::vector<double> M2_CsI_E_t;
  std::vector<double> M2_CsI_E_a;
  std::vector<double> M2_ExNoBeam;
  std::vector<double> M2_ExNoProton;
  std::vector<double> M2_EDC;
  std::vector<double> M2_ELab;
  std::vector<double> M2_ThetaLab;
  std::vector<double> M2_ThetaCM;
  std::vector<double> M2_X;
  std::vector<double> M2_Y;
  std::vector<double> M2_Z;
  std::vector<double> M2_dE;
  
  TTreeReaderValue<unsigned short>*M2_TelescopeM_;
  TTreeReaderValue<std::vector<double>>* M2_Ex_p_;
  TTreeReaderValue<std::vector<double>>* M2_Ex_d_;
  TTreeReaderValue<std::vector<double>>* M2_Ex_t_;
  TTreeReaderValue<std::vector<double>>* M2_Ex_a_;
  TTreeReaderValue<std::vector<double>>* M2_CsI_E_p_;
  TTreeReaderValue<std::vector<double>>* M2_CsI_E_d_;
  TTreeReaderValue<std::vector<double>>* M2_CsI_E_t_;
  TTreeReaderValue<std::vector<double>>* M2_CsI_E_a_;
  TTreeReaderValue<std::vector<double>>* M2_ExNoBeam_;
  TTreeReaderValue<std::vector<double>>* M2_ExNoProton_;
  TTreeReaderValue<std::vector<double>>* M2_EDC_;
  TTreeReaderValue<std::vector<double>>* M2_ELab_;
  TTreeReaderValue<std::vector<double>>* M2_ThetaLab_;
  TTreeReaderValue<std::vector<double>>* M2_ThetaCM_;
  TTreeReaderValue<std::vector<double>>* M2_X_;
  TTreeReaderValue<std::vector<double>>* M2_Y_;
  TTreeReaderValue<std::vector<double>>* M2_Z_;
  TTreeReaderValue<std::vector<double>>* M2_dE_;
  NPL::Reaction *Cr48_pd = new Reaction("48Cr(p,d)47Cr@1620MeV"); 
  NPL::Reaction *Cr48_pt = new Reaction("48Cr(p,t)46Cr@1620MeV"); 
  double TargetThickness = 53*micrometer;
  NPL::EnergyLoss deuteron_CH2 = NPL::EnergyLoss("deuteron_CH2.G4table","G4table",100); 
  //ClassDef(MCRMethod,0);
};












  


#endif
