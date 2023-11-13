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
#include "TCATSPhysics.h"



  
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
  TString path = "./RootR/";
  TTreeReaderValue<TCATSPhysics> *CATSPhysics_[2];
  TCATSPhysics *CATSPhysics;
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

  double tini;
  double tfinale;
  double pas;
  double lambda;
  double reponse;
  int n;
  //ClassDef(MCRMethod,0);
};












  


#endif
