#include "NPReaction.h"
#include "TChain.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TMust2Physics.h"
#include "TExogamPhysics.h"
#include "TZDDPhysics.h"
#include "TTACPhysics.h"
#include "TCATSPhysics.h"
  
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
  
  unsigned short IC_ZDDRaw[6];
  TTreeReaderArray<UShort_t>* IC_ZDDRaw_;
  unsigned long long IC_ZDDRawTS[6];
  TTreeReaderArray<ULong64_t>* IC_ZDDRaw_TS_;
  float ICCal[4];
  
  unsigned short TAC_CATS_PL;
  TTreeReaderValue<UShort_t>* TAC_CATS_PL_;
  unsigned long long TAC_CATS_PLTS;
  TTreeReaderValue<ULong64_t>* TAC_CATS_PL_TS_;
  
  unsigned short TAC_CATS_HF;
  TTreeReaderValue<UShort_t>* TAC_CATS_HF_;
  unsigned long long TAC_CATS_HFTS;
  TTreeReaderValue<ULong64_t>* TAC_CATS_HF_TS_;
  
  unsigned short TAC_CATS_EXOGAM;
  TTreeReaderValue<UShort_t>* TAC_CATS_EXOGAM_;
  unsigned long long TAC_CATS_EXOGAMTS;
  TTreeReaderValue<ULong64_t>* TAC_CATS_EXOGAM_TS_;
  
  unsigned short TAC_MMG_CATS2;
  TTreeReaderValue<UShort_t>* TAC_MMG_CATS2_;
  unsigned long long TAC_MMG_CATS2TS;
  TTreeReaderValue<ULong64_t>* TAC_MMG_CATS2_TS_;
  
  unsigned short TAC_MMG_CATS1;
  TTreeReaderValue<UShort_t>* TAC_MMG_CATS1_;
  unsigned long long TAC_MMG_CATS1TS;
  TTreeReaderValue<ULong64_t>* TAC_MMG_CATS1_TS_;
  
  unsigned short TAC_MMG_EXOGAM;
  TTreeReaderValue<UShort_t>* TAC_MMG_EXOGAM_;
  unsigned long long TAC_MMG_EXOGAMTS;
  TTreeReaderValue<ULong64_t>* TAC_MMG_EXOGAM_TS_;
  
  unsigned short TAC_CATS1_CATS2;
  TTreeReaderValue<UShort_t>* TAC_CATS1_CATS2_;
  unsigned long long TAC_CATS1_CATS2TS;
  TTreeReaderValue<ULong64_t>* TAC_CATS1_CATS2_TS_;
  
  unsigned short TAC_D4_CATS1;
  TTreeReaderValue<UShort_t>* TAC_D4_CATS1_;
  unsigned long long TAC_D4_CATS1TS;
  TTreeReaderValue<ULong64_t>* TAC_D4_CATS1_TS_;
  
  unsigned short TAC_PL_1;
  TTreeReaderValue<UShort_t>* TAC_PL_1_;
  unsigned long long TAC_PL_1TS;
  TTreeReaderValue<ULong64_t>* TAC_PL_1_TS_;
  unsigned short TAC_PL_2;
  TTreeReaderValue<UShort_t>* TAC_PL_2_;
  unsigned long long TAC_PL_2TS;
  TTreeReaderValue<ULong64_t>* TAC_PL_2_TS_;
  unsigned short TAC_PL_3;
  TTreeReaderValue<UShort_t>* TAC_PL_3_;
  unsigned long long TAC_PL_3TS;
  TTreeReaderValue<ULong64_t>* TAC_PL_3_TS_;
  unsigned short TAC_PL_4;
  TTreeReaderValue<UShort_t>* TAC_PL_4_;
  unsigned long long TAC_PL_4TS;
  TTreeReaderValue<ULong64_t>* TAC_PL_4_TS_;
  unsigned short TAC_PL_5;
  TTreeReaderValue<UShort_t>* TAC_PL_5_;
  unsigned long long TAC_PL_5TS;
  TTreeReaderValue<ULong64_t>* TAC_PL_5_TS_;
  
  std::vector<unsigned int> GATCONFMASTER;
  TTreeReaderValue<std::vector<unsigned int>>* GATCONFMASTER_;
  
  unsigned short DATATRIG_CATS;
  TTreeReaderValue<unsigned short>* DATATRIG_CATS_;
  
  unsigned short PlasticRaw[10];
  TTreeReaderArray<UShort_t>* PlasticRaw_;
  unsigned long long PlasticRawTS[10];
  TTreeReaderArray<ULong64_t>* PlasticRaw_TS_;
    
  TChain* c = new TChain("PhysicsTree");
  TTreeReader* TreeReader;
  
  TMust2Physics Must2Physics;
  TExogamPhysics ExogamPhysics;
  TCATSPhysics CATSPhysics;
  TZDDPhysics ZDDPhysics;
  TTACPhysics TACPhysics;
  TTreeReaderValue<TCATSPhysics> *CATSPhysics_;
  TTreeReaderValue<TExogamPhysics> *ExogamPhysics_;
  TTreeReaderValue<TTACPhysics> *TACPhysics_;
  TTreeReaderValue<TZDDPhysics> *ZDDPhysics_;
  TTreeReaderValue<TMust2Physics> *Must2Physics_;
  
  
  TFile * f_cut_deuton = new TFile("./CUT_deuton.root");
  TCutG* cut_deuton =  (TCutG*) f_cut_deuton->FindObjectAny("CUT_deuton");
  TFile * f_cut_triton = new TFile("./CUT_triton.root");
  TCutG* cut_triton =  (TCutG*) f_cut_triton->FindObjectAny("CUT_triton");

  TFile * f_cut_Cr = new TFile("./CUT_Cr.root");
  TCutG* cut_Cr =  (TCutG*) f_cut_Cr->FindObjectAny("CUTCr");
  
  
  //NPL::Reaction Cr48_pd("48Cr(p,d)47Cr@1511MeV"); 
  //NPL::Reaction Cr48_pt("48Cr(p,t)46Cr@1511MeV"); 
  NPL::Reaction Cr48_pd("48Cr(p,d)47Cr@1620MeV"); 
  NPL::Reaction Cr48_pt("48Cr(p,t)46Cr@1620MeV"); 
  double TargetThickness = 53*micrometer;
  NPL::EnergyLoss deuteron_CH2 = NPL::EnergyLoss("deuteron_CH2.G4table","G4table",100); 
  NPL::EnergyLoss triton_CH2 = NPL::EnergyLoss("triton_CH2.G4table","G4table",100); 