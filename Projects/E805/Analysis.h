#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  MUST_AND_ZDD analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include"Geometry_Clover_Exogam.h"
#include "TMust2Physics.h"
#include "TCATSPhysics.h"
#include "TMust2PhysicsReader.h"
#include"NPVAnalysis.h"
#include"TZDDPhysics.h"
#include"NPEnergyLoss.h"
#include"NPFunction.h"
#include"NPReaction.h"
#include"NPOptionManager.h"
#include"RootInput.h"
#include"RootOutput.h"
#include"TInitialConditions.h"
#include "TReactionConditions.h"
#include"NPParticle.h"
#include"NPBeam.h"
#include "TCutG.h"
#include<random>
class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void SetParticles();
    void SetEnergyLoss();
    void Init();
    bool UnallocateBeforeBuild();
    bool FillOutputCondition();
    bool UnallocateBeforeTreat(){return true;};
    void TreatEvent();
    void End();
    void ReInit();
    void InitOutputBranch();
    void InitInputBranch();
    bool CheckGoodEvent();
    void UnallocateVariables();
    // CheckIC is checking that all IC are crossed (E>0) not less and not more than 1 time
    bool CheckIC(); 
    // CheckPlastics is checking that only and at leat 1 plastic is hit
    bool CheckPlastics();
    // CHecking multiplicity 1 in DC
    bool CheckDC();
    bool CheckExoAddBack(int ExoDet1, int ExoCristal1, int ExoSegment1, float ExoTime1, int ExoDet2, int ExoCristal2, int ExoSegment2, float ExoTime2);
    bool CheckExoDeltaTV(float ExoTime);

    static NPL::VAnalysis* Construct();
    TTreeReaderValue<TMust2Data>* r_ReaderEventData; 

  private:

  /// Currently only treating multiplicity 1 events
  // ZDD info 
  double ZDD_DC_X;
  double ZDD_DC_Y;
  double ZDD_ThetaIC = 30/deg;
  double ZDD_ThetaAfterTarget;
  double ZDD_ThetaAfterTarget_X;
  double ZDD_ThetaLab;
  double ZDD_E_tot;
  double ZDD_E_Plastic;
  double ZDD_dE_tot;
  std::vector<double> ZDD_Corrected_IC_E;
  
  // MUST2 info 
  unsigned short M2_TelescopeM;
  std::vector<double> M2_Ex;
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
  std::vector<double> M2_ECsI_from_deltaE ;
  std::vector<double> Beta_from_deltaE;
  std::vector<double> Beth_from_deltaE;

  double OriginalBeamEnergy ; // AMEV
  double FinalBeamEnergy; 


  float CATS1_X; 
  float CATS2_X; 
  float CATS1_Y; 
  float CATS2_Y; 
  
  float Xf; 
  float Yf; 
  float Tf; 
  float Pf; 

  unsigned long long MUGAST_TS[1];
  unsigned long long DATATRIG_CATSTS[1];
 
  int Inner6MVM;
  TTreeReaderValue<int>* Inner6MVM_;
  float Inner6MV[48];
  TTreeReaderArray<float>* Inner6MV_;
  unsigned short Inner6MVN[48];
  TTreeReaderArray<unsigned short>* Inner6MVN_;
  unsigned long long Inner6MVTS[48];
  
  
  int BGOVM;
  TTreeReaderValue<int>* BGOVM_;
  float BGOV[48];
  TTreeReaderArray<float>* BGOV_;
  unsigned short BGOVN[48];
  TTreeReaderArray<unsigned short>* BGOVN_;
  
  unsigned short GATCONFMASTER;
  TTreeReaderValue<unsigned short>* GATCONFMASTER_;
  
  unsigned short DATATRIG_CATS;
  TTreeReaderValue<unsigned short>* DATATRIG_CATS_;

//  int Strip_X_M;
//  float Strip_X_E[48];
//  float Strip_X_T[48];
//  unsigned int Strip_X_Nb[48];
//  unsigned short Strip_X_Det[48];
//  
//  int Strip_Y_M;
//  float Strip_Y_E[48];
//  float Strip_Y_T[48];
//  unsigned int Strip_Y_Nb[48];
//  unsigned short Strip_Y_Det[48];
  
  int DeltaTVM;
  float DeltaTV[48];
  unsigned short DeltaTVN[48];
  unsigned long long DeltaTVTS[48];
  
  float EnergyDoppler; 
  float EnergyAddBackDoppler; 
  float EnergyAddBack;
  int ExogamDetNb[3];
  int CristalNb[3];
  int SegmentNb[3];
  
  std::vector<int> event1;
  std::vector<int> event2;
  int highest_E;

  int OutersVM;
  TTreeReaderValue<int>* OutersVM_;
  float OutersV[192];
  TTreeReaderArray<float>* OutersV_;
  unsigned short OutersVN[192];
  TTreeReaderArray<unsigned short>* OutersVN_;

  int DCRawM;
  unsigned short DCRaw[4];
  unsigned short DCRawNr[4];
  unsigned long long DCRawTS[4];
  
  unsigned short PlasticRaw[10];
  unsigned long long PlasticRawTS[10];
  float PlasticCal[10];
  
  float PlasticEner[5];
  int PlasticEnerM;
  unsigned short PlasticEnerN[5];
  unsigned long long PlasticEnerTS[5];
  int PlasticCounter;
  float PlasticThreshold;
  float PlasticEner_tmp;
  
  unsigned short IC_ZDDRaw[6];
  unsigned long long IC_ZDDRawTS[6];
  float ICCal[4];
  
  unsigned short TAC_CATS_PL;
  unsigned long long TAC_CATS_PLTS;
  
  unsigned short TAC_CATS_HF;
  unsigned long long TAC_CATS_HFTS;
  
  unsigned short TAC_CATS_EXOGAM;
  unsigned long long TAC_CATS_EXOGAMTS;
  
  unsigned short TAC_MMG_CATS2;
  unsigned long long TAC_MMG_CATS2TS;
  
  unsigned short TAC_MMG_CATS1;
  unsigned long long TAC_MMG_CATS1TS;
  
  unsigned short TAC_MMG_EXOGAM;
  unsigned long long TAC_MMG_EXOGAMTS;
  
  unsigned short TAC_CATS1_CATS2;
  unsigned long long TAC_CATS1_CATS2TS;
  
  unsigned short TAC_D4_CATS1;
  unsigned long long TAC_D4_CATS1TS;
  
  unsigned short TAC_PL_1;
  unsigned long long TAC_PL_1TS;
  unsigned short TAC_PL_2;
  unsigned long long TAC_PL_2TS;
  unsigned short TAC_PL_3;
  unsigned long long TAC_PL_3TS;
  unsigned short TAC_PL_4;
  unsigned long long TAC_PL_4TS;
  unsigned short TAC_PL_5;
  unsigned long long TAC_PL_5TS;
  
  double xtarget;
  double ytarget;
  double IncidentTheta = 0;
  int DetectorNumber  ;
  double ThetaNormalTarget;
  double ThetaM2Surface ;
  double ThetaMGSurface ;
  double Si_E_M2 ;
  double CsI_E_M2  ;
  double Energy ;
  double BeamEnergy;
  double ThetaGDSurface ;
  
  double Beta;
  double Beta_light;
  double Gamma;
  double Velocity;
  
  double Drift_Speed;
  double TargetThickness;
  double CorrectedBeamEnergy;
  std::vector<double> IC_Energy;
  
  TVector3 ZDir;
  TVector3 BeamDirection;
  TVector3 BeamImpact;


  
  NPL::Reaction* reaction = new Reaction;
  NPL::Particle* Co_57 = new Particle;
  NPL::Particle* Ni_57 = new Particle;
  NPL::EnergyLoss Beam_Target;
  NPL::EnergyLoss Heavy_Target;
  NPL::EnergyLoss LightTarget;
  NPL::EnergyLoss ProtonSi;
  std::vector<NPL::EnergyLoss> Heavy_IC_Gas;
  std::vector<NPL::EnergyLoss> Heavy_IC_Windows;
  std::vector<NPL::EnergyLoss> Heavy_IC_Mylar;
  std::vector<NPL::EnergyLoss> Heavy_DC_Gas;
  NPL::EnergyLoss LightAl;
  NPL::EnergyLoss LightSi;

  string BeamName = "48Cr";
  
  int Nb_Mylar_After_IC = 1;
  int Nb_Mylar_Before_IC = 4;
  int Nb_IC = 10;
  double Drift_Chamber_Length = 180*mm;
  double Drift_Chamber_Width = 180*mm;
  double ZDD_R = 1027*mm;

  private:
  TMust2Data* M2_Raw;
  TMust2Physics* M2;
  TCATSPhysics* CATS;
  TZDDPhysics* ZDD;
  TMust2PhysicsReader* M2_Reader;
  TInitialConditions* InitialConditions;
  TReactionConditions* ReactionConditions;
  
  std::ifstream ExogamTopo;
  std::ifstream DC_calib;
  std::ifstream Plastic_calib;
  std::ifstream IC_calib;
  std::string ExogamLine;
  std::string DCLine;
  std::string PlasticLine;
  std::string ICLine;

  Int_t Plastic_Nb_tmp;
  Double_t Plastic_peak_tmp,Plastic_pedestal_tmp;
  Double_t Plastic_pedestal[10];
  Double_t Plastic_peak[10];
  
  Int_t IC_Nb_tmp;
  Double_t IC_peak_tmp,IC_pedestal_tmp;
  Double_t IC_pedestal[4];
  Double_t IC_peak[4];
  
  int DC_Nr;
  long long dt;
  char DC_XY[1];
  Int_t DC_numb_tmp;
  Double_t off_tmp,cff_tmp,sqr_tmp;
  Double_t off[4];
  Double_t cff[4];
  Double_t sqr[4];

  double DC_X;
  double DC_Y;

  //Char_t* ExogamLine_char[100];


  double Theta_seg;
  double Phi_seg;

  Char_t Exogam[100];
  Int_t ExoNumb;
  Int_t Flange_tmp;
  int FlangeNumb[12];



  /////////////Exogam
  Clover_struc Exogam_Clovers_struc[12];


  TCutG* proton_cut[4];
  TCutG* deuteron_cut[4];
  TCutG* triton_cut[4];
  
  TFile* proton_cut_file[4];
  TFile* deuteron_cut_file[4];
  TFile* triton_cut_file[4];
  
  Beam* BeamPart = new Beam;
  Particle* HeavyEjectile = new Particle;
  
  default_random_engine generator;
  normal_distribution<double> distribution = normal_distribution<double>(0.0,1.0);

};
#endif