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

#include "TMust2Physics.h"
#include"NPVAnalysis.h"
#include"TZDDPhysics.h"
#include"NPEnergyLoss.h"
#include"NPFunction.h"
#include"NPReaction.h"
#include"NPOptionManager.h"
#include"RootInput.h"
#include"TInitialConditions.h"
#include "TReactionConditions.h"
#include"NPParticle.h"
#include"NPBeam.h"
class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void SetParticles();
    void SetEnergyLoss();
    void Init();
    void TreatEvent();
    void End();
    void ReInit();
    void InitOutputBranch();
    void InitInputBranch();
    bool CheckGoodEvent();
    // CheckIC is checking that all IC are crossed (E>0) not less and not more than 1 time
    bool CheckIC(); 
    // CheckPlastics is checking that only and at leat 1 plastic is hit
    bool CheckPlastics();
    // CHecking multiplicity 1 in DC
    bool CheckDC(); 

   static NPL::VAnalysis* Construct();

  private:

  /// Currently only treating multiplicity 1 events
  // ZDD info 
  double ZDD_DC_X;
  double ZDD_DC_Y;
  double ZDD_ThetaIC = 30/deg;
  double ZDD_ThetaAfterTarget;
  double ZDD_ThetaAfterTarget_X;
  double ZDD_E_tot;
  double ZDD_E_Plastic;
  double ZDD_dE_tot;
  
  // MUST2 info 
  double M2_Ex;
  double M2_ExNoBeam;
  double M2_ExNoProton;
  double M2_EDC;
  double M2_ELab;
  double M2_ThetaLab;
  double M2_ThetaCM;
  double M2_X;
  double M2_Y;
  double M2_Z;
  double M2_dE;

  double OriginalBeamEnergy ; // AMEV
  double FinalBeamEnergy; 


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
  NPL::EnergyLoss Beam_Target;
  NPL::EnergyLoss Heavy_Target;
  NPL::EnergyLoss LightTarget;
  std::vector<NPL::EnergyLoss> Heavy_IC_Gas;
  std::vector<NPL::EnergyLoss> Heavy_IC_Windows;
  std::vector<NPL::EnergyLoss> Heavy_IC_Mylar;
  std::vector<NPL::EnergyLoss> Heavy_DC_Gas;
  NPL::EnergyLoss LightAl;
  NPL::EnergyLoss LightSi;

  string BeamName = "48Cr";
  
  int Nb_Mylar_After_IC = 1;
  int Nb_Mylar_Before_IC = 4;
  double Drift_Chamber_Length = 400*mm;
  double Drift_Chamber_Width = 400*mm;
  double ZDD_R = 1*m;

  private:
  TMust2Physics* M2;
  TZDDPhysics* ZDD;
  TInitialConditions* InitialConditions;
  TReactionConditions* ReactionConditions;
  
  Beam* BeamPart = new Beam;
  Particle* HeavyEjectile = new Particle;
};
#endif