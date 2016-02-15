#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2025                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include"NPVAnalysis.h"
#include"NPEnergyLoss.h"
#include"NPReaction.h"
#include"RootOutput.h"
#include"RootInput.h"
#include "TMust2Physics.h"
#include "GaspardTracker.h"
#include "TInitialConditions.h"
#include <TRandom3.h>
#include <TVector3.h>
#include <TMath.h>

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();

  void InitOutputBranch();
  void InitInputBranch();
  void ReInitValue();
  static NPL::VAnalysis* Construct();
 
  private:
  double Ex;
  double ELab;
  double ThetaLab;
  double ThetaCM;
  NPL::Reaction* myReaction;
  TInitialConditions* myInit ;
  //	Energy loss table: the G4Table are generated by the simulation
  EnergyLoss LightCD2;
  EnergyLoss LightAl;
  EnergyLoss LightSi;
  EnergyLoss BeamCD2;
  EnergyLoss* BeamWindow;
  EnergyLoss* LightWindow;

  TVector3 BeamImpact;
  
  double TargetThickness ;
  double WindowsThickness;
  // Beam Energy
  double OriginalBeamEnergy ; // AMEV
  double FinalBeamEnergy; 
  
  // intermediate variable
  TRandom3 Rand ;
  int Run;
  int DetectorNumber  ;
  double ThetaNormalTarget;
  double ThetaM2Surface ;
  double Si_E_M2 ;
  double CsI_E_M2  ;
  double Energy ;
  double E_M2 ;
  
  double ThetaGDSurface ;
  double X_GD ;
  double Y_GD ;
  double Z_GD  ;
  double Si_E_GD ;
  double E_GD ;
  double Si_X_GD ;
  double Si_Y_GD ;
  TMust2Physics* M2;
  GaspardTracker* GD;
};
#endif
