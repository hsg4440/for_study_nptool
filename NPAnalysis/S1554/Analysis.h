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
#include "TSharcPhysics.h"
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
  TVector3 BeamDirection;
  TVector3 TargetPosition;


  double TargetThickness ;
  // Beam Energy
  double OriginalBeamEnergy ; // AMEV
                                                           // intermediate variable
  TRandom3 Rand ;
  int DetectorNumber  ;
  int RunNumber;
  int RunNumberMinor;
  double ThetaNormalTarget;
  double ThetaM2Surface ;
  double Si_E_M2 ;
  double CsI_E_M2  ;
  double Energy ;
  double E_M2 ;
  
  double ThetaSharcSurface ;
  double X_Sharc ;
  double Y_Sharc ;
  double Z_Sharc  ;
  double Si_E_Sharc ;
  double E_Sharc ;
  double Si_X_Sharc ;
  double Si_Y_Sharc ;
  TSharcPhysics* Sharc;
};
#endif
