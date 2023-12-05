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
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "TMust2Physics.h"
#include "TMugastPhysics.h"
#include "TCATSPhysics.h"
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
  double ExNoBeam;
  double ExNoProton;
  double EDC;
  double ELab;
  double ThetaLab;
  double ThetaCM;
  double OriginalELab;
  double OriginalThetaLab;

  NPL::Reaction reaction;
    //	Energy loss table: the G4Table are generated by the simulation
  NPL::EnergyLoss LightTarget;
  NPL::EnergyLoss LightAl;
  NPL::EnergyLoss LightSi;
  NPL::EnergyLoss BeamTarget;
  NPL::EnergyLoss* BeamWindow;
  NPL::EnergyLoss* LightWindow;

  double TargetThickness ;
  double WindowsThickness;
  // Beam Energy
  double OriginalBeamEnergy ; // AMEV
  double FinalBeamEnergy; 
  
  // intermediate variable
  TVector3 BeamDirection;
  TVector3 BeamImpact;
  TRandom3 Rand ;
  int Run;
  int DetectorNumber  ;
  double ThetaNormalTarget;
  double ThetaM2Surface ;
  double ThetaMGSurface ;
  double Si_E_M2 ;
  double CsI_E_M2  ;
  double Energy ;
  double BeamEnergy;
  
  double ThetaGDSurface ;
  double X ;
  double Y ;
  double Z ;
  // Vamos Branches
  unsigned long long int LTS;
  
  // Agata branches
   double agata_zShift;
   unsigned long long int TStrack;
   int nbHits;
   int nbTrack;
   float *trackE= new float(100);
   float *trackX1= new float(100);
   float *trackY1= new float(100);
   float *trackZ1= new float(100);
   float *trackT= new float(100);
   int *trackCrystalID = new int(100);
   int nbCores;
   int *coreId= new int(100);
   ULong64_t *coreTS= new ULong64_t(100);
   float *coreE0= new float(100);
   //
  double dE;
  double dTheta;
  // Branches and detectors
  TMust2Physics* M2;
  TMugastPhysics* MG;
  TCATSPhysics* CATS;
  bool simulation;
  TInitialConditions* Initial;
  TReactionConditions* ReactionConditions;
};
#endif
