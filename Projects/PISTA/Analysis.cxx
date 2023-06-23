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
 *  This class describe  PISTA analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"RootOutput.h"
#include"RootInput.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  PISTA= (TPISTAPhysics*) m_DetectorManager->GetDetector("PISTA");
  InitialConditions = new TInitialConditions();
  ReactionConditions = new TReactionConditions();
  InteractionCoordinates = new TInteractionCoordinates();
  InitOutputBranch();
  InitInputBranch();
  Rand = TRandom3();

  TargetThickness = m_DetectorManager->GetTargetThickness();

  Transfer = new NPL::Reaction("238U(12C,12C)238U@1417");
  //Transfer = new NPL::Reaction("238U(12C,14C)236U@1428");

  // Energy loss table
  Be10C = EnergyLoss("EnergyLossTable/Be10_C.G4table","G4Table",100);
  C12C = EnergyLoss("EnergyLossTable/C12_C.G4table","G4Table",100);
  U238C = EnergyLoss("EnergyLossTable/U238_C.G4table","G4Table",100);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  OriginalThetaLab = ReactionConditions->GetTheta(0);
  OriginalElab = ReactionConditions->GetKineticEnergy(0);
  OriginalBeamEnergy = ReactionConditions->GetBeamEnergy();
  OriginalEx = ReactionConditions->GetExcitation4();

  int mult = InteractionCoordinates->GetDetectedMultiplicity();
  if(mult==1){
    for(int i=0; i<mult; i++){
      Xpista = InteractionCoordinates->GetDetectedPositionX(i);
      Ypista = InteractionCoordinates->GetDetectedPositionY(i);
      Zpista = InteractionCoordinates->GetDetectedPositionZ(i);
      R = sqrt(Xpista*Xpista + Ypista*Ypista + Zpista*Zpista);
    }
  }
  XTarget = InitialConditions->GetIncidentPositionX();
  YTarget = InitialConditions->GetIncidentPositionY();
  ZTarget = InitialConditions->GetIncidentPositionZ();

  TVector3 BeamDirection = InitialConditions->GetBeamDirection();
  TVector3 BeamPosition(XTarget,YTarget,ZTarget);

  //TVector3 PositionOnTarget(0,0,0);
  TVector3 PositionOnTarget(Rand.Gaus(XTarget, 0.6/2.35), Rand.Gaus(YTarget, 0.6/2.35), 0);
  //TVector3 PositionOnTarget(XTarget, YTarget, 0);
  BeamEnergy = 1417.;//InitialConditions->GetIncidentInitialKineticEnergy();
  BeamEnergy = U238C.Slow(BeamEnergy,TargetThickness*0.5,0);
  Transfer->SetBeamEnergy(BeamEnergy);
  //Transfer->SetBeamEnergy(OriginalBeamEnergy);
  
  //cout << PISTA->EventMultiplicity << endl;
  if(PISTA->EventMultiplicity==1){
    for(unsigned int i = 0; i<PISTA->EventMultiplicity; i++){
      double Energy = PISTA->DE[i] + PISTA->E[i];
      DeltaE = PISTA->DE[i];
      Eres = PISTA->E[i];
      int det = PISTA->DetectorNumber[i];
      int strip_E = PISTA->E_StripNbr[i];
      int strip_DE = PISTA->DE_StripNbr[i];
      
      TVector3 PISTA_pos = PISTA->GetPositionOfInteraction(det,strip_E,strip_DE);
      TVector3 HitDirection = PISTA_pos;// - PositionOnTarget;
      Xcalc = PISTA_pos.X();
      Ycalc = PISTA_pos.Y();
      Zcalc = PISTA_pos.Z();
      //ThetaLab = HitDirection.Angle(BeamDirection);
      
      ThetaLab = HitDirection.Angle(TVector3(0,0,1));
      ThetaDetectorSurface = HitDirection.Angle(-PISTA->GetDetectorNormal(i));
      
      DeltaE = DeltaE/cos(ThetaDetectorSurface);
      //PID = pow(Energy,1.78)-pow(PISTA->E[i],1.78);
      PID = pow(DeltaE+Eres,1.78)-pow(Eres,1.78);

      ThetaNormalTarget = HitDirection.Angle(TVector3(0,0,1));
      //Elab = Be10C.EvaluateInitialEnergy(Energy,TargetThickness*0.5,ThetaNormalTarget);
      Elab = C12C.EvaluateInitialEnergy(Energy,TargetThickness*0.5,ThetaNormalTarget);
      OptimumEx = Transfer->ReconstructRelativistic(OriginalElab, OriginalThetaLab*deg);
      Ex = Transfer->ReconstructRelativistic(Elab, ThetaLab);
      ThetaCM = Transfer->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
      ThetaLab = ThetaLab/deg;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("OriginalBeamEnergy",&OriginalBeamEnergy,"OriginalBeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("OriginalEx",&OriginalEx,"OriginalEx/D");
  RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("XTarget",&XTarget,"XTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("YTarget",&YTarget,"YTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZTarget",&ZTarget,"ZTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("OptimumEx",&OptimumEx,"OptimumEx/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaE",&DeltaE,"DeltaE/D");
  RootOutput::getInstance()->GetTree()->Branch("Eres",&Eres,"Eres/D");
  RootOutput::getInstance()->GetTree()->Branch("PID",&PID,"PID/D");
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab,"Elab/D");
  RootOutput::getInstance()->GetTree()->Branch("OriginalElab",&OriginalElab,"OriginalElab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("OriginalThetaLab",&OriginalThetaLab,"OriginalThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("R",&R,"R/D");
  RootOutput::getInstance()->GetTree()->Branch("Xpista",&Xpista,"Xpista/D");
  RootOutput::getInstance()->GetTree()->Branch("Ypista",&Ypista,"Ypista/D");
  RootOutput::getInstance()->GetTree()->Branch("Zpista",&Zpista,"Zpista/D");
  RootOutput::getInstance()->GetTree()->Branch("Xcalc",&Xcalc,"Xcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Ycalc",&Ycalc,"Ycalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Zcalc",&Zcalc,"Zcalc/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);
  RootInput::getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fRC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&ReactionConditions);
  RootInput::getInstance()->GetChain()->SetBranchStatus("InteractionCoordinates",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fDetected_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",&InteractionCoordinates);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  OriginalBeamEnergy = -1000;
  OriginalEx = -1000;
  BeamEnergy = -1000;
  OptimumEx = -1000;
  Ex = -1000;
  DeltaE = -1000;
  Eres = -1000;
  Elab = -1000;
  OriginalElab = -1000;
  OriginalThetaLab = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
  XTarget = -1000;
  YTarget = -1000;
  ZTarget = -1000;
  OriginalThetaLab = -1000;
  R = -1000;
  Xpista = -1000;
  Ypista = -1000;
  Zpista = -1000;
  Xcalc = -1000;
  Ycalc = -1000;
  Zcalc = -1000;
  PID = -1000;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy{
    public:
      proxy(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy p;
}

