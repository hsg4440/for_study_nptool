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
  InitOutputBranch();
  Rand = TRandom3();

  TargetThickness = 0.044*micrometer;

  Transfer = new NPL::Reaction("238U(12C,10Be)240Pu@1428");

  // Energy loss table
  Be10C = EnergyLoss("EnergyLossTable/Be10_C.G4table","G4Table",100);
  U238C = EnergyLoss("EnergyLossTable/U238_C.G4table","G4Table",100);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  XTarget = 0;
  YTarget = 0;
  ZTarget = 0;

  TVector3 BeamPosition(XTarget,YTarget,ZTarget);
  TVector3 PositionOnTarget(0,0,0);
  BeamEnergy = 1428.;
  //BeamEnergy = U238C.Slow(BeamEnergy,TargetThickness*0.5,0);
  Transfer->SetBeamEnergy(BeamEnergy);
  
  //cout << PISTA->EventMultiplicity << endl;
  /*if(PISTA->EventMultiplicity==1){
    for(unsigned int i = 0; i<PISTA->EventMultiplicity; i++){
      double Energy = PISTA->DE[i] + PISTA->E[i];
      DeltaE = PISTA->DE[i];
      Eres = PISTA->E[i];
      
      TVector3 HitDirection = PISTA->GetPositionOfInteraction(i)-PositionOnTarget;
      Xcalc = PISTA->GetPositionOfInteraction(i).X();
      Ycalc = PISTA->GetPositionOfInteraction(i).Y();
      Zcalc = PISTA->GetPositionOfInteraction(i).Z();
      //ThetaLab = HitDirection.Angle(BeamDirection);
      
      ThetaLab = HitDirection.Angle(TVector3(0,0,1));
      ThetaDetectorSurface = HitDirection.Angle(-PISTA->GetDetectorNormal(i));
      
      DeltaE = DeltaE/cos(ThetaDetectorSurface);
      //PID = pow(Energy,1.78)-pow(PISTA->E[i],1.78);
      PID = pow(DeltaE+Eres,1.78)-pow(Eres,1.78);

      ThetaNormalTarget = HitDirection.Angle(TVector3(0,0,1));
      Elab = Be10C.EvaluateInitialEnergy(Energy,TargetThickness*0.5,ThetaNormalTarget);
      Ex = Transfer->ReconstructRelativistic(Elab, ThetaLab);
      ThetaCM = Transfer->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
      ThetaLab = ThetaLab/deg;
    }
  }*/
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("XTarget",&XTarget,"XTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("YTarget",&YTarget,"YTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZTarget",&ZTarget,"ZTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaE",&DeltaE,"DeltaE/D");
  RootOutput::getInstance()->GetTree()->Branch("Eres",&Eres,"Eres/D");
  RootOutput::getInstance()->GetTree()->Branch("PID",&PID,"PID/D");
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab,"Elab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Xcalc",&Xcalc,"Xcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Ycalc",&Ycalc,"Ycalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Zcalc",&Zcalc,"Zcalc/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  BeamEnergy = -1000;
  Ex = -1000;
  DeltaE = -1000;
  Eres = -1000;
  Elab = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
  XTarget = -1000;
  YTarget = -1000;
  ZTarget = -1000;
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

