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
  InitInputBranch();
  InitOutputBranch();
  Rand = TRandom3();

  TargetThickness = 0.44*micrometer;

  Transfer10Be = new NPL::Reaction("238U(12C,10Be)240Pu@1417");
  Transfer14C  = new NPL::Reaction("238U(12C,14C)236U@1417");
  Elastic      = new NPL::Reaction("238U(12C,12C)238U@1417");
  
  chain = RootInput::getInstance()->GetChain();

  // Energy loss table
  C12C = EnergyLoss("EnergyLossTable/C14_C.G4table","G4Table",100);
  Be10C = EnergyLoss("EnergyLossTable/Be10_C.G4table","G4Table",100);
  U238C = EnergyLoss("EnergyLossTable/U238_C.G4table","G4Table",100);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  double Xoffset = -1.0;
  double Yoffset = 1.0;
  if(fTP_X>-1500 && fTP_Y>-1500){
    XTarget = 0.5*fTP_X + Xoffset;
    YTarget = fTP_Y + Yoffset;
    ZTarget = 0;
  }
  else{  
    XTarget = -1.0;
    YTarget = -2 + Yoffset;
    ZTarget = 0;
  }


  TVector3 BeamPosition(XTarget,YTarget,ZTarget);
  TVector3 PositionOnTarget(XTarget,YTarget,ZTarget);
  BeamEnergy = 1417.;
  //BeamEnergy = U238C.Slow(BeamEnergy,TargetThickness*0.5,0);
  Transfer10Be->SetBeamEnergy(BeamEnergy);
  Transfer14C->SetBeamEnergy(BeamEnergy);
  Elastic->SetBeamEnergy(BeamEnergy);


  if(fIC[1]>0 && fIC[5]>0){
    Chio_DE = 0.5*(fIC[1]+fIC[2]+fIC[3]+fIC[4]);
    Chio_E = fIC[5]+fIC[6]+fIC[7];
  }

  //cout << PISTA->EventMultiplicity << endl;
  if(PISTA->EventMultiplicity==1){
    for(unsigned int i = 0; i<PISTA->EventMultiplicity; i++){
      double Energy = PISTA->DE[i] + PISTA->back_E[i];
      DeltaE = PISTA->DE[i];
      Eres = PISTA->back_E[i];

      int strip_DE = PISTA->DE_StripNbr[i];
      int strip_E = PISTA->E_StripNbr[i];
      if(strip_DE>0 && strip_DE<92 && strip_E>0 && strip_E<58){
        TVector3 HitDirection = PISTA->GetPositionOfInteraction(i)-PositionOnTarget;
        PhiLab = PISTA->GetPositionOfInteraction(i).Phi();
        Xcalc = PISTA->GetPositionOfInteraction(i).X();
        Ycalc = PISTA->GetPositionOfInteraction(i).Y();
        Zcalc = PISTA->GetPositionOfInteraction(i).Z();
        //ThetaLab = HitDirection.Angle(BeamDirection);

        ThetaLab = HitDirection.Angle(TVector3(0,0,1));
        ThetaDetectorSurface = HitDirection.Angle(PISTA->GetDetectorNormal(i));

        DeltaEcorr = DeltaE*cos(ThetaDetectorSurface);
        //PID = pow(Energy,1.78)-pow(PISTA->E[i],1.78);
        PID = pow(DeltaEcorr+Eres,1.78)-pow(Eres,1.78);

        ThetaNormalTarget = HitDirection.Angle(TVector3(0,0,1));
        Elab = Energy;//Be10C.EvaluateInitialEnergy(Energy,TargetThickness*0.5,ThetaNormalTarget);
        Ex240Pu = Transfer10Be->ReconstructRelativistic(Elab, ThetaLab);
        Ex236U  = Transfer14C->ReconstructRelativistic(Elab, ThetaLab);
        Ex238U  = Elastic->ReconstructRelativistic(Elab, ThetaLab);
        ThetaCM = Transfer10Be->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
        ThetaLab = ThetaLab/deg;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("XTarget",&XTarget,"XTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("YTarget",&YTarget,"YTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZTarget",&ZTarget,"ZTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex240Pu",&Ex240Pu,"Ex240Pu/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex236U",&Ex236U,"Ex236U/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex238U",&Ex238U,"Ex238U/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaE",&DeltaE,"DeltaE/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaEcorr",&DeltaEcorr,"DeltaEcorr/D");
  RootOutput::getInstance()->GetTree()->Branch("Eres",&Eres,"Eres/D");
  RootOutput::getInstance()->GetTree()->Branch("PID",&PID,"PID/D");
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab,"Elab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("PhiLab",&PhiLab,"PhiLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Xcalc",&Xcalc,"Xcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Ycalc",&Ycalc,"Ycalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Zcalc",&Zcalc,"Zcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Chio_DE",&Chio_DE,"Chio_DE/D");
  RootOutput::getInstance()->GetTree()->Branch("Chio_E",&Chio_E,"Chio_E/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("IC",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("IC",fIC);

  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_X",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_X",&fTP_X);
  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_Y",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_Y",&fTP_Y);

  RootInput::getInstance()->GetChain()->SetBranchStatus("TMWPat_00TS",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TMWPat_00TS",&fTS_TMW);
  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_PISTA_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_PISTA_C",&fTAC_MW1_PISTA);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW0_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW0_C",&fTAC_TMW1_FPMW0);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW1_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW1_C",&fTAC_TMW1_FPMW1);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW2_FPMW0_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW2_FPMW0_C",&fTAC_TMW2_FPMW0);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW2_FPMW1_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW2_FPMW1_C",&fTAC_TMW2_FPMW1);




  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner6MVM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InnerM6VM",&Inner6MVM);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner6MV",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner6MV",Inner6MV);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner6MVN",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner6MVN",Inner6MVN);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner6MVTS",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner6MVTS",Inner6MVTS);

  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner20MVM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner20MVM",&Inner20MVM);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner20MV",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner20MV",Inner20MV);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner6MVN",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner20MVN",Inner20MVN);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Inner20MVTS",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Inner20MVTS",Inner20MVTS);

  RootInput::getInstance()->GetChain()->SetBranchStatus("DeltaTVM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("DeltaTVM",&DeltaTVM);
  RootInput::getInstance()->GetChain()->SetBranchStatus("DeltaTV",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("DeltaTV",DeltaTV);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Delta6MVN",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("DeltaTVN",DeltaTVN);
  RootInput::getInstance()->GetChain()->SetBranchStatus("DeltaTVTS",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("DeltaTVTS",DeltaTVTS);

  RootInput::getInstance()->GetChain()->SetBranchStatus("OutersVM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("OutersVM",&OutersVM);
  RootInput::getInstance()->GetChain()->SetBranchStatus("OutersV",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("OutersV",OutersV);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Delta6MVN",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("OutersVN",OutersVN);

  // Output //
  RootOutput::getInstance()->GetTree()->Branch("fTP_X",&fTP_X,"fTP_X/F");
  RootOutput::getInstance()->GetTree()->Branch("fTP_Y",&fTP_Y,"fTP_Y/F");
  RootOutput::getInstance()->GetTree()->Branch("fTS_TMW",&fTS_TMW,"fTS_TMW/l");
  RootOutput::getInstance()->GetTree()->Branch("fTAC_TMW1_FPMW0",&fTAC_TMW1_FPMW0,"fTAC_TMW1_FPMW0/F");
  RootOutput::getInstance()->GetTree()->Branch("fTAC_TMW1_FPMW1",&fTAC_TMW1_FPMW0,"fTAC_TMW1_FPMW1/F");
  RootOutput::getInstance()->GetTree()->Branch("fTAC_TMW2_FPMW0",&fTAC_TMW2_FPMW0,"fTAC_TMW2_FPMW0/F");
  RootOutput::getInstance()->GetTree()->Branch("fTAC_TMW2_FPMW1",&fTAC_TMW2_FPMW1,"fTAC_TMW2_FPMW1/F");
  RootOutput::getInstance()->GetTree()->Branch("fTAC_MW1_PISTA",&fTAC_MW1_PISTA,"fTAC_MW1_PISTA/F");

  RootOutput::getInstance()->GetTree()->Branch("Inner6MVM",&Inner6MVM,"Inner6MVM/I");
  RootOutput::getInstance()->GetTree()->Branch("Inner6MV",Inner6MV);
  RootOutput::getInstance()->GetTree()->Branch("Inner6MVN",Inner6MVN);
  RootOutput::getInstance()->GetTree()->Branch("Inner6MVTS",Inner6MVTS);

  RootOutput::getInstance()->GetTree()->Branch("Inner20MVM",&Inner20MVM,"Inner20MVM/I");
  RootOutput::getInstance()->GetTree()->Branch("Inner20MV",Inner20MV);
  RootOutput::getInstance()->GetTree()->Branch("Inner20MVN",Inner20MVN);
  RootOutput::getInstance()->GetTree()->Branch("Inner20MVTS",Inner20MVTS);

  RootOutput::getInstance()->GetTree()->Branch("DeltaTVM",&DeltaTVM,"DeltaTVM/I");
  RootOutput::getInstance()->GetTree()->Branch("DeltaTV",DeltaTV);
  RootOutput::getInstance()->GetTree()->Branch("DeltaTVN",DeltaTVN);
  RootOutput::getInstance()->GetTree()->Branch("DeltaTVTS",DeltaTVTS);

  RootOutput::getInstance()->GetTree()->Branch("OutersVM",&OutersVM,"OutersVM/I");
  RootOutput::getInstance()->GetTree()->Branch("OutersV",OutersV);
  RootOutput::getInstance()->GetTree()->Branch("OutersVN",OutersVN);


}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  BeamEnergy = -1000;
  Ex240Pu = -1000;
  Ex236U = -1000;
  Ex238U = -1000;
  DeltaE = -1000;
  DeltaEcorr = -1000;
  Eres = -1000;
  Elab = -1000;
  ThetaLab = -1000;
  PhiLab = -1000;
  ThetaCM = -1000;
  XTarget = -1000;
  YTarget = -1000;
  ZTarget = -1000;
  Xcalc = -1000;
  Ycalc = -1000;
  Zcalc = -1000;
  PID = -1000;

  Chio_DE = -1000;
  Chio_E = -1000;
  

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

