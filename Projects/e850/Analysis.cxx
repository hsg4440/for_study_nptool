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

  TFile *Fspline=new TFile("./macro/spline_Chio.root");
  for(int i=0;i<34;i++){
    gSpline3[i]=(TSpline3*)Fspline->Get(Form("fspline_%d",i+1));
    Reference[i]=10000-i/34.*(10000-6539);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  double Xoffset = -1.0;
  double Yoffset = 1.0;
  if(fTP_X>-1500 && fTP_Y>-1500){
    XTarget = 0.5*fTP_X + Xoffset;
    YTarget = fTP_Y + Yoffset;
    ZTarget = 1.35;
  }
  else{  
    XTarget = Xoffset;
    YTarget = -2 + Yoffset;
    ZTarget = 1.35;
  }


  TVector3 BeamPosition(XTarget,YTarget,ZTarget);
  TVector3 PositionOnTarget(XTarget,YTarget,ZTarget);
  TVector3 PositionOnTargetNoCorr(Xoffset,-2+Yoffset,0);
  BeamEnergy = 1417.;
  //BeamEnergy = U238C.Slow(BeamEnergy,TargetThickness*0.5,0);
  Transfer10Be->SetBeamEnergy(BeamEnergy);
  Transfer14C->SetBeamEnergy(BeamEnergy);
  Elastic->SetBeamEnergy(BeamEnergy);

  // *** VAMOS *** //
  if(fIC[1]>0 && fIC[5]>0){
    Chio_DE = 0.5*(fIC[0]+fIC[1]+fIC[2]+fIC[3])+fIC[4];
    Chio_E = fIC[5]+fIC[6]+fIC[7]+fIC[8]+fIC[9];
  }
  if(Chio_E>4000 && Chio_E<15000) 
  {
    int Zid=FindClosestZ(Chio_E,Chio_DE); //IdentifyZ(Chio_E,Chio_DE);

    if(Zid==-1) Z = -1000;
    double reference=Reference[Zid];
    double delta=gSpline3[Zid]->Eval(Chio_E) - gSpline3[Zid]->Eval(reference);
    Z= Chio_DE-delta;
  }

  double D1 = fPath-(16.27)/cos(fTP_Theta/1000.)/cos(fTP_Phi/1000.)+(796.56-760.)/cos(fTf/1000.);
  double Toff[20] = {0,0,0,0
		     ,1.4,1.8,1.2,0.6,0.5,2.5,2.1,0.9,1.0,1.6,1.,.9
		     ,0,0,0,0};

  FPMWPat=fFPMWPat_0RawNr[0];

// (//0.0*(fFPMWPat_0RawNr[0]==0)+0.0*(fFPMWPat_0RawNr[0]==1)-0.0*(fFPMWPat_0RawNr[0]==2)+0.0*(fFPMWPat_0RawNr[0]==3)
// 		 +1.4*(fFPMWPat_0RawNr[0]==4)+1.8*(fFPMWPat_0RawNr[0]==5)+1.2*(fFPMWPat_0RawNr[0]==6)+0.6*(fFPMWPat_0RawNr[0]==7) +0.5*(fFPMWPat_0RawNr[0]==8) +2.5*(fFPMWPat_0RawNr[0]==9)+2.1*(fFPMWPat_0RawNr[0]==10)+ 0.9*(fFPMWPat_0RawNr[0]==11)+ 1.0*(fFPMWPat_0RawNr[0]==12)+ 1.6*(fFPMWPat_0RawNr[0]==13)+1.*(fFPMWPat_0RawNr[0]==14)+0.9*(fFPMWPat_0RawNr[0]==15)
// 		 //+0.0*(fFPMWPat_0RawNr[0]==16)+0.0*(fFPMWPat_0RawNr[0]==17)+0.0*(fFPMWPat_0RawNr[0]==18)+0.0*(fFPMWPat_0RawNr[0]==19)
// 		 )*(fFPMWPat_0RawM==1);

  double T1;
  if(fFPMWPat_0RawM==1)
    {
      // cout << "FPMWPat=" << FPMWPat << endl;
      T1 = fT_TMW1_FPMW0_C-21.+Toff[FPMWPat];
      //  else T1=fT_TMW1_FPMW0_C-21.;
      double V1 = D1/T1;
      double Beta1 = V1/29.9792458;
      double Gamma1 = 1./sqrt(1.0-Beta1*Beta1);
      M_Q1 = 1.0*(fBrho/3.105/Beta1/Gamma1);
      E1 = 0.02411*(0.8686*fIC[0]+0.7199*fIC[1]+0.6233*fIC[2]+0.4697*fIC[3]+0.9787*fIC[4]+0.9892*fIC[5]+2.1038*fIC[6]+1.9429*fIC[7]+1.754*fIC[8]+2.5*fIC[9]);
      M1 = E1/931.5016/(Gamma1-1.);
      Q1 = M1/M_Q1;

      Mass = int(Q1+0.5)*M_Q1;
    }

  //cout << PISTA->EventMultiplicity << endl;
  double Energy = 0;
  int strip_DE = 0;
  int strip_E = 0;
  if(PISTA->EventMultiplicity==1){
    DeltaE = PISTA->DE[0];
    Eres = PISTA->back_E[0];
    Energy = DeltaE + Eres;

    strip_DE = PISTA->DE_StripNbr[0];
    strip_E = PISTA->E_StripNbr[0];
    Telescope = PISTA->DetectorNumber[0];
  }
  else if(PISTA->EventMultiplicity==2 && abs(PISTA->DE_StripNbr[0]-PISTA->DE_StripNbr[1])==1){
    DeltaE = PISTA->DE[0] + PISTA->DE[1];
    Eres = PISTA->back_E[0];
    Energy = DeltaE + Eres;

    Telescope = PISTA->DetectorNumber[0];
    if(PISTA->DE[0]>PISTA->DE[1])
      strip_DE = PISTA->DE_StripNbr[0];
    else
      strip_DE = PISTA->DE_StripNbr[1];

    strip_E = PISTA->E_StripNbr[0];
  }
  if(strip_DE>0 && strip_DE<92 && strip_E>0 && strip_E<58){
    TVector3 PISTA_pos = PISTA->GetPositionOfInteraction(Telescope, strip_E, strip_DE);
    TVector3 HitDirection = PISTA_pos - PositionOnTarget;
    PhiLab = PISTA_pos.Phi();
    Xcalc  = PISTA_pos.X();
    Ycalc  = PISTA_pos.Y();
    Zcalc  = PISTA_pos.Z();

    double ThetaLabNoCorr = (PISTA_pos-PositionOnTargetNoCorr).Theta();
    ThetaLab = HitDirection.Angle(TVector3(0,0,1));
    ThetaDetectorSurface = HitDirection.Angle(PISTA->GetDetectorNormal(0));

    DeltaEcorr = DeltaE*cos(ThetaDetectorSurface);
    PID = pow(DeltaEcorr+Eres,1.78)-pow(Eres,1.78);

    ThetaNormalTarget = HitDirection.Angle(TVector3(0,0,1));
    Elab = Energy;//Be10C.EvaluateInitialEnergy(Energy,TargetThickness*0.5,ThetaNormalTarget);
    Ex240Pu = Transfer10Be->ReconstructRelativistic(Elab, ThetaLab);
    Ex240PuNoCorr = Transfer10Be->ReconstructRelativistic(Elab, ThetaLabNoCorr);
    Ex236U  = Transfer14C->ReconstructRelativistic(Elab, ThetaLab);
    Ex238U  = Elastic->ReconstructRelativistic(Elab, ThetaLab);
    ThetaCM = Transfer10Be->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
    ThetaLab = ThetaLab/deg;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("XTarget",&XTarget,"XTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("YTarget",&YTarget,"YTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZTarget",&ZTarget,"ZTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex240Pu",&Ex240Pu,"Ex240Pu/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex240PuNoCorr",&Ex240PuNoCorr,"Ex240PuNoCorr/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex236U",&Ex236U,"Ex236U/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex238U",&Ex238U,"Ex238U/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaE",&DeltaE,"DeltaE/D");
  RootOutput::getInstance()->GetTree()->Branch("DeltaEcorr",&DeltaEcorr,"DeltaEcorr/D");
  RootOutput::getInstance()->GetTree()->Branch("Eres",&Eres,"Eres/D");
  RootOutput::getInstance()->GetTree()->Branch("Telescope",&Telescope,"Telescope/I");
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

  RootOutput::getInstance()->GetTree()->Branch("FPMWPat",&FPMWPat,"FPMWPat/s");
  RootOutput::getInstance()->GetTree()->Branch("M1",&M1,"M1/D");
  RootOutput::getInstance()->GetTree()->Branch("M_Q1",&M_Q1,"M_Q1/D");
  RootOutput::getInstance()->GetTree()->Branch("E1",&E1,"E1/D");
  RootOutput::getInstance()->GetTree()->Branch("Q1",&Q1,"Q1/D");
  RootOutput::getInstance()->GetTree()->Branch("Mass",&Mass,"Mass/D");
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z,"Z/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("IC",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("IC",fIC);

  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_X",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_X",&fTP_X);
  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_Y",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_Y",&fTP_Y);

  RootInput::getInstance()->GetChain()->SetBranchStatus("PISTA_TS",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("PISTA_TS",&fPISTA_TS);
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

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW0_C",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW0_C",&fT_TMW1_FPMW0_C);
  RootInput::getInstance()->GetChain()->SetBranchStatus("FPMWPat_0RawNr",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FPMWPat_0RawNr",fFPMWPat_0RawNr);
  RootInput::getInstance()->GetChain()->SetBranchStatus("FPMWPat_0RawM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FPMWPat_0RawM",&fFPMWPat_0RawM);

  RootInput::getInstance()->GetChain()->SetBranchStatus("Path",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Path",&fPath);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Brho",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Brho",&fBrho);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Pf",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Pf",&fPf);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Tf",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Tf",&fTf);
  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_Theta",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_Theta",&fTP_Theta);
  RootInput::getInstance()->GetChain()->SetBranchStatus("TP_Phi",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("TP_Phi",&fTP_Phi);



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
  RootOutput::getInstance()->GetTree()->Branch("fPISTA_TS",&fPISTA_TS,"fPISTA_TS/l");

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
  Ex240PuNoCorr = -1000;
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
  Telescope = -1;
  Chio_DE = -1000;
  Chio_E = -1000;

  M1 = -1000;
  E1 = -1000;
  M_Q1 = -1000;
  Q1 = -1000;
  Mass = -1000;
  Z = -1000;

}

////////////////////////////////////////////////////////////////////////////////
int Analysis::FindClosestZ(double Chio_E, double Chio_DE)
{
  double DEmin=100000;
  int Zmin=-1;

  for(int i=2;i<34;i++)
  {
    //cout << "i= " << i << endl;
    double test_DE=gSpline3[i]->Eval(Chio_E);
    //cout << test_DE << endl;
    if(fabs(test_DE-Chio_DE)<DEmin && test_DE>1000)
    {
      //cout << "DEmin=" << DEmin << " / test_DE-Chio_DE=" << fabs(test_DE-Chio_DE) << " / test_DE=" << test_DE << " /Chio_DE=" << Chio_DE << endl;
      DEmin=fabs(test_DE-Chio_DE);
      Zmin=i;
    }
    
  }
  return Zmin;
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

