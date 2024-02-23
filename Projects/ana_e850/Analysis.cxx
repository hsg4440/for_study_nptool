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

  PISTA = (TPISTAPhysics*) m_DetectorManager->GetDetector("PISTA");
  FPMW  = (TFPMWPhysics*) m_DetectorManager->GetDetector("FPMW");
  IC  = (TICPhysics*) m_DetectorManager->GetDetector("IC");
  Tracking = new TVamosReconstruction();
  Tracking->ReadMatrix("MRec.mat");

  InitInputBranch();
  InitOutputBranch();
  Rand = TRandom3();
  LoadCalibParameter();

  TargetThickness = 0.44*micrometer;

  Transfer10Be = new NPL::Reaction("238U(12C,10Be)240Pu@1417");
  Transfer14C  = new NPL::Reaction("238U(12C,14C)236U@1417");
  Elastic      = new NPL::Reaction("238U(12C,12C)238U@1417");

  chain = RootInput::getInstance()->GetChain();

  // Energy loss table
  C12C = EnergyLoss("EnergyLossTable/C14_C.G4table","G4Table",100);
  Be10C = EnergyLoss("EnergyLossTable/Be10_C.G4table","G4Table",100);
  U238C = EnergyLoss("EnergyLossTable/U238_C.G4table","G4Table",100);

  Exo_Energy = new vector<float>();
  Exo_Crystal = new vector<int>();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadCalibParameter(){
  string input_path = "./Calibration/VAMOS/CHIO/charge.cal";

  ifstream ifile;
  ifile.open(input_path.c_str());
  string token;
  double p0, p1;
  int i=0;
  while(ifile>>token>>p0>>p1){
    cout << token << " " << p0 << " " << p1 << endl;
    m_Q_p0[i] = p0 + 0.2;
    m_Q_p1[i] = p1;
    i++;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();

  if(Exo_Mult==1){
    Exogam_Crystal = Exo_Crystal->at(0);
    Exogam_Energy = Exo_Energy->at(0);
  }
  TVector3 PositionOnTarget(0,0,0);
  BeamEnergy = 1417.;
  //BeamEnergy = U238C.Slow(BeamEnergy,TargetThickness*0.5,0);
  Transfer10Be->SetBeamEnergy(BeamEnergy);
  Transfer14C->SetBeamEnergy(BeamEnergy);
  Elastic->SetBeamEnergy(BeamEnergy);

  double Brho_ref = 1.1;
  UShort_t FPMWPat = FPMWPat_0RawNr[0];
  FPMW_Section = FPMWPat;
  IC->SetFPMWSection(FPMW_Section);
  IC->BuildSimplePhysicalEvent();
  //double Toff[20] = {0,0,0,0,1.0,1.3,0.8,0.2,0.0,2.0,1.7,0.4,0.6,1.1,0.5,0.4,3.9,0,0,0};
  double Toff[20] = {0,0,0,0,1.3,1.7,1.2,0.6,0.5,2.5,2.1,0.9,1.2,1.7,1.1,1.1,1.2,0,0,0};
  double Theta = -1000;
  if(FPMW->Xf!=-1000){
    Tracking->CalculateReconstruction(FPMW->Xf, 1000*FPMW->Thetaf, Brho_ref, FF_Brho, Theta, FF_Path);
    // FF_Path is in cm !
    double path1 = FPMW->GetDetectorPositionZ(0)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
    double path2 = (FPMW->GetDetectorPositionZ(2)-7600)/10./cos(FPMW->Thetaf);
    //double path2 = (7965.6-7600)/cos(FPMW->Thetaf)/10;
    FF_D = FF_Path - path1 + path2;
    //FF_T = T12*0.99 + 73.5 + Toff[FPMWPat];// - 0.5;
    FF_T = T02 - 22.6 + Toff[FPMWPat] - 0.3;
    FF_V = FF_D/FF_T;
    FF_Beta = FF_V/29.9792458;
    FF_Gamma = 1./sqrt(1.0 - FF_Beta*FF_Beta);
    FF_AoQ = 1.0*(FF_Brho/3.105/FF_Beta/FF_Gamma);
    FF_Etot = IC->Etot;
    FF_M1 = FF_Etot/931.5016/(FF_Gamma-1);
    FF_Q = FF_M1/FF_AoQ;
    FF_Q = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q;
    int iQ = (int) round(FF_Q);
    FF_Mass = iQ*FF_AoQ;
    //FF_Mass = int(FF_Q+0.5)*FF_AoQ;

 
    double Toff14[20] = {0,0,0,3.6,3.7,4.0,4.1,3.7,1.5,3.6,4.4,4.0,4.8,5.2,4.9,5.1,3.5,2,4.2,2};
    double path3 = FPMW->GetDetectorPositionZ(0)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
    double path4 = (FPMW->GetDetectorPositionZ(3)-7600)/10./cos(FPMW->Thetaf);
    FF_D14 = FF_Path - path3 + path4;
    FF_T14 = (T03 - 12.5)*0.99 + Toff14[FPMWPat];
    FF_V14 = FF_D14/FF_T14;
    double FF_Beta14 = FF_V14/29.9792458;
    double FF_Gamma14 = 1./sqrt(1.0 - FF_Beta14*FF_Beta14);
    FF_AoQ14 = 1.0*(FF_Brho/3.105/FF_Beta14/FF_Gamma14);
    FF_M114 = FF_Etot/931.5016/(FF_Gamma14-1);
    FF_Q14 = FF_M114/FF_AoQ14;
    FF_Q14 = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q14;
    int iQ14 = (int) round(FF_Q14);
    FF_Mass14 = iQ14*FF_AoQ14;
    //FF_Mass14 = int(FF_Q14+0.5)*FF_AoQ14;

    double FF_Dav = 0.5*(FF_D+FF_D14);
    double FF_Tav = 0.5*(FF_T+FF_T14);
    double FF_Vav = FF_Dav/FF_Tav;
    double FF_Betaav = FF_Vav/29.9792458;
    double FF_Gammaav = 1./sqrt(1.0 - FF_Betaav*FF_Betaav);
    double FF_AoQav = 1.0*(FF_Brho/3.105/FF_Betaav/FF_Gammaav);
    double FF_M1av = FF_Etot/931.5016/(FF_Gammaav-1);
    FF_Qav = FF_M1av/FF_AoQav;
    FF_Qav = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Qav;
    int iQav = (int) round(FF_Qav);
    FF_Massav = iQav*FF_AoQav;
    //FF_Massav = int(FF_Qav+0.5)*FF_AoQav;
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

    ThetaLab = HitDirection.Angle(TVector3(0,0,1));
    ThetaDetectorSurface = HitDirection.Angle(PISTA->GetDetectorNormal(0));

    DeltaEcorr = DeltaE*cos(ThetaDetectorSurface);
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
  RootOutput::getInstance()->GetTree()->Branch("Telescope",&Telescope,"Telescope/I");
  RootOutput::getInstance()->GetTree()->Branch("PID",&PID,"PID/D");
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab,"Elab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("PhiLab",&PhiLab,"PhiLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Xcalc",&Xcalc,"Xcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Ycalc",&Ycalc,"Ycalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Zcalc",&Zcalc,"Zcalc/D");
  
  RootOutput::getInstance()->GetTree()->Branch("FF_Brho",&FF_Brho,"FF_Brho/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Path",&FF_Path,"FF_Path/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_D",&FF_D,"FF_D/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T",&FF_T,"FF_T/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V",&FF_V,"FF_V/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ",&FF_AoQ,"FF_AoQ/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Beta",&FF_Beta,"FF_Beta/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Gamma",&FF_Gamma,"FF_Gamma/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q",&FF_Q,"FF_Q/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M1",&FF_M1,"FF_M1/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass",&FF_Mass,"FF_Mass/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_D14",&FF_D14,"FF_D14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T14",&FF_T14,"FF_T14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V14",&FF_V14,"FF_V14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ14",&FF_AoQ14,"FF_AoQ14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q14",&FF_Q14,"FF_Q14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M114",&FF_M114,"FF_M114/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass14",&FF_Mass14,"FF_Mass14/D");
 
  RootOutput::getInstance()->GetTree()->Branch("FF_Qav",&FF_Qav,"FF_Qav/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Massav",&FF_Massav,"FF_Massav/D");
 
  RootOutput::getInstance()->GetTree()->Branch("FF_Etot",&FF_Etot,"FF_Etot/D");
  RootOutput::getInstance()->GetTree()->Branch("FPMW_Section",&FPMW_Section,"FPMW_Section/I");
  
  RootOutput::getInstance()->GetTree()->Branch("Exogam_Energy",&Exogam_Energy,"Exogam_Energy/D");
  RootOutput::getInstance()->GetTree()->Branch("Exogam_Crystal",&Exogam_Crystal,"Exogam_Crystal/I");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW0_FPMW0",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW0_FPMW0",&T02);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW0_FPMW1",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW0_FPMW1",&T03);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW0",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW0",&T12);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW1",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW1",&T13);

  RootInput::getInstance()->GetChain()->SetBranchStatus("FPMWPat_0RawNr",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FPMWPat_0RawNr",FPMWPat_0RawNr);

  RootInput::getInstance()->GetChain()->SetBranchStatus("Exo_Mult",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Exo_Mult",&Exo_Mult);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Exo_Energy",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Exo_Energy",&Exo_Energy);
  RootInput::getInstance()->GetChain()->SetBranchStatus("Exo_Crystal",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Exo_Crystal",&Exo_Crystal);

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
  Telescope = -1;

  FF_Brho = -1;
  FF_Path = -1;
  FF_D = -1;
  FF_T = -1;
  FF_V = -1;
  FF_AoQ = -1;
  FF_Beta = -1;
  FF_Gamma = -1;
  FF_Q = -1;
  FF_M1 = -1;
  FF_Mass = -1;
  FF_Etot = -1;
  FPMW_Section = -1;

  FF_D14 = -1;
  FF_T14 = -1;
  FF_V14 = -1;
  FF_Q14 = -1;
  FF_M114 = -1;
  FF_AoQ14 = -1;
  FF_Mass14 = -1;

  FF_Qav = -1;
  FF_Massav = -1;

  Exogam_Crystal = -1;
  Exogam_Energy = -1;
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

