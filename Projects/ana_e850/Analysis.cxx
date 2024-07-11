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
  EXOGAM  = (TExogamPhysics*) m_DetectorManager->GetDetector("Exogam");
  Tracking = new TVamosReconstruction();
  Tracking->ReadMatrix("MRec.mat");

  // Energy loss table
  C12C = EnergyLoss("EnergyLossTable/C12_C.G4table","G4Table",100);
  C12Al = EnergyLoss("EnergyLossTable/C12_Al.G4table","G4Table",100);
  Be10C = EnergyLoss("EnergyLossTable/Be10_C.G4table","G4Table",100);
  U238C = EnergyLoss("EnergyLossTable/U238_C.G4table","G4Table",100);
  geloss_C12C   = new TGraph("EnergyLossTable/C12_C.dat");
  geloss_C12Al  = new TGraph("EnergyLossTable/C12_Al.dat");
  geloss_Be10C  = new TGraph("EnergyLossTable/Be10_C.dat");
  geloss_Be10Al = new TGraph("EnergyLossTable/Be10_Al.dat");



  m_BeamEnergy = 5.955*238;
  m_XTarget_offset = 0;
  m_YTarget_offset = 0;
  m_ZTarget_offset = 0;
  m_Beam_ThetaX = 0;
  m_Beam_ThetaY = 0;
  


  InitInputBranch();
  InitOutputBranch();
  Rand = TRandom3();
  LoadCalibParameter();
  ReadAnalysisConfig();
  LoadTimeOffset();

  TargetThickness = 0.44*micrometer;

  m_BeamEnergy = U238C.Slow(m_BeamEnergy,0.5*TargetThickness,0);
  Transfer10Be = new NPL::Reaction("238U(12C,10Be)240Pu@1417");
  Transfer8Be  = new NPL::Reaction("238U(12C,8Be)242Pu@1417");
  Transfer14C  = new NPL::Reaction("238U(12C,14C)236U@1417");
  Elastic      = new NPL::Reaction("238U(12C,12C)238U@1417");

  Elastic->SetBeamEnergy(m_BeamEnergy);
  Transfer10Be->SetBeamEnergy(m_BeamEnergy);
  Transfer8Be->SetBeamEnergy(m_BeamEnergy);
  Transfer14C->SetBeamEnergy(m_BeamEnergy);

  C12  = new NPL::Particle("12C");
  Be10 = new NPL::Particle("10Be");

  chain = RootInput::getInstance()->GetChain();

  Xmean = 0;
  Ymean = 0;
  Xmean_iter = 0;
  Ymean_iter = 0;
  iteration = 0;
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
    m_Q_p0[i] = p0;
    m_Q_p1[i] = p1;
    i++;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadTimeOffset(){
  ifstream ifile;
  // T13
  string filename = "./macro/AoQ/T13_offset.cal";
  ifile.open(filename.c_str());
  string token;
  double offset;
  int section=0;
  while(ifile>>token>>offset){
    m_T13_Offset[section] = offset;
    section++;
  }
  ifile.close();

  // T14
  filename = "./macro/AoQ/T14_offset.cal";
  ifile.open(filename.c_str());
  section=0;
  while(ifile>>token>>offset){
    m_T14_Offset[section] = offset;
    section++;
  }
  ifile.close();

}


////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();

  VAMOS_TS_hour = fVAMOS_TS_sec/3600.;
  PISTA_TS_hour = fPISTA_TS_sec/3600.;

  VamosAnalysis();
  PistaAnalysis();
  ExogamAnalysis();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ExogamAnalysis(){
  unsigned int size = EXOGAM->E_AB.size();

  //for(unsigned int i=0; i<size; i++){
  if(size==1){
    Exo_Theta = EXOGAM->Theta[0];
    Exo_Phi   = EXOGAM->Phi[0];
    Exo_E     = EXOGAM->E_AB[0];

    if(Exo_Theta!=-1000){
      if(FF_Theta!=-100){
        Exo_cosa  = sin(FF_Theta)*cos(FF_Phi)*sin(Exo_Theta)*cos(Exo_Phi) + sin(FF_Theta)*sin(FF_Phi)*sin(Exo_Theta)*sin(Exo_Phi) + cos(FF_Theta)*cos(Exo_Theta);
        //Exo_cosa  = sin(0)*cos(0)*sin(Exo_Theta)*cos(Exo_Phi) + sin(0)*sin(0)*sin(Exo_Theta)*sin(Exo_Phi) + cos(0)*cos(Exo_Theta);
        Exo_EDC_vamos = EXOGAM->CorrectionDoppler(Exo_Theta,Exo_Phi,FF_Theta,FF_Phi,FF_Beta13,Exo_E);
      }
      if(ThetaLab!=-1000){
        Exo_EDC_pista = EXOGAM->CorrectionDoppler(Exo_Theta,Exo_Phi,ThetaLab,PhiLab,Beta_pista,Exo_E);
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::PistaAnalysis(){  

  //cout << PISTA->EventMultiplicity << endl;
  double Energy = 0;
  if(PISTA->EventMultiplicity==1){
    DeltaE = PISTA->DE[0];
    //Eres = PISTA->E[0];
    Eres = PISTA->back_E[0];
    Energy = DeltaE + Eres;

    strip_DE = PISTA->DE_StripNbr[0];
    strip_E = PISTA->E_StripNbr[0];
    Telescope = PISTA->DetectorNumber[0];
    Time_E = PISTA->back_E_Time[0];
  }

  else if(PISTA->EventMultiplicity==2 && abs(PISTA->DE_StripNbr[0]-PISTA->DE_StripNbr[1])==1){
    DeltaE = PISTA->DE[0] + PISTA->DE[1];
    //Eres = PISTA->E[0];
    Eres = PISTA->back_E[0];
    Energy = DeltaE + Eres;

    Telescope = PISTA->DetectorNumber[0];
    if(PISTA->DE[0]>PISTA->DE[1])
      strip_DE = PISTA->DE_StripNbr[0];
    else
      strip_DE = PISTA->DE_StripNbr[1];

    strip_E = PISTA->E_StripNbr[0];
    Time_E = PISTA->back_E_Time[0];
  }

  if(strip_DE>0 && strip_DE<92 && strip_E>0 && strip_E<58){
    //TVector3 PISTA_pos = PISTA->GetPositionOfInteraction(Telescope, 58-strip_E, strip_DE);
    TVector3 PISTA_pos = PISTA->GetPositionOfInteraction(Telescope, strip_E, strip_DE);
    TVector3 HitDirection = PISTA_pos - PositionOnTarget;
    PhiLab = PISTA_pos.Phi();
    Xcalc  = PISTA_pos.X();
    Ycalc  = PISTA_pos.Y();
    Zcalc  = PISTA_pos.Z();

    TVector3 BeamDirection = TVector3(0,0,1);
    BeamDirection.RotateX(m_Beam_ThetaX*3.1415/180);
    BeamDirection.RotateY(m_Beam_ThetaY*3.1415/180);
    ThetaLab = HitDirection.Angle(BeamDirection);
    ThetaDetectorSurface = HitDirection.Angle(PISTA->GetDetectorNormal(0));

    DeltaEcorr = DeltaE*cos(ThetaDetectorSurface);
    PID = pow(DeltaEcorr+Eres,1.78)-pow(Eres,1.78);

    ThetaNormalTarget = HitDirection.Angle(TVector3(0,0,1));
    Elab = Energy;

    // 12C case //
    double Elab_12C;
    Elab_12C = geloss_C12Al->Eval(Eres);
    Elab_12C = geloss_C12Al->Eval(Elab_12C);
    Elab_12C += DeltaE;
    Elab_12C = geloss_C12Al->Eval(Elab_12C);
    Elab_12C = geloss_C12C->Eval(Elab_12C);
    C12->SetKineticEnergy(Elab_12C);

    // 120Be case //
    double Elab_10Be;
    Elab_10Be = geloss_Be10Al->Eval(Eres);
    Elab_10Be = geloss_Be10Al->Eval(Elab_10Be);
    Elab_10Be += DeltaE;
    Elab_10Be = geloss_Be10Al->Eval(Elab_10Be);
    Elab_10Be = geloss_Be10C->Eval(Elab_10Be);
    Be10->SetKineticEnergy(Elab_10Be);
    Beta_pista = Be10->GetBeta();

    Ex240Pu = Transfer10Be->ReconstructRelativistic(Elab_10Be, ThetaLab);
    Ex236U  = Transfer14C->ReconstructRelativistic(Elab, ThetaLab);
    Ex238U  = Elastic->ReconstructRelativistic(Elab_12C, ThetaLab);
    ThetaCM = Transfer10Be->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
    //ThetaLab = ThetaLab/deg;

    double distance_to_hit = HitDirection.Mag()/1000;
    double C12_tof = C12->GetTimeOfFlight();
    double time_to_hit = distance_to_hit*C12_tof*1e9;
    Pista_Time_Target = Time_E - time_to_hit;
  }

  if(PISTA->EventMultiplicity==4){
    TwoAlphaAnalysis();
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::TwoAlphaAnalysis(){
  m_2alpha = 1;

  double DE1 = PISTA->DE[0];
  double DE2 = PISTA->DE[1];
  double DE3 = PISTA->DE[2];
  double DE4 = PISTA->DE[3];
  double DE = (DE1+DE2+DE3+DE4)/4.;

  double E1 = PISTA->back_E[0];
  double E2 = PISTA->back_E[1];
  double E3 = PISTA->back_E[2];
  double E4 = PISTA->back_E[3];
  double E = (E1+E2+E3+E4)/4.;

  // couple1 // 
  Elab1.push_back(DE1+E1);
  Elab2.push_back(DE3+E2);

  // couple2 //
  Elab1.push_back(DE1+E2);
  Elab2.push_back(DE3+E1);

  int det1 = PISTA->DetectorNumber[0];
  int det2 = PISTA->DetectorNumber[1];
  int strip_DE1 = PISTA->DE_StripNbr[0];
  int strip_DE2 = PISTA->DE_StripNbr[2];

  int strip_E1  = PISTA->E_StripNbr[0];
  int strip_E2  = PISTA->E_StripNbr[1];

  Telescope=det1;

  TVector3 HitDirection;
  if(det1==det2){
    TVector3 PISTA_pos1 = PISTA->GetPositionOfInteraction(det1,strip_E1,strip_DE1);
    TVector3 PISTA_pos2 = PISTA->GetPositionOfInteraction(det1,strip_E1,strip_DE2);
    TVector3 PISTA_pos3 = PISTA->GetPositionOfInteraction(det1,strip_E2,strip_DE1);
    TVector3 PISTA_pos4 = PISTA->GetPositionOfInteraction(det1,strip_E2,strip_DE2);

    TVector3 HitDirection1 = PISTA_pos1 - PositionOnTarget;
    TVector3 HitDirection2 = PISTA_pos2 - PositionOnTarget;
    TVector3 HitDirection3 = PISTA_pos3 - PositionOnTarget;
    TVector3 HitDirection4 = PISTA_pos4 - PositionOnTarget;
    HitDirection = 1./4*(HitDirection1+HitDirection2+HitDirection3+HitDirection4);

    double ThetaLab1 = HitDirection1.Angle(TVector3(0,0,1));
    double ThetaLab2 = HitDirection2.Angle(TVector3(0,0,1));
    double ThetaLab3 = HitDirection3.Angle(TVector3(0,0,1));
    double ThetaLab4 = HitDirection4.Angle(TVector3(0,0,1));

    ThetaLab = (ThetaLab1 + ThetaLab2 + ThetaLab3 + ThetaLab4)/4.;
  }

  /*else{
    TVector3 PISTA_pos1 = PISTA->GetPositionOfInteraction(det1,strip_E1,strip_DE1);
    TVector3 PISTA_pos2 = PISTA->GetPositionOfInteraction(det2,strip_E2,strip_DE2);

    TVector3 HitDirection1 = PISTA_pos1 - PositionOnTarget;
    TVector3 HitDirection2 = PISTA_pos2 - PositionOnTarget;

    double ThetaLab1 = HitDirection1.Angle(TVector3(0,0,1));
    double ThetaLab2 = HitDirection2.Angle(TVector3(0,0,1));

    ThetaLab = (ThetaLab1 + ThetaLab2)/2.;

    }*/

  // *** 8Be construction ** //
  DeltaE = DE;
  ThetaDetectorSurface = HitDirection.Angle(PISTA->GetDetectorNormal(0));
  DeltaEcorr = DeltaE*cos(ThetaDetectorSurface);
  Eres = E;
  Elab = DE + E;
  Ex242Pu = Transfer8Be->ReconstructRelativistic(Elab, ThetaLab);
  ThetaCM = Transfer8Be->EnergyLabToThetaCM(Elab, ThetaLab)/deg;
  ThetaLab = ThetaLab/deg;

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::VamosAnalysis(){

  if(abs(FPMW->Xt)<5 && abs(FPMW->Yt)<5){
    XTarget = FPMW->Xt*cos(FPMW->Theta_in)/cos(FPMW->Theta_in+m_Vamos_Angle*deg);
    YTarget = FPMW->Yt;
    ZTarget = 0;

    FF_Theta = FPMW->Theta_in+m_Vamos_Angle*deg;
    FF_Phi = FPMW->Phi_in;

    Xmean_iter += XTarget;
    Ymean_iter += YTarget;
    iteration++;
    if(iteration%100==0){
      Xmean = Xmean_iter/iteration;
      Ymean = Ymean_iter/iteration;

      iteration = 0;
      Xmean_iter = 0;
      Ymean_iter = 0;
    }
  }
  else{
    XTarget = Xmean;
    YTarget = Ymean;
    ZTarget = 0;
  }

  XTarget += m_XTarget_offset;
  YTarget += m_YTarget_offset;
  ZTarget += m_ZTarget_offset;

  PositionOnTarget = TVector3(XTarget,YTarget,ZTarget);
  if(FPMWPat_0RawM==1){

    UShort_t FPMWPat = FPMWPat_0RawNr[0];
    FPMW_Section = FPMWPat;
    IC->SetFPMWSection(FPMW_Section);
    IC->BuildSimplePhysicalEvent();
    double Theta = -1000;
    if(FPMW->Xf!=-1000){
      FF_DE   = IC->DE;
      FF_Eres = IC->Eres;

      Tracking->CalculateReconstruction(FPMW->Xf, 1000*FPMW->Thetaf, m_Brho_ref, FF_Brho, Theta, FF_Path);
      // FF_Path is in cm ! 

      // T13 //
      double path1 = FPMW->GetDetectorPositionZ(0)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
      double path2 = (FPMW->GetDetectorPositionZ(2)-7600)/10./cos(FPMW->Thetaf);
      FF_D13 = FF_Path - path1 + path2;
      FF_T13 = T13 - 20 + m_T13_Offset[FPMWPat];
      FF_V13 = FF_D13/FF_T13;
      FF_Beta13 = FF_V13/29.9792458;
      FF_Gamma13 = 1./sqrt(1.0 - FF_Beta13*FF_Beta13);
      FF_AoQ13 = 1.0*(FF_Brho/3.10761/FF_Beta13/FF_Gamma13);
      FF_Etot13 = IC->Etot;
      FF_M113 = FF_Etot13/931.5016/(FF_Gamma13-1);
      FF_Q13 = FF_M113/FF_AoQ13;
      //FF_Q13 = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q13;
      int iQ13 = (int) round(FF_Q13);
      FF_Mass13 = iQ13*FF_AoQ13;
      //FF_Mass = int(FF_Q+0.5)*FF_AoQ;

      Vamos_Time_Target = FF_T13 - FPMW->GetDetectorPositionZ(0)/10./FF_V13;

      // T14 //
      path1 = FPMW->GetDetectorPositionZ(0)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
      path2 = (FPMW->GetDetectorPositionZ(3)-7600)/10./cos(FPMW->Thetaf);
      FF_D14 = FF_Path - path1 + path2;
      FF_T14 = T14 - 11 + m_T14_Offset[FPMWPat];
      FF_V14 = FF_D14/FF_T14;
      double FF_Beta14 = FF_V14/29.9792458;
      double FF_Gamma14 = 1./sqrt(1.0 - FF_Beta14*FF_Beta14);
      FF_AoQ14 = 1.0*(FF_Brho/3.10761/FF_Beta14/FF_Gamma14);
      double FF_Etot14 = IC->Etot;
      FF_M114 = FF_Etot14/931.5016/(FF_Gamma14-1);
      FF_Q14 = FF_M114/FF_AoQ14;
      //FF_Q14 = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q14;
      int iQ14 = (int) round(FF_Q14);
      FF_Mass14 = iQ14*FF_AoQ14;

      // T23 //
      path1 = FPMW->GetDetectorPositionZ(1)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
      path2 = (FPMW->GetDetectorPositionZ(2)-7600)/10./cos(FPMW->Thetaf);
      FF_D23 = FF_Path - path1 + path2;
      FF_T23 = T23 + 72;
      FF_V23 = FF_D23/FF_T23;
      double FF_Beta23 = FF_V23/29.9792358;
      double FF_Gamma23 = 1./sqrt(1.0 - FF_Beta23*FF_Beta23);
      FF_AoQ23 = 1.0*(FF_Brho/3.10761/FF_Beta23/FF_Gamma23);
      double FF_Etot23 = IC->Etot;
      FF_M123 = FF_Etot23/931.5016/(FF_Gamma23-1);
      FF_Q23 = FF_M123/FF_AoQ23;
      //FF_Q23 = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q23;
      int iQ23 = (int) round(FF_Q23);
      FF_Mass23 = iQ23*FF_AoQ23;

      // T24 //
      path1 = FPMW->GetDetectorPositionZ(1)/10./cos(FPMW->Theta_in)/cos(FPMW->Phi_in);
      path2 = (FPMW->GetDetectorPositionZ(3)-7600)/10./cos(FPMW->Thetaf);
      FF_D24 = FF_Path - path1 + path2;
      FF_T24 = T24 + 81.7;
      FF_V24 = FF_D24/FF_T24;
      double FF_Beta24 = FF_V24/29.9792458;
      double FF_Gamma24 = 1./sqrt(1.0 - FF_Beta24*FF_Beta24);
      FF_AoQ24 = 1.0*(FF_Brho/3.10761/FF_Beta24/FF_Gamma24);
      double FF_Etot24 = IC->Etot;
      FF_M124 = FF_Etot24/931.5016/(FF_Gamma24-1);
      FF_Q24 = FF_M124/FF_AoQ24;
      //FF_Q24 = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Q24;
      int iQ24 = (int) round(FF_Q24);
      FF_Mass24 = iQ24*FF_AoQ24;


      double FF_Dav = 0.5*(FF_D13+FF_D14);
      double FF_Tav = 0.5*(FF_T13+FF_T14);
      double FF_Vav = FF_Dav/FF_Tav;
      double FF_Betaav = FF_Vav/29.9792458;
      double FF_Gammaav = 1./sqrt(1.0 - FF_Betaav*FF_Betaav);
      double FF_AoQav = 1.0*(FF_Brho/3.10761/FF_Betaav/FF_Gammaav);
      double FF_Etotav = 0.5*(FF_Etot13 + FF_Etot14);
      double FF_M1av = FF_Etotav/931.5016/(FF_Gammaav-1);
      FF_Qav = FF_M1av/FF_AoQav;
      FF_Qav = m_Q_p0[FPMW_Section] + m_Q_p1[FPMW_Section]*FF_Qav;
      int iQav = (int) round(FF_Qav);
      FF_Massav = iQav*FF_AoQav;
      //FF_Massav = int(FF_Qav+0.5)*FF_AoQav;
    }
  }
  else{
    FPMW_Section = -1;
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",&BeamEnergy,"BeamEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("XTarget",&XTarget,"XTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("YTarget",&YTarget,"YTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZTarget",&ZTarget,"ZTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex242Pu",&Ex242Pu,"Ex242Pu/D");
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
  RootOutput::getInstance()->GetTree()->Branch("Beta_pista",&Beta_pista,"Beta_pista/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaDetectorSurface",&ThetaDetectorSurface,"ThetaDetectorSurface/D");
  RootOutput::getInstance()->GetTree()->Branch("PhiLab",&PhiLab,"PhiLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Xcalc",&Xcalc,"Xcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Ycalc",&Ycalc,"Ycalc/D");
  RootOutput::getInstance()->GetTree()->Branch("Zcalc",&Zcalc,"Zcalc/D");
  RootOutput::getInstance()->GetTree()->Branch("strip_DE",&strip_DE,"strip_DE/I");
  RootOutput::getInstance()->GetTree()->Branch("strip_E",&strip_E,"strip_E/I");
  RootOutput::getInstance()->GetTree()->Branch("Time_E",&Time_E,"Time_E/D");

  RootOutput::getInstance()->GetTree()->Branch("Pista_Time_Target",&Pista_Time_Target,"Pista_Time_Target/D");
  RootOutput::getInstance()->GetTree()->Branch("Vamos_Time_Target",&Vamos_Time_Target,"Vamos_Time_Target/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_DE",&FF_DE,"FF_DE/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Eres",&FF_Eres,"FF_Eres/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Z",&FF_Z,"FF_Z/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Theta",&FF_Theta,"FF_Theta/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Phi",&FF_Phi,"FF_Phi/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Brho",&FF_Brho,"FF_Brho/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Path",&FF_Path,"FF_Path/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_D13",&FF_D13,"FF_D13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T13",&FF_T13,"FF_T13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V13",&FF_V13,"FF_V13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ13",&FF_AoQ13,"FF_AoQ13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Beta13",&FF_Beta13,"FF_Beta13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Gamma13",&FF_Gamma13,"FF_Gamma13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q13",&FF_Q13,"FF_Q13/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M113",&FF_M113,"FF_M113/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass13",&FF_Mass13,"FF_Mass13/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_D14",&FF_D14,"FF_D14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T14",&FF_T14,"FF_T14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V14",&FF_V14,"FF_V14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ14",&FF_AoQ14,"FF_AoQ14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q14",&FF_Q14,"FF_Q14/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M114",&FF_M114,"FF_M114/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass14",&FF_Mass14,"FF_Mass14/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_D23",&FF_D23,"FF_D23/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T23",&FF_T23,"FF_T23/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V23",&FF_V23,"FF_V23/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ23",&FF_AoQ23,"FF_AoQ23/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q23",&FF_Q23,"FF_Q23/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M123",&FF_M123,"FF_M123/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass23",&FF_Mass23,"FF_Mass23/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_D24",&FF_D24,"FF_D24/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_T24",&FF_T24,"FF_T24/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_V24",&FF_V24,"FF_V24/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_AoQ24",&FF_AoQ24,"FF_AoQ24/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Q24",&FF_Q24,"FF_Q24/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_M124",&FF_M124,"FF_M124/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Mass24",&FF_Mass24,"FF_Mass24/D");


  RootOutput::getInstance()->GetTree()->Branch("FF_Qav",&FF_Qav,"FF_Qav/D");
  RootOutput::getInstance()->GetTree()->Branch("FF_Massav",&FF_Massav,"FF_Massav/D");

  RootOutput::getInstance()->GetTree()->Branch("FF_Etot13",&FF_Etot13,"FF_Etot13/D");
  RootOutput::getInstance()->GetTree()->Branch("FPMW_Section",&FPMW_Section,"FPMW_Section/I");

  RootOutput::getInstance()->GetTree()->Branch("Elab1",&Elab1);
  RootOutput::getInstance()->GetTree()->Branch("Elab2",&Elab2);
  RootOutput::getInstance()->GetTree()->Branch("m_2alpha",&m_2alpha,"m_2alpha/I");

  RootOutput::getInstance()->GetTree()->Branch("Exo_cosa",&Exo_cosa,"Exo_cosa/D");
  RootOutput::getInstance()->GetTree()->Branch("Exo_E",&Exo_E,"Exo_E/D");
  RootOutput::getInstance()->GetTree()->Branch("Exo_EDC_vamos",&Exo_EDC_vamos,"Exo_EDC_vamos/D");
  RootOutput::getInstance()->GetTree()->Branch("Exo_EDC_pista",&Exo_EDC_pista,"Exo_EDC_pista/D");
  RootOutput::getInstance()->GetTree()->Branch("Exo_Theta",&Exo_Theta,"Exo_Theta/D");
  RootOutput::getInstance()->GetTree()->Branch("Exo_Phi",&Exo_Phi,"Exo_Phi/D");

  RootOutput::getInstance()->GetTree()->Branch("VAMOS_TS_hour",&VAMOS_TS_hour,"VAMOS_TS_hour/D");
  RootOutput::getInstance()->GetTree()->Branch("PISTA_TS_hour",&PISTA_TS_hour,"PISTA_TS_hour/D");


}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW0_FPMW0",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW0_FPMW0",&T13);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW0_FPMW1",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW0_FPMW1",&T14);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW0",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW0",&T23);

  RootInput::getInstance()->GetChain()->SetBranchStatus("T_TMW1_FPMW1",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("T_TMW1_FPMW1",&T24);

  RootInput::getInstance()->GetChain()->SetBranchStatus("FPMWPat_0RawNr",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FPMWPat_0RawNr",FPMWPat_0RawNr);
  RootInput::getInstance()->GetChain()->SetBranchStatus("FPMWPat_0RawM",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FPMWPat_0RawM",&FPMWPat_0RawM);

  RootInput::getInstance()->GetChain()->SetBranchStatus("fVAMOS_TS_sec",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("fVAMOS_TS_sec",&fVAMOS_TS_sec);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fPISTA_TS_sec",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("fPISTA_TS_sec",&fPISTA_TS_sec);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  BeamEnergy = -1000;
  Ex242Pu = -1000;
  Ex240Pu = -1000;
  Ex236U = -1000;
  Ex238U = -1000;
  DeltaE = -1000;
  DeltaEcorr = -1000;
  Eres = -1000;
  Elab = -1000;
  ThetaLab = -1000;
  Beta_pista = -1;
  ThetaDetectorSurface = -1000;
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
  strip_DE = -1;
  strip_E = -1;
  Time_E = -1000;
  Pista_Time_Target = -1000;
  Vamos_Time_Target = -1000;

  Exo_cosa = -100;
  Exo_E = -100;
  Exo_EDC_vamos = -100;
  Exo_EDC_pista = -100;
  Exo_Theta = -100;
  Exo_Phi = -100;

  FF_Theta = -100;
  FF_Phi = -100;
  FF_Brho = -1;
  FF_Path = -1;
  FF_D13 = -1;
  FF_T13 = -1;
  FF_V13 = -1;
  FF_AoQ13 = -1;
  FF_Beta13 = -1;
  FF_Gamma13 = -1;
  FF_Q13 = -1;
  FF_M113 = -1;
  FF_Mass13 = -1;
  FF_Etot13 = -1;
  FPMW_Section = -1;

  FF_D14 = -1;
  FF_T14 = -1;
  FF_V14 = -1;
  FF_Q14 = -1;
  FF_M114 = -1;
  FF_AoQ14 = -1;
  FF_Mass14 = -1;

  FF_D23 = -1;
  FF_T23 = -1;
  FF_V23 = -1;
  FF_Q23 = -1;
  FF_M123 = -1;
  FF_AoQ23 = -1;
  FF_Mass23 = -1;

  FF_D24 = -1;
  FF_T24 = -1;
  FF_V24 = -1;
  FF_Q24 = -1;
  FF_M124 = -1;
  FF_AoQ24 = -1;
  FF_Mass24 = -1;


  FF_Qav = -1;
  FF_Massav = -1;

  m_2alpha = 0;
  Elab1.clear();
  Elab2.clear();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReadAnalysisConfig(){
  bool ReadingStatus = false;

  string filename = "AnalysisConfig.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(filename.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No AnalysisConfig.dat found: Default parameter loaded for Analayis " << filename << endl;
    return;
  }
  cout << "**** Loading user parameter for Analysis from AnalysisConfig.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% AnalysisConfig.dat %%%");
  asciiConfig->Append(filename.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "AnalysisConfig";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="XTARGET_OFFSET") {
        AnalysisConfigFile >> DataBuffer;
        m_XTarget_offset = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_XTarget_offset << endl;
      }
      else if (whatToDo=="YTARGET_OFFSET") {
        AnalysisConfigFile >> DataBuffer;
        m_YTarget_offset = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_YTarget_offset << endl;
      }
      else if (whatToDo=="ZTARGET_OFFSET") {
        AnalysisConfigFile >> DataBuffer;
        m_ZTarget_offset = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_ZTarget_offset << endl;
      }
      else if (whatToDo=="BEAM_THETAX") {
        AnalysisConfigFile >> DataBuffer;
        m_Beam_ThetaX = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_Beam_ThetaX << endl;
      }
      else if (whatToDo=="BEAM_THETAY") {
        AnalysisConfigFile >> DataBuffer;
        m_Beam_ThetaY = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_Beam_ThetaY << endl;
      }
      else if (whatToDo=="BRHO_REF") {
        AnalysisConfigFile >> DataBuffer;
        m_Brho_ref = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_Brho_ref << endl;
      }
      else if (whatToDo=="VAMOS_ANGLE") {
        AnalysisConfigFile >> DataBuffer;
        m_Vamos_Angle = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << m_Vamos_Angle << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
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

