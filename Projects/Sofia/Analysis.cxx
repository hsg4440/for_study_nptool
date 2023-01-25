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
 *  This class describe  Sofia analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include "Analysis.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "RootInput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  //Sofia= (TSofiaPhysics*) m_DetectorManager->GetDetector("Sofia");

  beam_ana = true;
  InitialConditions = new TInitialConditions();
  InteractionCoordinates = new TInteractionCoordinates();
  
  /*if(beam_ana==false)
    FissionConditions = new TFissionConditions();
*/

  InitInputBranch();
  InitOutputBranch();

  m_GladField = new GladFieldMap();
  m_GladField->SetCurrent(2185);
  m_GladField->SetGladEntrance(0,0,-1113.5*mm);
  m_GladField->SetGladTiltAngle(-0*deg);
  m_GladField->LoadMap("GladFieldMap_50mm.dat");
  m_GladField->SetBin(50);
  m_GladField->SetTimeStep(0.8);
  m_GladField->SetCentralTheta(-20*deg);

  double R_MW3 = 4440*mm;
  double Z_MW3 = R_MW3 * cos(20*deg);
  double X_MW3 = -R_MW3 * sin(20*deg);

  m_GladField->Set_MWPC3_Position(X_MW3,0,Z_MW3);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

  ReInit();

  if(beam_ana==true)
    BeamAnalysis();
  else 
    FFAnalysis();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::FFAnalysis(){

  /*int mult = InteractionCoordinates->GetDetectedMultiplicity();
  if(mult==2){
    int A_CN = FissionConditions->GetA_CN();
    int Z_CN = FissionConditions->GetZ_CN();
    double Elab_CN = FissionConditions->GetELab_CN();
    double z_start_beam = InitialConditions->GetIncidentPositionZ();
    double z_target = -4093;
    TVector3 v1 = TVector3(0,0,z_target-z_start_beam);
    TVector3 v2 = TVector3(0,0,z_start_beam);
    NPL::Particle* beam = new NPL::Particle(Z_CN,A_CN);
    beam->SetKineticEnergy(Elab_CN);
    double beta = beam->GetBeta();
    double BeamTimeOffset = v1.Mag()*mm/(beta*NPUNITS::c_light);
    for(int i=0; i<mult; i++){
      double Time = ran.Gaus(InteractionCoordinates->GetTime(i) - BeamTimeOffset,0.);
      m_XE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionX(i),0));
      m_YE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionY(i),0));
      m_ZE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionZ(i),0));
      TVector3 vE = TVector3(m_XE[i],m_YE[i],m_ZE[i]);
      //TVector3 vE = TVector3(m_XE[i],0,m_ZE[i]);
      TVector3 vA = TVector3(0,0,z_target);
      int Asim = FissionConditions->GetFragmentA(i);
      int Adet = InteractionCoordinates->GetA(i);
      //int Z = FissionConditions->GetFragmentZ(i);
      int Z = InteractionCoordinates->GetZ(i);
      double ThetaIn = FissionConditions->GetFragmentTheta(i);
      if(FissionConditions->GetFragmentMomentumX(i)<0)
        ThetaIn = -ThetaIn;
      double Phi = FissionConditions->GetFragmentPhi(i);
      double Brho = FissionConditions->GetFragmentBrho(i);
 
      TVector3 dir = TVector3(sin(ThetaIn*deg)*cos(Phi*deg), sin(ThetaIn*deg)*sin(Phi*deg), cos(ThetaIn*deg));
      //TVector3 dir = TVector3(sin(ThetaIn*deg), 0, cos(ThetaIn*deg));
      m_TOF.push_back(Time);
      m_Adet.push_back(Adet);
      m_Asim.push_back(Asim);
      m_Brho.push_back(Brho);
      m_ThetaIn.push_back(ThetaIn);
 

      double Brho_calc = m_GladField->FindBrho(vA,dir,vE);
      double PathLength = m_GladField->GetFlightPath(vA,Brho_calc,vA,dir)/1000;
      double velocity = PathLength/Time;
      double beta = velocity * m/ns / NPUNITS::c_light;
      double gamma = 1. / sqrt(1 - beta*beta);
      double AoQ = Brho_calc / (3.107 * beta * gamma);
      double A_calc = AoQ*Z;

      m_Brho_calc.push_back(Brho_calc);
      m_FlightPath.push_back(PathLength);
      m_A_calc.push_back(A_calc);
    }
  }*/
    
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamAnalysis(){

  int mult = InteractionCoordinates->GetDetectedMultiplicity();
  if(mult==1){
    double Elab_CN = InitialConditions->GetIncidentInitialKineticEnergy();
    double x_start_beam = InitialConditions->GetIncidentPositionX();
    double y_start_beam = InitialConditions->GetIncidentPositionY();
    double z_start_beam = InitialConditions->GetIncidentPositionZ();
    
    int Adet = InteractionCoordinates->GetA(0);
    int Zdet = InteractionCoordinates->GetZ(0);
    NPL::Particle* beam = new NPL::Particle(Zdet,Adet);
    beam->SetKineticEnergy(Elab_CN);
    //double beta = beam->GetBeta();
    double Brho = beam->GetBrho();
    double Time = ran.Gaus(InteractionCoordinates->GetTime(0),0.);
    m_XE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionX(0),0));
    m_YE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionY(0),0));
    m_ZE.push_back(ran.Gaus(InteractionCoordinates->GetDetectedPositionZ(0),0));
    TVector3 vE = TVector3(m_XE[0],m_YE[0],m_ZE[0]);
    TVector3 vA = TVector3(x_start_beam,y_start_beam,z_start_beam);

    double ThetaIn = InitialConditions->GetIncidentEmittanceTheta();
    double Phi = InitialConditions->GetIncidentEmittancePhi();

    TVector3 dir = TVector3(sin(ThetaIn*deg)*cos(Phi*deg), sin(ThetaIn*deg)*sin(Phi*deg), cos(ThetaIn*deg));
    m_TOF.push_back(Time);
    m_Adet.push_back(Adet);
    m_Brho.push_back(Brho);
    m_ThetaIn.push_back(ThetaIn);

    double Brho_calc = m_GladField->FindBrho(vA,dir,vE);
    double PathLength = m_GladField->GetFlightPath(vA,Brho_calc,vA,dir)/1000;
    double velocity = PathLength/Time;
    double beta = velocity * m/ns / NPUNITS::c_light;
    double gamma = 1. / sqrt(1 - beta*beta);
    double AoQ = Brho_calc / (3.107 * beta * gamma);
    double A_calc = AoQ*Zdet;

    m_Brho_calc.push_back(Brho_calc);
    m_FlightPath.push_back(PathLength);
    m_A_calc.push_back(A_calc);
  }

}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);

  RootInput::getInstance()->GetChain()->SetBranchStatus("InteractionCoordinates",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fDetecetd_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InteractionCoordinates",&InteractionCoordinates);

  /*if(beam_ana==false){
    RootInput::getInstance()->GetChain()->SetBranchStatus("FissionConditions",true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("fFC_*",true);
    RootInput::getInstance()->GetChain()->SetBranchAddress("FissionConditions",&FissionConditions);
  }*/

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("TOF",&m_TOF);
  RootOutput::getInstance()->GetTree()->Branch("ThetaIn",&m_ThetaIn);
  RootOutput::getInstance()->GetTree()->Branch("Brho",&m_Brho);
  RootOutput::getInstance()->GetTree()->Branch("Asim",&m_Asim);
  RootOutput::getInstance()->GetTree()->Branch("Adet",&m_Adet);
  RootOutput::getInstance()->GetTree()->Branch("Brho_calc",&m_Brho_calc);
  RootOutput::getInstance()->GetTree()->Branch("A_calc",&m_A_calc);
  RootOutput::getInstance()->GetTree()->Branch("FlightPath",&m_FlightPath);
  RootOutput::getInstance()->GetTree()->Branch("XE",&m_XE);
  RootOutput::getInstance()->GetTree()->Branch("YE",&m_YE);
  RootOutput::getInstance()->GetTree()->Branch("ZE",&m_ZE);

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInit(){
  m_TOF.clear();
  m_Brho.clear();
  m_Asim.clear();
  m_Adet.clear();
  m_ThetaIn.clear();
  m_XE.clear();
  m_YE.clear();
  m_ZE.clear();

  m_Brho_calc.clear();
  m_A_calc.clear();
  m_FlightPath.clear();
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

