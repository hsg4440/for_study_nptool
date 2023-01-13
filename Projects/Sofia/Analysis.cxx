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

  InitialConditions = new TInitialConditions();
  InteractionCoordinates = new TInteractionCoordinates();
  FissionConditions = new TFissionConditions();

  InitInputBranch();
  InitOutputBranch();

  m_GladField = new GladFieldMap();
  m_GladField->SetCurrent(2185);
  m_GladField->SetGladEntrance(0,0,-1.1135*m);
  m_GladField->SetGladTiltAngle(-14*deg);
  m_GladField->LoadMap("GladFieldMap_50mm.dat");
  m_GladField->SetBin(50);
  m_GladField->SetTimeStep(0.8);
  m_GladField->SetCentralTheta(-20*deg);

  double Z_MW3 = 4440 * cos(20*deg);
  double X_MW3 = Z_MW3 * tan(20*deg);
  m_GladField->Set_MWPC3_Position(X_MW3,0,Z_MW3);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

  ReInit();

  int mult = InteractionCoordinates->GetDetectedMultiplicity();
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
      double Time = ran.Gaus(InteractionCoordinates->GetTime(i) - BeamTimeOffset,0.02);
      double XE = ran.Gaus(InteractionCoordinates->GetDetectedPositionX(i),0.1);
      double YE = 0;//ran.Gaus(InteractionCoordinates->GetDetectedPositionY(i),0.1);
      double ZE = ran.Gaus(InteractionCoordinates->GetDetectedPositionZ(i),0.1);
      TVector3 vE = TVector3(XE,YE,ZE);
      TVector3 vA = TVector3(0,0,z_target);
      int A = FissionConditions->GetFragmentA(i);
      int Z = FissionConditions->GetFragmentZ(i);
      double Theta = FissionConditions->GetFragmentTheta(i);
      double Phi = FissionConditions->GetFragmentPhi(i);
      double Brho = FissionConditions->GetFragmentBrho(i);
 
      TVector3 dir = TVector3(sin(Theta*deg)*cos(Phi*deg), sin(Theta*deg)*sin(Phi*deg), cos(Theta*deg));
      m_TOF.push_back(Time);
      m_A.push_back(A);
      m_Brho.push_back(Brho);
  
      double Brho_calc = m_GladField->FindBrho(vA,dir,vE);
      double PathLength = m_GladField->GetFlightPath(vA,Brho_calc,vA,dir);
      double velocity = PathLength/Time;
      double beta = velocity * m/ns / NPUNITS::c_light;
      double gamma = 1. / sqrt(1 - beta*beta);
      double AoQ = Brho_calc / (3.107 * beta * gamma);
      double A_calc = AoQ*Z;

      m_Brho_calc.push_back(Brho_calc);
      m_FlightPath.push_back(PathLength);
      m_A.push_back(A_calc);
    }
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

  RootInput::getInstance()->GetChain()->SetBranchStatus("FissionConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fFC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("FissionConditions",&FissionConditions);


}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("TOF",&m_TOF);
  RootOutput::getInstance()->GetTree()->Branch("Brho",&m_Brho);
  RootOutput::getInstance()->GetTree()->Branch("A",&m_A);
  RootOutput::getInstance()->GetTree()->Branch("Brho_calc",&m_Brho_calc);
  RootOutput::getInstance()->GetTree()->Branch("A_calc",&m_A_calc);
  RootOutput::getInstance()->GetTree()->Branch("FlightPath",&m_FlightPath);

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInit(){
  m_TOF.clear();
  m_Brho.clear();
  m_A.clear();

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

