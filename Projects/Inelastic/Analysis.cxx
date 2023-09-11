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
  m_Vendeta= (TVendetaPhysics*) m_DetectorManager->GetDetector("Vendeta");
  InitialConditions = new TInitialConditions();
  InteractionCoordinates = new TInteractionCoordinates();
  ReactionConditions = new TReactionConditions();

  InitInputBranch();
  InitOutputBranch();

  my_Reaction = new NPL::Reaction("1n(238U,1n)238U@2.1");

  neutron = new NPL::Nucleus("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  Einit = InitialConditions->GetKineticEnergy(0);
  double init_ThetaLab = ReactionConditions->GetTheta(0)*deg;
  double init_BeamEnergy = ReactionConditions->GetBeamEnergy();
  //my_Reaction->SetBeamEnergy(init_BeamEnergy);

  neutron->SetKineticEnergy(Einit);
  double beam_TOF = neutron->GetTimeOfFlight();

  double Xtarget = InitialConditions->GetIncidentPositionX();
  double Ytarget = InitialConditions->GetIncidentPositionY();
  double Ztarget = 0;//InitialConditions->GetIncidentPositionZ();
  TVector3 TargetPos = TVector3(Xtarget,Ytarget,Ztarget);

  m_Vendeta->SetAnodeNumber(0);
  m_Vendeta->BuildPhysicalEvent();
  unsigned int size = m_Vendeta->LG_DetectorNumber.size();
  for(int i=0; i<size; i++){
    double Rdet, R;
    Rdet = m_Vendeta->GetDistanceFromTarget(m_Vendeta->LG_DetectorNumber[i]);
    TVector3 DetPos = m_Vendeta->GetVectorDetectorPosition(m_Vendeta->LG_DetectorNumber[i]);
    //TVector3 HitPos = DetPos-TargetPos;
    
    R= Rdet*mm;
    Distance.push_back(R);	
    Det.push_back(m_Vendeta->LG_DetectorNumber[i]); 
    T.push_back(m_Vendeta->LG_Time[i]);
    double T_stop = (m_Vendeta->LG_Time[i])*1e-9;
    neutron->SetTimeOfFlight((T_stop-beam_TOF)/(Rdet*1e-3));
    //neutron->SetTimeOfFlight((T_stop)/(Rdet*1e-3-8e-3));
    Elab.push_back(neutron->GetEnergy());


    double DeltaTheta = atan(63.5/Rdet);
    double exp_ThetaLab = m_Vendeta->GetVectorDetectorPosition(m_Vendeta->LG_DetectorNumber[i]).Theta();
    double random_ThetaLab = ra.Uniform(exp_ThetaLab-DeltaTheta, exp_ThetaLab+DeltaTheta);
    double dEx = my_Reaction->ReconstructRelativistic(Elab[i], random_ThetaLab);

    ThetaLab.push_back(random_ThetaLab/deg);
    //ThetaLab.push_back(exp_ThetaLab/deg);
    Ex.push_back(dEx);
  }
}

///////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Einit",&Einit,"Einit/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("Elab",&Elab);   
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex);   
  RootOutput::getInstance()->GetTree()->Branch("T",&T);   
  RootOutput::getInstance()->GetTree()->Branch("Distance",&Distance);   
  RootOutput::getInstance()->GetTree()->Branch("Det",&Det);   
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);

  RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&ReactionConditions);
}

////////////////////////////////////////////////////////////////////////////////     
void Analysis::ReInitValue(){
  Einit      = -100;
  Ex.clear();
  ThetaLab.clear();
  Elab.clear();
  T.clear();
  Distance.clear();
  Det.clear();
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

