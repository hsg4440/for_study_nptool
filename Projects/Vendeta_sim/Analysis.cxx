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
 *  This class describe  Vendeta analysis project                       *
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
#include"RootInput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){

  Vendeta= (TVendetaPhysics*) m_DetectorManager->GetDetector("Vendeta");
  FC= (TFissionChamberPhysics*) m_DetectorManager->GetDetector("FissionChamber");
  InitialConditions = new TInitialConditions();

  InitInputBranch();
  InitOutputBranch();
  
  neutron = new NPL::Particle("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

  ReInitValue();

  double incomingE = InitialConditions->GetKineticEnergy(0);  
  inEnergy = incomingE;

  unsigned int FC_mult = FC->AnodeNumber.size();
  unsigned int HF_mult = FC->Time_HF.size();

  double GammaOffset[11]={0,0,0,0,0,0,0,0,0,0,0};
  double FakeFission_Offset = 0;

  double incomingDT=0;
  double flight_path = 21500.;
  if(FC_mult>0){
    int EventMax=0;

    int anode = FC->AnodeNumber[EventMax];
    double Time_FC = FC->Time[EventMax];
    bool isFake = FC->isFakeFission[EventMax];
    double Time_HF = FC->Time_HF[EventMax];
    incomingDT = FC->Time[EventMax] - FC->Time_HF[EventMax];
    TVector3 AnodePos = FC->GetVectorAnodePosition(6);

    if(anode>0){ 
      incomingDT -= GammaOffset[anode-1];
      AnodePos = FC->GetVectorAnodePosition(anode);
    }
    else if(anode ==0){
      incomingDT -= FakeFission_Offset;
    }

    if(incomingDT<0){
      incomingDT += 1790;
    }

    double length = flight_path;// + 6*FC->AnodeNumber[i];	
    neutron->SetBeta((length/incomingDT) / c_light);

    inToF.push_back(incomingDT);

    FC_Q1.push_back(FC->Q1[EventMax]);
    FC_Q2.push_back(FC->Q2[EventMax]);
    FC_Qmax.push_back(FC->Qmax[EventMax]);
    FC_DT.push_back(FC->DT_FC[EventMax]);
    FC_FakeFission.push_back(FC->isFakeFission[EventMax]);
    FC_Anode_ID.push_back(anode);
  }
  Vendeta->SetAnodeNumber(6);
  Vendeta->BuildPhysicalEvent();

  // VENDETA LG 
  unsigned int Vendeta_LG_mult = Vendeta->LG_DetectorNumber.size();
  for(unsigned int i=0; i<Vendeta_LG_mult; i++){

    int DetNbr          = Vendeta->LG_DetectorNumber[i];
    double Time_Vendeta = Vendeta->LG_Time[i];
    TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);
    double Rdet         = DetPos.Mag();
    //double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendlg[DetNbr-1];
    double DT = Time_Vendeta;

    double DeltaTheta = atan(63.5/Rdet);
    double Theta_Vendeta = DetPos.Theta();
    double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
    //neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
    //neutron->SetTimeOfFlight(DT*1e-9/(0.55));
    neutron->SetBeta(  (Rdet/DT) / c_light); 

    double En = neutron->GetEnergy();

    // Filling output tree
    LG_Tof.push_back(DT);
    LG_ID.push_back(DetNbr);
    LG_ELab.push_back(En);
    LG_ThetaLab.push_back(Theta_random);
    LG_Q1.push_back(Vendeta->LG_Q1[i]);
    LG_Q2.push_back(Vendeta->LG_Q2[i]);
    LG_Qmax.push_back(Vendeta->LG_Qmax[i]);

  }

  // VENDETA HG 
  unsigned int Vendeta_HG_mult = Vendeta->HG_DetectorNumber.size();
  for(unsigned int i=0; i<Vendeta_HG_mult; i++){
    int DetNbr          = Vendeta->HG_DetectorNumber[i];
    double Time_Vendeta = Vendeta->HG_Time[i];
    TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);
    double Rdet         = DetPos.Mag();
    //double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendhg[DetNbr-1];
    double DT = Time_Vendeta;


    double DeltaTheta = atan(63.5/Rdet);

    double Theta_Vendeta = DetPos.Theta();
    double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
    //neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
    //neutron->SetTimeOfFlight(DT*1e-9/(0.55));
    neutron->SetBeta( (Rdet/DT) / c_light); 
    double En = neutron->GetEnergy();
    // Filling output tree
    HG_ID.push_back(DetNbr);
    HG_Tof.push_back(DT);
    HG_ELab.push_back(En);
    HG_ThetaLab.push_back(Theta_random);
    HG_Q1.push_back(Vendeta->HG_Q1[i]);
    HG_Q2.push_back(Vendeta->HG_Q2[i]);
    HG_Qmax.push_back(Vendeta->HG_Qmax[i]);
  }

  //Highlight saturated detectors
  static vector<int> LG_Saturated, HG_Saturated, LG_T_sat, HG_T_sat; 
  LG_Saturated.clear(), HG_Saturated.clear(), LG_T_sat.clear(), HG_T_sat.clear();
  for(int j = 0; j < LG_Tof.size();j++){
    int lgID = LG_ID[j];
    if( (lgID %4 < 2) && LG_Q1[j]>500e3){ // 50 Ohm
      LG_Saturated.push_back(lgID);
      LG_T_sat.push_back(LG_Tof[j]);
    }  
    else if( lgID %4 > 1 && LG_Q1[j]>590e3){ // 70 Ohm
      LG_Saturated.push_back(lgID);
      LG_T_sat.push_back(LG_Tof[j]);
    }
  }
  for(int j = 0; j < HG_Tof.size();j++){
    int hgID = HG_ID[j] ;
    if( (hgID %4 ==0 || hgID %4 ==1) && HG_Qmax[j] > 23800){
      HG_Saturated.push_back(hgID);
      HG_T_sat.push_back(HG_Tof[j]);
    }  
    else if( hgID %4 !=0 && hgID %4 !=1 && HG_Qmax[j] > 24300){
      HG_Saturated.push_back(hgID);
      HG_T_sat.push_back(HG_Tof[j]);
    }       
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions",true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*",true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  // Incoming neutron
  RootOutput::getInstance()->GetTree()->Branch("inToF",&inToF);
  RootOutput::getInstance()->GetTree()->Branch("inEnergy",&inEnergy,"inEnergy/D");

  // FissionChamber
  RootOutput::getInstance()->GetTree()->Branch("FC_Q1",&FC_Q1);
  RootOutput::getInstance()->GetTree()->Branch("FC_Q2",&FC_Q2);
  RootOutput::getInstance()->GetTree()->Branch("FC_Qmax",&FC_Qmax);
  RootOutput::getInstance()->GetTree()->Branch("FC_DT",&FC_DT);
  RootOutput::getInstance()->GetTree()->Branch("FC_FakeFission",&FC_FakeFission);
  RootOutput::getInstance()->GetTree()->Branch("FC_Anode_ID",&FC_Anode_ID);
  /* RootOutput::getInstance()->GetTree()->Branch("FC_Time",&FC_Time); */
  /* RootOutput::getInstance()->GetTree()->Branch("HF_Time",&HF_Time); */

  // LG 
  RootOutput::getInstance()->GetTree()->Branch("LG_ID",&LG_ID);
  RootOutput::getInstance()->GetTree()->Branch("LG_ThetaLab",&LG_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("LG_ELab",&LG_ELab);
  RootOutput::getInstance()->GetTree()->Branch("LG_Tof",&LG_Tof);
  RootOutput::getInstance()->GetTree()->Branch("LG_Q1",&LG_Q1);
  RootOutput::getInstance()->GetTree()->Branch("LG_Q2",&LG_Q2);
  RootOutput::getInstance()->GetTree()->Branch("LG_Qmax",&LG_Qmax);
  /* RootOutput::getInstance()->GetTree()->Branch("LG_Time",&LG_Time); */

  // HG
  RootOutput::getInstance()->GetTree()->Branch("HG_ID",&HG_ID);
  RootOutput::getInstance()->GetTree()->Branch("HG_ThetaLab",&HG_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("HG_ELab",&HG_ELab);
  RootOutput::getInstance()->GetTree()->Branch("HG_Tof",&HG_Tof);
  RootOutput::getInstance()->GetTree()->Branch("HG_Q1",&HG_Q1);
  RootOutput::getInstance()->GetTree()->Branch("HG_Q2",&HG_Q2);
  RootOutput::getInstance()->GetTree()->Branch("HG_Qmax",&HG_Qmax);
  /* RootOutput::getInstance()->GetTree()->Branch("HG_Time",&HG_Time); */
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  inToF.clear();
  inEnergy = 0;

  LG_ThetaLab.clear();
  LG_ELab.clear();
  LG_Tof.clear();
  LG_ID.clear();
  LG_Q1.clear();
  LG_Q2.clear();
  LG_Qmax.clear();
  /* LG_Time.clear(); */

  HG_ThetaLab.clear();
  HG_ELab.clear();
  HG_Tof.clear();
  HG_ID.clear();
  HG_Q1.clear();
  HG_Q2.clear();
  HG_Qmax.clear();
  /* HG_Time.clear(); */

  FC_Q1.clear();
  FC_Q2.clear();
  FC_Qmax.clear();
  FC_DT.clear();
  FC_FakeFission.clear();
  FC_Anode_ID.clear();
  /* FC_Time.clear(); */
  /* HF_Time.clear(); */
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

