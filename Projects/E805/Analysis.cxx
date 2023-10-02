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
 *  This class describe  MUST_AND_ZDD analysis project                       *
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
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  InitInputBranch();
  InitOutputBranch();
  //CATS = (TCATSPhysics*)  m_DetectorManager -> GetDetector("CATSDetector");
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
   
  reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction->GetBeamEnergy();


  string Path = "../../Inputs/EnergyLoss/";
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string beam=  NPL::ChangeNameToG4Standard(reaction->GetNucleus1()->GetName());
  string heavy_ejectile=  NPL::ChangeNameToG4Standard(reaction->GetNucleus4()->GetName());
  string light=NPL::ChangeNameToG4Standard(reaction->GetNucleus3()->GetName());



//
  ProtonSi = NPL::EnergyLoss(Path+ "proton_Si.G4table", "G4Table", 100);
}
  ///////////////////////////// Initialize some important parameters //////////////////////////////////


bool Analysis::UnallocateBeforeBuild(){
  //std::cout << "test unallocate" << std::endl;
  GATCONFMASTER = **GATCONFMASTER_;
  return (GATCONFMASTER > 0); 
  //DATATRIG_CATS = **DATATRIG_CATS_;
  //return (DATATRIG_CATS > 0); 
  //return true;
}

bool Analysis::FillOutputCondition(){
  return true;
  //return (CATS->MapX.size() > 0 && CATS->MapY.size()>0);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  
    ReInit();
  
  //for(unsigned int countMust2 = 0 ; countMust2 < M2->Si_E.size() ; countMust2++){
  //Si_E_M2 = M2->Si_E[countMust2];
  //CsI_E_M2 = M2->CsI_E[countMust2];
  //ThetaM2Surface = 0;
  //  if(Si_E_M2 > 0 && CsI_E_M2 > 8192){
  //    //double EfromdeltaE = ProtonSi.EvaluateEnergyFromDeltaE(
  //    //  Si_E_M2, 300*um, ThetaM2Surface, 6.0 * MeV, 300.0 * MeV,
  //    //  0.001 * MeV, 10000);
  //    double EfromdeltaE = (CsI_E_M2-8192)*0.1;
  //    M2_ECsI_from_deltaE.push_back(EfromdeltaE);
  //    if(EfromdeltaE > 0){
  //      Beta_light = sqrt(1./(1.+1./(pow(EfromdeltaE/911. + 1,2)-1)));
  //      Beta_from_deltaE.push_back(Beta_light);
  //      if(Beta_light>0){
  //        double Beth = log(2*511.*Beta_light*Beta_light/(0.174*(1-Beta_light*Beta_light))) - Beta_light*Beta_light; 
  //        Beth_from_deltaE.push_back(Beth);
  //      }
  //    }
  //  }
  //}
}



void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("M2_ECsI_from_deltaE",&M2_ECsI_from_deltaE);
  RootOutput::getInstance()->GetTree()->Branch("Beta_from_deltaE",&Beta_from_deltaE);
  RootOutput::getInstance()->GetTree()->Branch("Beth_from_deltaE",&Beth_from_deltaE);
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){
  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
  GATCONFMASTER_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"GATCONFMASTER");
  //DATATRIG_CATS_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"DATATRIG_CATS");
}
////////////////////////////////////////////////////////////////////////////////

void Analysis::ReInit(){
  //M2_ECsI_from_deltaE.clear();
  //Beta_from_deltaE.clear();
  //Beth_from_deltaE.clear();
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

