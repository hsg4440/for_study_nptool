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
  CATS = (TCATSPhysics*)  m_DetectorManager -> GetDetector("CATSDetector");
   
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
  //GATCONFMASTER = **GATCONFMASTER_;
  //return (GATCONFMASTER > 0); 
  //DATATRIG_CATS = **DATATRIG_CATS_;
  //return (DATATRIG_CATS > 0); 
   return true;
}

bool Analysis::FillOutputCondition(){
  return true;
  //return (CATS->MapX.size() > 0 && CATS->MapY.size()>0);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
}



void Analysis::InitOutputBranch() {
  std::cout << "Test output branch /////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){
  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
  //DATATRIG_CATS_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"DATATRIG_CATS");
}
////////////////////////////////////////////////////////////////////////////////

void Analysis::ReInit(){

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

