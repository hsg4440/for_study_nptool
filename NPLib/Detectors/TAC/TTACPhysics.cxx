/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TTACPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TTACPhysics)


///////////////////////////////////////////////////////////////////////////
TTACPhysics::TTACPhysics()
   : m_EventData(new TTACData),
     m_PreTreatedData(new TTACData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TTACPhysics::AddDetector(TVector3 , string ){
} 

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::AddDetector(double R, double Theta, double Phi, string shape){
} 
  
///////////////////////////////////////////////////////////////////////////
void TTACPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::PreTreat() {
  ClearPreTreatedData();

}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ReadAnalysisConfig() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::Clear() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ReadConfiguration(NPL::InputParser parser) {
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::FillSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::CheckSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TTACPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("TAC", "D"+ NPL::itoa(i+1)+"_ENERGY","TAC_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("TAC", "D"+ NPL::itoa(i+1)+"_TIME","TAC_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("TAC",  true );
  inputChain->SetBranchAddress("TAC", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("TAC", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("TAC", "TTACPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TTACPhysics::Construct() {
  return (NPL::VDetector*) new TTACPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_TAC{
  public:
    proxy_TAC(){
      NPL::DetectorFactory::getInstance()->AddToken("TAC","TAC");
      NPL::DetectorFactory::getInstance()->AddDetector("TAC",TTACPhysics::Construct);
    }
};

proxy_TAC p_TAC;
}

