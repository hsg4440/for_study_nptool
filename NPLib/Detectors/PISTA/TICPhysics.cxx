/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace  contact address: pierre.morfouace@cea.fr   *
 *                                                                           *
 * Creation Date  : October 2023                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold IC Treated  data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TICPhysics.h"

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

ClassImp(TICPhysics);

///////////////////////////////////////////////////////////////////////////
TICPhysics::TICPhysics()
  : m_EventData(new TICData),
  m_PreTreatedData(new TICData),
  m_EventPhysics(this),
  m_FPMW_Section(-1),
  m_NumberOfDetectors(0){
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TICPhysics::AddDetector(TVector3 Pos){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
} 

///////////////////////////////////////////////////////////////////////////
void TICPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TICPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::BuildPhysicalEvent() {

  if(m_FPMW_Section<0)
    return;

  Clear();
  PreTreat();

  static CalibrationManager* Cal = CalibrationManager::getInstance();

  int size = m_PreTreatedData->GetICMult();
  for(int i=0; i<size; i++){
    int segment = m_PreTreatedData->GetIC_Section(i);
    double gain = Cal->GetValue("IC/SEC"+NPL::itoa(m_FPMW_Section)+"_SEG"+NPL::itoa(segment)+"_ALIGN",0);

    fIC_raw[i] = m_PreTreatedData->GetIC_Charge(i);
    fIC[i] = gain*m_PreTreatedData->GetIC_Charge(i);

  }

  if(fIC[1]>0 && fIC[5]>0){
    DE = 0.5*(fIC_raw[0] + fIC_raw[1] + fIC_raw[2] + fIC_raw[3]) + fIC_raw[4];
    Eres = fIC_raw[5] + fIC_raw[6] + fIC_raw[7] + fIC_raw[8] + fIC_raw[9];

    double scalor = Cal->GetValue("IC/ETOT_SCALING",0);

    for(int i=0; i<10; i++){
      Etot += fIC[i];
    }

    //DE = 0.5*(fIC[0] + fIC[1] + fIC[2] + fIC[3]) + fIC[4];
    //Eres = fIC[5] + fIC[6] + fIC[7] + fIC[8] + fIC[9];
    //Etot = 0.02411*(0.8686*fIC_raw[0]+0.7199*fIC_raw[1]+0.6233*fIC_raw[2]+0.4697*fIC_raw[3]+0.9787*fIC_raw[4]+0.9892*fIC_raw[5]+2.1038*fIC_raw[6]+1.9429*fIC_raw[7]+1.754*fIC_raw[8]+2.5*fIC_raw[9]); 

  }
  else{
    DE = -100;
    Eres = -100;
    Etot = -100;
  }

  m_FPMW_Section = -1;
}

///////////////////////////////////////////////////////////////////////////
void TICPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  //static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetICMult();
  for (unsigned int i = 0; i < mysize ; ++i) {
    int segment = m_EventData->GetIC_Section(i);
    //cout << section << " " << gain << endl;
    //double charge = gain*m_EventData->GetIC_Charge(i);
    double charge = m_EventData->GetIC_Charge(i);

    m_PreTreatedData->SetIC_Charge(charge);
    m_PreTreatedData->SetIC_Section(segment);
  }
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigIC.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigIC.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigIC.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigIC.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigIC";
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

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
void TICPhysics::Clear() {
  DE = 0;
  Eres = 0;
  Etot = 0;
  for(int i=0; i<11; i++){
    fIC[i] = 0;
    fIC_raw[i] = 0;
  }
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("IC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  IC " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  IC " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  ReadAnalysisConfig();
}


///////////////////////////////////////////////////////////////////////////
void TICPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int section = 0; section<20; section++){
    for(int segment = 0; segment<11; segment++){
      Cal->AddParameter("IC","SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_ALIGN","IC_SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_ALIGN");
    }
  }
  Cal->AddParameter("IC","ETOT_SCALING","IC_ETOT_SCALING");
}

///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("IC",  true );
  inputChain->SetBranchAddress("IC", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("IC", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("IC", "TICPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TICPhysics::Construct() {
  return (NPL::VDetector*) new TICPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_IC{
    public:
      proxy_IC(){
        NPL::DetectorFactory::getInstance()->AddToken("IC","IC");
        NPL::DetectorFactory::getInstance()->AddDetector("IC",TICPhysics::Construct);
      }
  };

  proxy_IC p_IC;
}

