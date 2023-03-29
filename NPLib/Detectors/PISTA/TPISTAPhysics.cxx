/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TPISTAPhysics.h"

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

ClassImp(TPISTAPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TPISTAPhysics::TPISTAPhysics(){
    EventMultiplicity = 0;
    m_EventData = new TPISTAData;
    m_PreTreatedData = new TPISTAData;
    m_EventPhysics = this;
    m_Spectra = NULL;
    m_E_RAW_Threshold = 0; // adc channels
    m_E_Threshold = 0;     // MeV
    m_NumberOfDetectors = 0;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatching = 0.050;
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TPISTAPhysics::AddDetector(TVector3){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddDetector(double R, double Theta, double Phi){
  m_NumberOfDetectors++;

  double Height = 61.8; // mm
  double Base = 78.1; // mm
  double NumberOfStripsX = 62;
  double NumberOfStripsY = 97;
  
  double StripPitchHeight = Height / NumberOfStripsY; // mm
  double StripPitchBase = Base / NumberOfStripsX; // mm


  // Vector U on detector face (parallel to Y strips) Y strips are along X axis
  TVector3 U;
  // Vector V on detector face (parallel to X strips)
  TVector3 V;
  // Vector W normal to detector face (pointing to the back)
  TVector3 W;
  // Vector C position of detector face center
  TVector3 C;
  C = TVector3(R*sin(Theta)*cos(Phi),
        R*sin(Theta)*sin(Phi),
        Height*0.5+R*cos(Theta));

  TVector3 P = TVector3(cos(Theta)*cos(Phi),
      cos(Theta)*sin(Phi),
      -sin(Theta));

  W = C.Unit();
  U = W.Cross(P);
  V = W.Cross(U);

  U = U.Unit();
  V = V.Unit();

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 Corner;
  TVector3 Strip_1_1;
  Strip_1_1 = C - (0.5*Base*U + 0.5*Height*V) + U*(StripPitchBase / 2.) + V*(StripPitchHeight / 2.);

  TVector3 StripPos;
  for(int i=0; i<NumberOfStripsX; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<NumberOfStripsY; j++){
      StripPos = Strip_1_1 + i*U*StripPitchBase + j*V*StripPitchHeight;
      lineX.push_back(StripPos.X());
      lineY.push_back(StripPos.Y());
      lineZ.push_back(StripPos.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
} 

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetPositionOfInteraction(const int i){
  TVector3 Position = TVector3(GetStripPositionX(DetectorNumber[i], E_Strip[i], DE_Strip[i]),
      GetStripPositionY(DetectorNumber[i], E_Strip[i], DE_Strip[i]),
      GetStripPositionZ(DetectorNumber[i], E_Strip[i], DE_Strip[i]));
  
  /*TVector3 Position = TVector3(GetStripPositionX(DetectorNumber[i], DE_StripX[i], E_StripY[i]),
      GetStripPositionY(DetectorNumber[i], DE_StripX[i], E_StripY[i]),
      GetStripPositionZ(DetectorNumber[i], DE_StripX[i], E_StripY[i]));
*/

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetDetectorNormal(const int i){
  TVector3 U = TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1))

    -TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1));

  TVector3 V = TVector3(GetStripPositionX(DetectorNumber[i],62,97),
      GetStripPositionY(DetectorNumber[i],62,97),
      GetStripPositionZ(DetectorNumber[i],62,97))

    -TVector3(GetStripPositionX(DetectorNumber[i],62,1),
      GetStripPositionY(DetectorNumber[i],62,1),
      GetStripPositionZ(DetectorNumber[i],62,1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}
///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  if(1 /*CheckEvent() == 1*/){
    //vector<TVector2> couple = Match_X_Y();
    //EventMultiplicity = couple.size();

    int DEMult = m_PreTreatedData->GetPISTADEMult();
    int EMult = m_PreTreatedData->GetPISTAEMult();

    for(unsigned int i=0; i<DEMult; i++){
      for(unsigned int j=0; j<EMult; j++){
        int DE_DetNbr = m_PreTreatedData->GetPISTA_DE_DetectorNbr(i);
        int E_DetNbr = m_PreTreatedData->GetPISTA_E_DetectorNbr(j);
        if(DE_DetNbr==E_DetNbr){
          int DE_StripNbr = m_PreTreatedData->GetPISTA_DE_StripNbr(i);
          int E_StripNbr = m_PreTreatedData->GetPISTA_E_StripNbr(j);
          
          // Taking Strip energy for DE
          double DE_Energy = m_PreTreatedData->GetPISTA_DE_StripEnergy(i);
          // Taking BAck Energy for E
          double E_Energy = m_PreTreatedData->GetPISTA_E_BackEnergy(j);

          double E_Time = m_PreTreatedData->GetPISTA_E_BackTime(j);

          DetectorNumber.push_back(DE_DetNbr);
          DE_Strip.push_back(DE_StripNbr);
          E_Strip.push_back(E_StripNbr);
          DE.push_back(DE_Energy);
          E.push_back(E_Energy);
          Time.push_back(E_Time);

          PosX.push_back(GetPositionOfInteraction(i).x());
          PosY.push_back(GetPositionOfInteraction(i).y());
          PosZ.push_back(GetPositionOfInteraction(i).z());

        }
      }
    }
  }
  EventMultiplicity = DetectorNumber.size();
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TPISTAPhysics::Match_X_Y(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_DEMult, m_EMult;
  m_DEMult = m_PreTreatedData->GetPISTADEMult();
  m_EMult = m_PreTreatedData->GetPISTAEMult();

  if(m_DEMult>m_MaximumStripMultiplicityAllowed || m_EMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
int TPISTAPhysics::CheckEvent(){
  // Check the size of the different elements
  if(m_PreTreatedData->GetPISTADEMult() == m_PreTreatedData->GetPISTAEMult() )
    return 1;

  else
    return -1;
}


///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // DE
  unsigned int sizeDE = m_EventData->GetPISTADEMult();
  for (UShort_t i = 0; i < sizeDE ; ++i) {
    if (m_EventData->GetPISTA_DE_StripEnergy(i) > m_E_RAW_Threshold) {
      int DetNbr = m_EventData->GetPISTA_DE_DetectorNbr(i);
      int StripNbr = m_EventData->GetPISTA_DE_StripNbr(i);
      double StripE = m_EventData->GetPISTA_DE_StripEnergy(i);
      double BackE  = m_EventData->GetPISTA_DE_BackEnergy(i);
      double StripT  = m_EventData->GetPISTA_DE_StripTime(i);
      double BackT  = m_EventData->GetPISTA_DE_BackTime(i);
      //StripE = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetFirstStage_XE_DetectorNbr(i)),m_EventData->GetFirstStage_XE_Energy(i));
      if (StripE > m_E_Threshold) {
        m_PreTreatedData->SetPISTA_DE(DetNbr, StripNbr, StripE, BackE, StripT, BackT);
      }
    }
  }
 
  // E
  unsigned int sizeE = m_EventData->GetPISTAEMult();
  for (UShort_t i = 0; i < sizeE ; ++i) {
    if (m_EventData->GetPISTA_E_StripEnergy(i) > m_E_RAW_Threshold) {
      int DetNbr = m_EventData->GetPISTA_E_DetectorNbr(i);
      int StripNbr = m_EventData->GetPISTA_E_StripNbr(i);
      double StripE = m_EventData->GetPISTA_E_StripEnergy(i);
      double BackE  = m_EventData->GetPISTA_E_BackEnergy(i);
      double StripT  = m_EventData->GetPISTA_E_StripTime(i);
      double BackT  = m_EventData->GetPISTA_E_BackTime(i);
      //StripE = Cal->ApplyCalibration("PISTA/ENERGY"+NPL::itoa(m_EventData->GetFirstStage_XE_DetectorNbr(i)),m_EventData->GetFirstStage_XE_Energy(i));
      if (StripE > m_E_Threshold) {
        m_PreTreatedData->SetPISTA_E(DetNbr, StripNbr, StripE, BackE, StripT, BackT);
      }
    }
  }
 
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigPISTA.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigPISTA.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigPISTA.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigPISTA.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigPISTA";
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

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
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
void TPISTAPhysics::Clear() {
  EventMultiplicity = 0;

  // Position Information
  PosX.clear();
  PosY.clear();
  PosZ.clear();

  DetectorNumber.clear();
  E.clear();
  DE.clear();
  DE_Strip.clear();
  E_Strip.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("PISTA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");

      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");

      AddDetector(R, Theta, Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitSpectra() {
  m_Spectra = new TPISTASpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TPISTAPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("PISTA", "D"+ NPL::itoa(i+1)+"_ENERGY","PISTA_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("PISTA", "D"+ NPL::itoa(i+1)+"_TIME","PISTA_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("PISTA",  true );
  inputChain->SetBranchAddress("PISTA", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("PISTA", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("PISTA", "TPISTAPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TPISTAPhysics::Construct() {
  return (NPL::VDetector*) new TPISTAPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_PISTA{
    public:
      proxy_PISTA(){
        NPL::DetectorFactory::getInstance()->AddToken("PISTA","PISTA");
        NPL::DetectorFactory::getInstance()->AddDetector("PISTA",TPISTAPhysics::Construct);
      }
  };

  proxy_PISTA p_PISTA;
}

