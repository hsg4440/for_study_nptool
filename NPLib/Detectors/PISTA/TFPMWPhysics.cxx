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
 *  This class hold FPMW Treated  data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TFPMWPhysics.h"

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

ClassImp(TFPMWPhysics);

///////////////////////////////////////////////////////////////////////////
TFPMWPhysics::TFPMWPhysics()
  : m_EventData(new TFPMWData),
  m_PreTreatedData(new TFPMWData),
  m_EventPhysics(this),
  m_NumberOfDetectors(0){
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TFPMWPhysics::AddDetector(TVector3 Pos){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  DetPosX.push_back(Pos.X());
  DetPosY.push_back(Pos.Y());
  DetPosZ.push_back(Pos.Z());

  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::BuildPhysicalEvent() {
  PreTreat();

  unsigned int sizeX = m_PreTreatedData->GetFPMWMultX();
  for(unsigned int i=0; i<sizeX; i++){
    DetectorHitX.insert(m_PreTreatedData->GetFPMW_DetX(i));
  }

  unsigned int sizeY = m_PreTreatedData->GetFPMWMultY();
  for(unsigned int i=0; i<sizeY; i++){
    if(DetectorHitX.find(m_PreTreatedData->GetFPMW_DetY(i)) != DetectorHitX.end()){
      DetectorHit.insert(m_PreTreatedData->GetFPMW_DetY(i));
    }
  }
  unsigned int sizeDet = DetectorHit.size();


  for(unsigned int i=0; i<sizeX; i++){
    int StrX = m_PreTreatedData->GetFPMW_StripX(i);
    int NX = m_PreTreatedData->GetFPMW_DetX(i);
    double QX = m_PreTreatedData->GetFPMW_ChargeX(i);
    if(DetectorHit.find(NX) != DetectorHit.end()){
      MapX[NX].push_back(std::make_pair(StrX,QX));
      QSumX[NX] += QX;
      if(MaxQX.find(NX)==MaxQX.end() || MaxQX[NX].second<QX){
        MaxQX[NX] = make_pair(StrX,QX);
      }
    }
  }
  for(unsigned int i=0; i<sizeY; i++){
    int StrY = m_PreTreatedData->GetFPMW_StripY(i);
    int NY = m_PreTreatedData->GetFPMW_DetY(i);
    double QY = m_PreTreatedData->GetFPMW_ChargeY(i);
    if(DetectorHit.find(NY) != DetectorHit.end()){
      MapY[NY].push_back(std::make_pair(StrY,QY));
      QSumY[NY] += QY;
      if(MaxQY.find(NY)==MaxQY.end() || MaxQY[NY].second<QY){
        MaxQY[NY] = make_pair(StrY,QY);
      }
    }
  }

  for(auto &DetN : DetectorHit){
    double PosX = AnalyticHyperbolicSecant(MaxQX[DetN],MapX[DetN]);
    double PosY = AnalyticHyperbolicSecant(MaxQY[DetN],MapY[DetN]);

    int sx0 = (int) PosX;
    int sx1 = sx0+1;
    int sy0 = (int) PosY;
    int sy1 = sy0+1;

    DetectorNbr.push_back(DetN);
    ChargeX.push_back(QSumX[DetN]);
    ChargeY.push_back(QSumY[DetN]);
    PositionX.push_back(PosX);
    PositionY.push_back(PosY);
  }

}
  
/////////////////////////////////////////////////////////////////////
double TFPMWPhysics::AnalyticHyperbolicSecant(std::pair<int,double>& MaxQ,std::vector<std::pair<int,double>>& Map){
  double sech = -1000 ;

  // std::cout << "test AnH 1" << std::endl;
  if(MaxQ.second > 0 && MaxQ.first > 2 && MaxQ.first<996){	
    //      if(Buffer_Q[MaxQ.first-1+1]==0||Buffer_Q[MaxQ.first-1-1]==0)
    //        return sech;
    // std::cout << "test AnH 2" << std::endl;
    double q2 = MaxQ.second;
    double q1 = 0,q3 = 0;
    for(auto &strip : Map){
      if(strip.first == MaxQ.first - 1){
        q1 = strip.second;
      }
      else if(strip.first == MaxQ.first + 1){
        q3 = strip.second;
      }
    }
    //std::cout << "test q " << q1 << " " << q2 << " " << q3 << std::endl;
    double vs[6];
    // std::cout << "test AnH 3" << std::endl;
    if(q1 > 0 && q3 > 0)    
    {
      // QsumSample[DetNum].push_back(QSum);
      vs[0] = sqrt(q2/q3);
      vs[1] = sqrt(q2/q1);
      vs[2] = 0.5*(vs[0] + vs[1]);
      vs[3] = log( vs[2] + sqrt(vs[2]*vs[2]-1.0) );
      vs[4] = abs((vs[0] - vs[1])/(2.0*sinh(vs[3])));	
      vs[5] = 0.5*log( (1.0+vs[4])/(1.0-vs[4]) ) ;
      // std::cout << "test AnH 4" << std::endl;

      if ( q3>q1 ) 
        sech = MaxQ.first + vs[5]/vs[3] ;
      else 
        sech = MaxQ.first - vs[5]/vs[3] ;
      //std::cout << "test sech " << sech << std::endl;
    }
  }

  return sech ;
}


///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysizeX = m_EventData->GetFPMWMultX();
  for (unsigned int i = 0; i < mysizeX ; ++i) {
    int det = m_EventData->GetFPMW_DetX(i);
    int strip = m_EventData->GetFPMW_StripX(i);
    double QX = m_EventData->GetFPMW_ChargeX(i);
    double Qcal = Cal->ApplyCalibration("FPMW/DET"+NPL::itoa(det)+"_STRIP"+NPL::itoa(strip),QX);
    m_PreTreatedData->SetFPMW_DetX(det);
    m_PreTreatedData->SetFPMW_StripX(strip);
    m_PreTreatedData->SetFPMW_ChargeX(Qcal);
  }

  unsigned int mysizeY = m_EventData->GetFPMWMultY();
  for (unsigned int i = 0; i < mysizeY ; ++i) {
    int det = m_EventData->GetFPMW_DetY(i);
    int strip = m_EventData->GetFPMW_StripY(i);
    double QY = m_EventData->GetFPMW_ChargeY(i);
    double Qcal = Cal->ApplyCalibration("FPMW/DET"+NPL::itoa(det)+"_STRIP"+NPL::itoa(strip),QY);
    m_PreTreatedData->SetFPMW_DetY(det);
    m_PreTreatedData->SetFPMW_StripY(strip);
    m_PreTreatedData->SetFPMW_ChargeY(Qcal);
  }

}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigFPMW.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigFPMW.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigFPMW.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigFPMW.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigFPMW";
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
void TFPMWPhysics::Clear() {
  DetectorNbr.clear();
  PositionX.clear();
  PositionY.clear();
  ChargeX.clear();
  ChargeY.clear();

  DetectorHitX.clear();
  DetectorHit.clear();

  MapX.clear();
  MapY.clear();
  MaxQX.clear();
  MaxQY.clear();
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("FPMW");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FPMW " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FPMW " << i+1 <<  endl;
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
void TFPMWPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int i = 0; i < m_NumberOfDetectors; i++){
    for(int s = 0; s < 996; s++){
      Cal->AddParameter("FPMW","DET"+NPL::itoa(i)+"_STRIP"+NPL::itoa(s),"FPMW_DET"+NPL::itoa(i)+"_STRIP"+NPL::itoa(s));
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("FPMW",  true );
  inputChain->SetBranchAddress("FPMW", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("FPMW", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("FPMW", "TFPMWPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TFPMWPhysics::Construct() {
  return (NPL::VDetector*) new TFPMWPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_FPMW{
    public:
      proxy_FPMW(){
        NPL::DetectorFactory::getInstance()->AddToken("FPMW","FPMW");
        NPL::DetectorFactory::getInstance()->AddDetector("FPMW",TFPMWPhysics::Construct);
      }
  };

  proxy_FPMW p_FPMW;
}

