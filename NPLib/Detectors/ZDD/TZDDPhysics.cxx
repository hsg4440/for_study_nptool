/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TZDDPhysics.h"

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

ClassImp(TZDDPhysics)


///////////////////////////////////////////////////////////////////////////
TZDDPhysics::TZDDPhysics()
   : m_EventData(new TZDDData),
     m_PreTreatedData(new TZDDData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_IC_E_RAW_Threshold(0), // adc channels
     m_PL_E_RAW_Threshold(0), // adc channels
     m_DC_E_RAW_Threshold(0), // adc channels
     m_EXO_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0)
      {
  //ICcounter = 0;
  //ACcounter = 0;
  //GGcounter = 0;
  //Entry_Exit_counter = 0;
  //Plasticcounter = 0;
}

  
///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::BuildPhysicalEvent() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  Match_IC();
  if(IC_Nbr.size() > 0)
    Match_PL();
  // Treat_DC();
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::Treat_DC(){
 /* // Dont have anything to modify for the moment
  unsigned int mysizeDC = m_PreTreatedData->GetMultDrift();
  for(int i = 0; i < mysizeDC; i++){
    DC_DetectorNumber.push_back(m_PreTreatedData->GetDrift_DetectorNbr(i));
    DC_DriftTime.push_back(m_PreTreatedData->Get_DriftTime(i));
  }*/
}
///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::Match_IC(){
//////////////////////////////// Currently Simply matching if mult = 5, could be improved to treat mult > 5
  if(m_PreTreatedData->GetZDD_ICMult() == 5){
    // CHecking that each IC is only encountered once and then sorting them in the right order with the map
    for(unsigned int i = 0; i < m_PreTreatedData->GetZDD_ICMult(); i++){
    if(SortIC.find(m_PreTreatedData->GetZDD_ICN(i)) != SortIC.end()){
      SortIC.clear();
      break;
    }
    SortIC[m_PreTreatedData->GetZDD_ICN(i)] = std::make_pair(m_PreTreatedData->GetZDD_ICE(i),m_PreTreatedData->GetZDD_ICTS(i));
    }
    // Adding the IC info to the std::vectors
    ICSum = 0;
    for(auto it = SortIC.begin(); it != SortIC.end(); ++it){
    ICSum+= (it->second).first;
    IC_Nbr.push_back(it->first);
    IC_E.push_back((it->second).first);
    IC_TS.push_back((it->second).second);
    }
  }
}

void TZDDPhysics::Match_IC1(){
//////////////////////////////// Currently Simply matching if mult = 5, could be improved to treat mult > 5
    // CHecking that each IC is only encountered once and then sorting them in the right order with the map
    for(unsigned int i = 0; i < m_PreTreatedData->GetZDD_ICMult(); i++){
    if(SortIC.find(m_PreTreatedData->GetZDD_ICN(i)) == SortIC.end()){
        SortIC[m_PreTreatedData->GetZDD_ICN(i)] = std::make_pair(m_PreTreatedData->GetZDD_ICE(i),m_PreTreatedData->GetZDD_ICTS(i));
      }
    }
    // Adding the IC info to the std::vectors
    ICSum = 0;
    for(auto it = SortIC.begin(); it != SortIC.end(); ++it){
    ICSum+= (it->second).first;
    IC_Nbr.push_back(it->first);
    IC_E.push_back((it->second).first);
    IC_TS.push_back((it->second).second);
    }
}

void TZDDPhysics::Match_PL(){
  for(unsigned int i = 0; i < m_PreTreatedData->GetZDD_PLMult(); i++){
    SortPL[m_PreTreatedData->GetZDD_PLN(i)] = std::make_pair(m_PreTreatedData->GetZDD_PLE(i),m_PreTreatedData->GetZDD_PLTS(i));
  }
  for(auto it = SortPL.begin(); it != SortPL.end(); ++it){
  PL_Nbr.push_back(it->first);
  PL_E.push_back((it->second).first);
  PL_TS.push_back((it->second).second);
  }
}
///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();
  
  m_IC_Mult = m_EventData->GetZDD_ICMult();
  m_PL_Mult = m_EventData->GetZDD_PLMult();
  m_DC_Mult = m_EventData->GetZDD_DCMult();
  m_EXO_Mult = m_EventData->GetZDD_EXOMult();
  
  for (unsigned int i = 0; i < m_IC_Mult; ++i) {
    if (m_EventData->GetZDD_ICE(i) > m_IC_E_RAW_Threshold) {
        // std::cout << m_EventData->GetZDD_ICN(i) << " " << Map_IC[m_EventData->GetZDD_ICN(i)] << std::endl;
        // std::cout << m_IC_E_RAW_Threshold<< std::endl;
        m_PreTreatedData->SetZDDIC(Map_IC[m_EventData->GetZDD_ICN(i)], m_EventData->GetZDD_ICE(i), m_EventData->GetZDD_ICTS(i));
    }
  }
  
  for (unsigned int i = 0; i < m_PL_Mult; ++i) {
    if (m_EventData->GetZDD_PLE(i) > m_PL_E_RAW_Threshold) {
        m_PreTreatedData->SetZDDPL(m_EventData->GetZDD_PLN(i), m_EventData->GetZDD_PLE(i), m_EventData->GetZDD_PLTS(i));
    }
  }
  
  for (unsigned int i = 0; i < m_DC_Mult; ++i) {
    if (m_EventData->GetZDD_DCE(i) > m_DC_E_RAW_Threshold) {
        m_PreTreatedData->SetZDDDC(m_EventData->GetZDD_DCN(i), m_EventData->GetZDD_DCE(i), m_EventData->GetZDD_DCTS(i));
    }
  }
  
  for (unsigned int i = 0; i < m_EXO_Mult; ++i) {
    if (m_EventData->GetZDD_EXOE(i) > m_EXO_E_RAW_Threshold) {
        m_PreTreatedData->SetZDDEXO(m_EventData->GetZDD_EXON(i), m_EventData->GetZDD_EXOE(i), m_EventData->GetZDD_EXOTS(i));
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreatEnergy(std::string Detector, CalibrationManager* Cal){
  /*unsigned int mysize = m_EventData->GetMultEnergy(Detector);
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->Get_Energy(Detector, i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration("ZDD/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(Detector, i)),m_EventData->Get_Energy(Detector, i));
      if (Energy > m_E_Threshold) {
        if(Detector == "IC")
          m_PreTreatedData->Set_IC_Energy(m_EventData->GetE_DetectorNbr(Detector, i), Energy);

        else if(Detector == "Plastic")
          m_PreTreatedData->Set_Plastic_Energy(m_EventData->GetE_DetectorNbr(Detector, i), Energy);
      }
    }
  }
*/
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreatTime(std::string Detector, CalibrationManager* Cal){
/*  unsigned int mysize = m_EventData->GetMultTime(Detector);
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("ZDD/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(Detector, i)),m_EventData->Get_Time(Detector,i));
    if(Detector == "IC")
      m_PreTreatedData->Set_IC_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
    else if(Detector == "Plastic")
      m_PreTreatedData->Set_Plastic_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
    else if(Detector == "DC")
      m_PreTreatedData->Set_DC_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
  }
*/
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigZDD.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigZDD.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigZDD.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigZDD.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigZDD";
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

      else if (whatToDo=="IC_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_IC_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_IC_E_RAW_Threshold << endl;
      }
      else if (whatToDo=="PL_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_PL_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_PL_E_RAW_Threshold << endl;
      }
      else if (whatToDo=="DC_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_DC_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_DC_E_RAW_Threshold << endl;
      }
      else if (whatToDo=="EXO_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_EXO_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_EXO_E_RAW_Threshold << endl;
      }
      else if (whatToDo=="MAP_IC") {
        AnalysisConfigFile >> DataBuffer;
        Map_IC[atoi(DataBuffer.substr(0,1).c_str())] = atoi(DataBuffer.substr(1,1).c_str());
        cout << whatToDo << " " << atoi(DataBuffer.substr(0,1).c_str()) << " " << atoi(DataBuffer.substr(1,1).c_str()) << endl;
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
void TZDDPhysics::Clear() {
  ICSum = 0;
  SortIC.clear();
  SortPL.clear();
  IC_Nbr.clear();
  IC_E.clear();
  IC_TS.clear();
  PL_Nbr.clear();
  PL_E.clear();
  PL_TS.clear();
  DC_DetectorNumber.clear();
  DC_DriftTime.clear();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ZDD");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> TokenZDD  = {"R", "Theta"};
  //vector<string> TokenDC = {"Z", "Gas","Thickness", "Pressure", "Temperature"};
  //vector<string> TokenAC = {"Z", "Thickness", "Material"};
  //vector<string> TokenEntryExit = {"Z", "Thickness", "Material"};
  //vector<string> TokenIC = {"Z", "Thickness", "Gas", "Pressure", "Temperature"};
  //vector<string> TokenPlastic = {"Material", "Width", "Length", "Thickness", "Pos"};



  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if (blocks[i]->HasTokenList(TokenZDD)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "////  ZDD " << i + 1 << endl;
        double R     = blocks[i]->GetDouble("R", "mm");
        double Theta = blocks[i]->GetDouble("Theta", "deg");

        Add_ZDD(R, Theta);
    }
//    else if (blocks[i]->GetMainValue() == "DC"
//            && blocks[i]->HasTokenList(TokenDC)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// DC" << i + 1 << endl;
//        double Z           = blocks[i]->GetDouble("Z", "mm");
//        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
//        string   Gas         = blocks[i]->GetString("Gas");
//        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
//        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
//        Add_Drift_Chamber(Z, Thickness, Gas, Pressure, Temperature);
//        }
//    else if (blocks[i]->GetMainValue() == "IC"
//            && blocks[i]->HasTokenList(TokenIC)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// IC" << ICcounter+1 << endl;
//        double Z           = blocks[i]->GetDouble("Z", "mm");
//        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
//        string   Gas         = blocks[i]->GetString("Gas");
//        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
//        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
//        Add_Ionisation_Chamber(Z, Thickness, Gas, Pressure, Temperature);
//        ICcounter++;
//    }
//    else if (blocks[i]->GetMainValue() == "GasGap"
//            && blocks[i]->HasTokenList(TokenIC)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// GasGap" << GGcounter+1 << endl;
//        double Z           = blocks[i]->GetDouble("Z", "mm");
//        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
//        string   Gas         = blocks[i]->GetString("Gas");
//        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
//        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
//        Add_Gas_Gap(Z, Thickness, Gas, Pressure, Temperature);
//        GGcounter++;
//    }
//    else if (blocks[i]->GetMainValue() == "AC"
//            && blocks[i]->HasTokenList(TokenAC)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// AC" << ACcounter+1 << endl;
//        double Z           = blocks[i]->GetDouble("Z", "mm");
//        double Thickness   = blocks[i]->GetDouble("Thickness", "um");
//        string   Material    = blocks[i]->GetString("Material");
//        Add_AC(Z, Thickness, Material);
//        ACcounter++;
//    }
//    else if (blocks[i]->GetMainValue() == "EntryExit"
//            && blocks[i]->HasTokenList(TokenEntryExit)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// AC" << Entry_Exit_counter+1 << endl;
//        double Z           = blocks[i]->GetDouble("Z", "mm");
//        double Thickness   = blocks[i]->GetDouble("Thickness", "um");
//        string   Material    = blocks[i]->GetString("Material");
//        Add_Entry_Exit(Z, Thickness, Material);
//        Entry_Exit_counter++;
//    }
//    else if (blocks[i]->GetMainValue() == "Plastic"
//            && blocks[i]->HasTokenList(TokenPlastic)) {
//        if (NPOptionManager::getInstance()->GetVerboseLevel())
//            cout << endl << "//// Plastic" << Plasticcounter + 1 << endl;
//        string        Material  = blocks[i]->GetString("Material");
//        double      Width     = blocks[i]->GetDouble("Width", "mm");
//        double      Length    = blocks[i]->GetDouble("Length", "mm");
//        double      Thickness = blocks[i]->GetDouble("Thickness", "mm");
//        int x = (blocks[i]->GetTVector3("Pos", "mm")).X();
//        int y = (blocks[i]->GetTVector3("Pos", "mm")).Y();
//        int z = (blocks[i]->GetTVector3("Pos", "mm")).Z();
//        std::vector<int> Pos{x, y, z}; 
//        Add_Plastic(Material, Width, Length, Thickness, Pos);
//        Plasticcounter++;
//    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  ReadAnalysisConfig();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitSpectra() {
  m_Spectra = new TZDDSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TZDDPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("ZDD", "D"+ NPL::itoa(i+1)+"_ENERGY","ZDD_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("ZDD", "D"+ NPL::itoa(i+1)+"_TIME","ZDD_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  inputChain->SetBranchStatus("ZDD",  true );
  inputChain->SetBranchAddress("ZDD", &m_EventData );
  }
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  inputChain->SetBranchAddress("ZDD", &m_EventPhysics);
  }
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("ZDD", "TZDDPhysics", &m_EventPhysics);
}

void TZDDPhysics::SetTreeReader(TTreeReader* TreeReader) {
   TZDDPhysicsReader::r_SetTreeReader(TreeReader);
 }


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TZDDPhysics::Construct() {
  return (NPL::VDetector*) new TZDDPhysics();
}

NPL::VTreeReader* TZDDPhysics::ConstructReader() { return (NPL::VTreeReader*)new TZDDPhysicsReader(); }



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_ZDD{
  public:
    proxy_ZDD(){
      NPL::DetectorFactory::getInstance()->AddToken("ZDD","ZDD");
      NPL::DetectorFactory::getInstance()->AddDetector("ZDD",TZDDPhysics::Construct);
      NPL::DetectorFactory::getInstance()->AddDetectorReader("ZDD", TZDDPhysics::ConstructReader);
    }
};

proxy_ZDD p_ZDD;
}

