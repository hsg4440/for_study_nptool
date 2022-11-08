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
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0)
      {
  ICcounter = 0;
  ACcounter = 0;
  GGcounter = 0;
  Entry_Exit_counter = 0;
  Plasticcounter = 0;
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TZDDPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  Match_E_T("IC");
  Match_E_T("Plastic");
  Treat_DC();
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::Treat_DC(){
  // Dont have anything to modify for the moment
  unsigned int mysizeDC = m_PreTreatedData->GetMultDrift();
  for(int i = 0; i < mysizeDC; i++){
    DC_DetectorNumber.push_back(m_PreTreatedData->GetDrift_DetectorNbr(i));
    DC_DriftTime.push_back(m_PreTreatedData->Get_DriftTime(i));
  }
}
///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::Match_E_T(std::string Detector){
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy(Detector);
  unsigned int mysizeT = m_PreTreatedData->GetMultTime(Detector);
  for (UShort_t e = 0; e < mysizeE ; e++) {
    for (UShort_t t = 0; t < mysizeT ; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(Detector, e) == m_PreTreatedData->GetT_DetectorNbr(Detector, t)) {
        if(Detector == "IC"){
          IC_DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(Detector, e));
          IC_Energy.push_back(m_PreTreatedData->Get_Energy(Detector, e));
          IC_Time.push_back(m_PreTreatedData->Get_Time(Detector, t));
        }
        else if(Detector == "Plastic"){
          Plastic_DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(Detector, e));
          Plastic_Energy.push_back(m_PreTreatedData->Get_Energy(Detector, e));
          Plastic_Time.push_back(m_PreTreatedData->Get_Time(Detector, t));
        }
        else{
          std::cout << "Detector should be IC or Plastic" << std::endl;
          return;
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // IC
  PreTreatEnergy("IC", Cal);
  PreTreatTime("IC", Cal);

  // Plastic
  PreTreatEnergy("Plastic", Cal);
  PreTreatTime("Plastic", Cal);

  //DC
  PreTreatTime("DC", Cal);
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreatEnergy(std::string Detector, CalibrationManager* Cal){
  unsigned int mysize = m_EventData->GetMultEnergy(Detector);
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
}

///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::PreTreatTime(std::string Detector, CalibrationManager* Cal){
  unsigned int mysize = m_EventData->GetMultTime(Detector);
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("ZDD/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(Detector, i)),m_EventData->Get_Time(Detector,i));
    if(Detector == "IC")
      m_PreTreatedData->Set_IC_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
    else if(Detector == "Plastic")
      m_PreTreatedData->Set_Plastic_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
    else if(Detector == "DC")
      m_PreTreatedData->Set_DC_Time(m_EventData->GetT_DetectorNbr(Detector, i), Time);
  }
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
void TZDDPhysics::Clear() {
  IC_DetectorNumber.clear();
  IC_Energy.clear();
  IC_Time.clear();
  Plastic_DetectorNumber.clear();
  Plastic_Energy.clear();
  Plastic_Time.clear();
  DC_DetectorNumber.clear();
  DC_DriftTime.clear();
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ZDD");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> TokenZDD  = {"R", "Theta"};
  vector<string> TokenDC = {"Z", "Gas","Thickness", "Pressure", "Temperature"};
  vector<string> TokenAC = {"Z", "Thickness", "Material"};
  vector<string> TokenEntryExit = {"Z", "Thickness", "Material"};
  vector<string> TokenIC = {"Z", "Thickness", "Gas", "Pressure", "Temperature"};
  vector<string> TokenPlastic = {"Material", "Width", "Length", "Thickness", "Pos"};



  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if (blocks[i]->HasTokenList(TokenZDD)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "////  ZDD " << i + 1 << endl;
        double R     = blocks[i]->GetDouble("R", "mm");
        double Theta = blocks[i]->GetDouble("Theta", "deg");
        Add_ZDD(R, Theta);
    }
    else if (blocks[i]->GetMainValue() == "DC"
            && blocks[i]->HasTokenList(TokenDC)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// DC" << i + 1 << endl;
        double Z           = blocks[i]->GetDouble("Z", "mm");
        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
        string   Gas         = blocks[i]->GetString("Gas");
        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
        Add_Drift_Chamber(Z, Thickness, Gas, Pressure, Temperature);
        }
    else if (blocks[i]->GetMainValue() == "IC"
            && blocks[i]->HasTokenList(TokenIC)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// IC" << ICcounter+1 << endl;
        double Z           = blocks[i]->GetDouble("Z", "mm");
        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
        string   Gas         = blocks[i]->GetString("Gas");
        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
        Add_Ionisation_Chamber(Z, Thickness, Gas, Pressure, Temperature);
        ICcounter++;
    }
    else if (blocks[i]->GetMainValue() == "GasGap"
            && blocks[i]->HasTokenList(TokenIC)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// GasGap" << GGcounter+1 << endl;
        double Z           = blocks[i]->GetDouble("Z", "mm");
        double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
        string   Gas         = blocks[i]->GetString("Gas");
        double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
        double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
        Add_Gas_Gap(Z, Thickness, Gas, Pressure, Temperature);
        GGcounter++;
    }
    else if (blocks[i]->GetMainValue() == "AC"
            && blocks[i]->HasTokenList(TokenAC)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// AC" << ACcounter+1 << endl;
        double Z           = blocks[i]->GetDouble("Z", "mm");
        double Thickness   = blocks[i]->GetDouble("Thickness", "um");
        string   Material    = blocks[i]->GetString("Material");
        Add_AC(Z, Thickness, Material);
        ACcounter++;
    }
    else if (blocks[i]->GetMainValue() == "EntryExit"
            && blocks[i]->HasTokenList(TokenEntryExit)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// AC" << Entry_Exit_counter+1 << endl;
        double Z           = blocks[i]->GetDouble("Z", "mm");
        double Thickness   = blocks[i]->GetDouble("Thickness", "um");
        string   Material    = blocks[i]->GetString("Material");
        Add_Entry_Exit(Z, Thickness, Material);
        Entry_Exit_counter++;
    }
    else if (blocks[i]->GetMainValue() == "Plastic"
            && blocks[i]->HasTokenList(TokenPlastic)) {
        if (NPOptionManager::getInstance()->GetVerboseLevel())
            cout << endl << "//// Plastic" << Plasticcounter + 1 << endl;
        string        Material  = blocks[i]->GetString("Material");
        double      Width     = blocks[i]->GetDouble("Width", "mm");
        double      Length    = blocks[i]->GetDouble("Length", "mm");
        double      Thickness = blocks[i]->GetDouble("Thickness", "mm");
        int x = (blocks[i]->GetTVector3("Pos", "mm")).X();
        int y = (blocks[i]->GetTVector3("Pos", "mm")).Y();
        int z = (blocks[i]->GetTVector3("Pos", "mm")).Z();
        std::vector<int> Pos{x, y, z}; 
        Add_Plastic(Material, Width, Length, Thickness, Pos);
        Plasticcounter++;
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  quicksort(&Entry_Exit_Z, &Entry_Exit_Thickness, &Entry_Exit_Material);
  quicksort(&AC_Z, &AC_Thickness, &AC_Material);
  quicksort(&m_Ionisation_Chamber_Z, &m_Ionisation_Chamber_Thickness, &m_Ionisation_Chamber_Gas, &m_Gas_Gap_Pressure, &m_Ionisation_Chamber_Temperature);
  quicksort(&m_Gas_Gap_Z, &m_Gas_Gap_Thickness, &m_Gas_Gap_Gas, &m_Gas_Gap_Pressure, &m_Gas_Gap_Temperature);
}

void TZDDPhysics::quicksort(std::vector<double>* Z, std::vector<double>* Thickness, std::vector<string>* Material){
  for(int i = 0; i < Z->size(); i++){
  }
  if(Z->size()>0){
    double minZ;
    int minindex;
    for(int i = 0; i < Z->size()-1; i++){
      minZ = (*Z)[i+1];
      minindex = i+1;
      for(int p = i; p < Z->size(); p++){
        if((*Z)[p] < minZ){
          minZ = (*Z)[p];
          minindex = p;
        }
      }
      double interZ = (*Z)[minindex];
      (*Z)[minindex] = (*Z)[i];
      (*Z)[i] = interZ; 
      
      double interThickness = (*Thickness)[minindex];
      (*Thickness)[minindex] = (*Thickness)[i];
      (*Thickness)[i] = interThickness; 
      
      string interMaterial = (*Material)[minindex];
      (*Material)[minindex] = (*Material)[i];
      (*Material)[i] = interMaterial; 
    }
  }
  for(int i = 0; i < Z->size(); i++){
  }
}


void TZDDPhysics::quicksort(std::vector<double>* Z, std::vector<double>* Thickness, std::vector<string>* Gas, std::vector<double>* Pressure, std::vector<double>* Temperature){  
  for(int i = 0; i < Z->size(); i++){
  }
  if(Z->size()>0){
    double minZ;
    int minindex;
    for(int i = 0; i < Z->size()-1; i++){
      minZ = (*Z)[i+1];
      minindex = i+1;
      for(int p = i; p < Z->size(); p++){
        if((*Z)[p] < minZ){
          minZ = (*Z)[p];
          minindex = p;
        }
      }
      double interZ = (*Z)[minindex];
      (*Z)[minindex] = (*Z)[i];
      (*Z)[i] = interZ; 
      
      double interThickness = (*Thickness)[minindex];
      (*Thickness)[minindex] = (*Thickness)[i];
      (*Thickness)[i] = interThickness; 
      
      string interGas = (*Gas)[minindex];
      (*Gas)[minindex] = (*Gas)[i];
      (*Gas)[i] = interGas; 

      double interPressure = (*Pressure)[minindex];
      (*Pressure)[minindex] = (*Pressure)[i];
      (*Pressure)[i] = interPressure; 

      double interTemperature = (*Temperature)[minindex];
      (*Temperature)[minindex] = (*Temperature)[i];
      (*Temperature)[i] = interTemperature; 
    }
  }
  for(int i = 0; i < Z->size(); i++){
  }
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
  inputChain->SetBranchStatus("ZDD",  true );
  inputChain->SetBranchAddress("ZDD", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("ZDD", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TZDDPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("ZDD", "TZDDPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TZDDPhysics::Construct() {
  return (NPL::VDetector*) new TZDDPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_ZDD{
  public:
    proxy_ZDD(){
      NPL::DetectorFactory::getInstance()->AddToken("ZDD","ZDD");
      NPL::DetectorFactory::getInstance()->AddDetector("ZDD",TZDDPhysics::Construct);
    }
};

proxy_ZDD p_ZDD;
}

