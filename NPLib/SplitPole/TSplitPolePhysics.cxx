/*****************************************************************************
 * Copyright (C) 2009-2013   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@ipno.in2p3.fr       *
 *                                                                           *
 * Creation Date  : november 2009                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *    This class hold SplitPole  Physics                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
//  NPL
#include "TSplitPolePhysics.h"

// C++ headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <stdlib.h>
using namespace std;

// NPTool headers
#include "RootOutput.h"
#include "RootInput.h"
#include "NPDetectorFactory.h"

//  ROOT
#include "TChain.h"

ClassImp(TSplitPolePhysics)



///////////////////////////////////////////////////////////////////////////
TSplitPolePhysics::TSplitPolePhysics()
   : m_EventData(new TSplitPoleData),
     m_PreTreatedData(new TSplitPoleData),
     m_EventPhysics(this),
     m_MagneticFieldCorrection(0),
     m_TimeDelay(6500),
     m_CalibP0(0.65),
     m_CalibP1(6e-6)
{    
}



///////////////////////////////////////////////////////////////////////////
TSplitPolePhysics::~TSplitPolePhysics()
{
   delete m_EventData;
   delete m_PreTreatedData;
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::Clear()
{
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::ReadConfiguration(string Path) 
{
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());
   string LineBuffer, DataBuffer;

   bool ReadingStatus = false;

   cout << "SP read configuration" << endl;

   while (!ConfigFile.eof()) {      
      getline(ConfigFile, LineBuffer);
      cout << LineBuffer << endl;

      // If SplitPole detector found, toggle Reading Block Status
      if (LineBuffer.compare(0, 9, "SplitPole") == 0) {
         cout << "Detector found: " << endl;
         ReadingStatus = true;
      }
      // else don't toggle to Reading Block Status
      else ReadingStatus = false;

      // Reading Block
      while (ReadingStatus) {
         // Pickup Next Word 
         ConfigFile >> DataBuffer;

         //  Comment Line 
         if (DataBuffer.compare(0, 1, "%") == 0) {ConfigFile.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );}

         //  Finding another telescope (safety), toggle out
         else if (DataBuffer.compare(0, 9, "SplitPole") == 0) {
            cout << "WARNING: Another Telescope is find before standard sequence of Token, Error may occured in Telecope definition" << endl;
            ReadingStatus = false;
         }

         //  If no Detector Token and no comment, toggle out
         else {
            ReadingStatus = false; 
            cout << "Wrong Token Sequence: Getting out " << DataBuffer << endl;
         }
      }
   }
          
   InitializeStandardParameters();
   ReadAnalysisConfig();
}


///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::AddParameterToCalibrationManager()
{
}



///////////////////////////////////////////////////////////////////////////
void  TSplitPolePhysics::InitializeRootInputRaw() 
{
   TChain* inputChain = RootInput::getInstance()->GetChain();
   inputChain->SetBranchStatus("SplitPole"   , true);
   inputChain->SetBranchStatus("fSplitPole_*", true);
   inputChain->SetBranchAddress("SplitPole"  , &m_EventData);
}



///////////////////////////////////////////////////////////////////////////
void  TSplitPolePhysics::InitializeRootInputPhysics() 
{
   TChain* inputChain = RootInput::getInstance()->GetChain();
   inputChain->SetBranchStatus("SplitPole",  true);
   inputChain->SetBranchStatus("fPosition",  true);
   inputChain->SetBranchStatus("fBrho",      true);
   inputChain->SetBranchStatus("fDeltaE",    true);
   inputChain->SetBranchStatus("fWire",      true);
   inputChain->SetBranchStatus("fPlasticP",  true);
   inputChain->SetBranchStatus("fPlasticG",  true);
   inputChain->SetBranchAddress("SplitPole", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::InitializeRootOutput()
{
   TTree* outputTree = RootOutput::getInstance()->GetTree();
   outputTree->Branch("SplitPole", "TSplitPolePhysics", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::BuildPhysicalEvent()
{
   BuildSimplePhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::BuildSimplePhysicalEvent()
{
   // Select active channels and apply thresholds
   PreTreat();

   // Correct for magnetic field variation
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::PreTreat()
{
   // Clear pre treated object
   ClearPreTreatedData();

   // pre treated object is the same as the raw data object
}



/*
bool TSplitPolePhysics::IsValidChannel(string Type, int detector, int channel)
{
   vector<bool>::iterator it;
   if (Type == "Front")
      return *(m_FrontChannelStatus[detector].begin()+channel);

   else if (Type == "Back")
      return *(m_BackChannelStatus[detector].begin()+channel);

   else 
      return false;
}
*/


///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::InitializeStandardParameters()
{
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::ReadAnalysisConfig()
{
   bool ReadingStatus = false;

   cout << "\t/////////// Reading ConfigSplitPole.dat file ///////////" << endl;

   // path to file
   string FileName = "./configs/ConfigSplitPole.dat";

   // open analysis config file
   ifstream AnalysisConfigFile;
   AnalysisConfigFile.open(FileName.c_str());

   if (!AnalysisConfigFile.is_open()) {
      cout << "\tNo ConfigSplitPole.dat found: Default parameter loaded for Analayis " << FileName << endl;
      return;
   }
   cout << "\tLoading user parameters from ConfigSplitPole.dat " << endl;

   // storing config file in the ROOT output file
   TAsciiFile *asciiFile = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
   asciiFile->AppendLine("%% ConfigSplitPole.dat %%");
   asciiFile->Append(FileName.c_str());
   asciiFile->AppendLine("");

   // read analysis config file
   string LineBuffer, DataBuffer, whatToDo;
   while (!AnalysisConfigFile.eof()) {
      // Pick-up next line
      getline(AnalysisConfigFile, LineBuffer);

      // search for "header"
      if (LineBuffer.compare(0, 15, "ConfigSplitPole") == 0) ReadingStatus = true;

      // loop on tokens and data
      while (ReadingStatus) {
         whatToDo = "";
         AnalysisConfigFile >> whatToDo;

         // Search for comment symbol (%)
         if (whatToDo.compare(0, 1, "%") == 0) {
            AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
         }

         else if (whatToDo == "CALIBRATION_P0") {
            AnalysisConfigFile >> DataBuffer;
            m_CalibP0 = atoi(DataBuffer.c_str());
            cout << "\t" << whatToDo << "\t" << m_CalibP0 << endl;
         }

         else if (whatToDo == "CALIBRATION_P1") {
            AnalysisConfigFile >> DataBuffer;
            m_CalibP1 = atof(DataBuffer.c_str());
            cout << "\t" << whatToDo << "\t" << m_CalibP1 << endl;
         }

         else if (whatToDo == "MAGNETIC_FIELD_CORRECTION") {
            AnalysisConfigFile >> DataBuffer;
            m_MagneticFieldCorrection = false;
            if (DataBuffer == "ON") m_MagneticFieldCorrection = true;
            cout << "\t" << whatToDo << "\t" << DataBuffer << endl;
         }

         else if (whatToDo == "TIME_DELAY") {
            AnalysisConfigFile >> DataBuffer;
            m_TimeDelay = atof(DataBuffer.c_str());
            cout << "\t" << whatToDo << "\t" << m_TimeDelay << endl;
         }

         else {
            ReadingStatus = false;
         }
      }
   }
   cout << "\t/////////////////////////////////////////////////" << endl;
}   



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::InitSpectra()
{
//   m_Spectra = new TSplitPoleSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::FillSpectra()
{
//   m_Spectra->FillRawSpectra(m_EventData);
//   m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
//   m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::CheckSpectra()
{
//   m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::ClearSpectra()
{
   // To be done
}



///////////////////////////////////////////////////////////////////////////
void TSplitPolePhysics::WriteSpectra()
{
//   m_Spectra->WriteSpectra();
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSplitPolePhysics::Construct()
{
   return (NPL::VDetector*) new TSplitPolePhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_splitpole {
   public:
      proxy_splitpole() {
         NPL::DetectorFactory::getInstance()->AddToken("SplitPole", "SplitPole");
         NPL::DetectorFactory::getInstance()->AddDetector("SplitPole", TSplitPolePhysics::Construct);
      }
};

proxy_splitpole p;
}

