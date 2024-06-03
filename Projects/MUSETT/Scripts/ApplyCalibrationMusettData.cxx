#include "ApplyCalibrationMusettData.hh"

namespace CALIB {
  //   DSSD
  //   X
  double fDSSD_X_E(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDXEDetectorNbr(i));
    //std::cout << "Det : " << m_EventData->GetDSSDXEDetectorNbr(i) << std::endl;
    name += "_DSSD_X";
    name += NPL::itoa(m_EventData->GetDSSDXEStripNbr(i));
    //std::cout << "Strip : " << m_EventData->GetDSSDXEStripNbr(i) << std::endl;
    name += "_E";
    //std::cout << name << std::endl;
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDXEEnergy(i),1);
  }

  double fDSSD_X_T(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDXTDetectorNbr(i));
    name += "_DSSD_X";
    name += NPL::itoa(m_EventData->GetDSSDXTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDXTTime(i),1);
  }

  //   Y
  double fDSSD_Y_E(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDYEDetectorNbr(i));
    name += "_DSSD_Y";
    name += NPL::itoa(m_EventData->GetDSSDYEStripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDYEEnergy(i),1);
  }

  double fDSSD_Y_T(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDYTDetectorNbr(i));
    name += "_DSSD_Y";
    name += NPL::itoa(m_EventData->GetDSSDYTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDYTTime(i),1);
  }
}

using namespace CALIB;

void ApplyCalibrationMusettData()
{
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  TTree *inputTree = static_cast<TTree*>(inputFile->Get("RD"));

  TMUSETTData* inputEvent = nullptr;
  inputTree->SetBranchAddress("MUSETT", &inputEvent);


  TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
  TTree *outputTree = new TTree("CalibMusettData", "Analysis Data");

  TMUSETTData* outputEvent = new TMUSETTData();
  outputTree->Branch("Event", &outputEvent);


  //const char* args[] = {"-C", "calibration/Calibration222Ra.txt"};
  const char* args[] = {" ", "-C", "calibration/Calibration222Ra.txt"};
  char** argv = const_cast<char**>(args);
  int argc = 3; // Since you have two arguments "-C" and the filename

  // Call getInstance with the prepared argc and argv
  NPOptionManager* myOptionManager = NPOptionManager::getInstance((int)2, argv);
  //NPOptionManager* myOptionManager = NPOptionManager::getInstance(argc,argv);

  myOptionManager->SetIsAnalysis();


  //std::cout << << std::endl;


  CalibrationManager* Cal = CalibrationManager::getInstance(myOptionManager->GetCalibrationFile());
  // Good for simulation, close to typical values
  vector<double> standardX    = {-63, 63. / 8192.};
  vector<double> standardY    = {63, -63. / 8192.};
  vector<double> standardSecondLayer  = {-63, 63. / 8192.};
  vector<double> standardT    = {-1000, 1000. / 8192.};


  //Cal->AddFile("/Users/lh270370/Software/nptool/Projects/MUSETT/calibration/Calibration222Ra.txt");

  for (int i = 0; i < 4; i++) {

    for (int j = 0; j < 128; j++) {
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_E",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_E",
          standardX);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_E",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_E",
          standardY);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_T",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_T",
          standardT);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_T",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_T",
          standardT);
    }

  }

  Cal->LoadParameterFromFile();



  Long64_t nEntries = inputTree->GetEntries();
  //Long64_t nEntries = 500;
  int multX, multY;
  double EX, EY;
  int det, strip;
  for (Long64_t i = 0; i < nEntries; ++i) {
    std::cout << "\r" // Move the cursor to the beginning of the line
          << "i/nEntries = " << i << " / " << nEntries
          << std::setw(20) << std::setfill(' ') // Fill the rest of the line if needed
          << std::flush; // Flush the stream to ensure output is updated immediately
      inputTree->GetEntry(i);
      multX = inputEvent->GetDSSDXEMult();
      multY = inputEvent->GetDSSDYEMult();

      for (int j = 0; j < multX; j++) {
        //std::cout << "before : " << inputEvent->GetDSSDXEEnergy(j);
        EX = fDSSD_X_E(inputEvent, j);
        //std::cout << "after : " << EX << std::endl;
        //std:cout << ", after EX = " << EX << std::endl;
        det = inputEvent->GetDSSDXEDetectorNbr(j);
        strip = inputEvent->GetDSSDXEStripNbr(j);
        //std::cout << "in Looop, (det, strip) = (" << det << "," << strip << ")" << std::endl;

        outputEvent->SetDSSDXE(true,det, strip, EX);
      }

      for (int j = 0; j < multY; j++) {
        EY = fDSSD_Y_E(inputEvent, j);
        det = inputEvent->GetDSSDYEDetectorNbr(j);
        strip = inputEvent->GetDSSDYEStripNbr(j);
        outputEvent->SetDSSDYE(true,det, strip, EY);
      }
      outputTree->Fill();
      outputEvent->Clear();
    }

    outputFile->Write();
    //sumEX->Draw();
    outputFile->Close();
    inputFile->Close();
}
