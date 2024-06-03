#include "IsoldeMerging3.hh"



void readConfig(const std::string& filename) {
    std::ifstream configFile(filename);
    if (!configFile.is_open()) {
        std::cerr << "Unable to open config file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(configFile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '=')) {
            std::string value;
            if (std::getline(iss, value)) {
                if (key == "t_max") {
                    t_max = std::stod(value);
                } else if (key == "coinc_time") {
                    coinc_time = std::stod(value);
                } else if (key == "T_ADC") {
                    T_ADC = std::stod(value);
                } else if (key == "FileName") {
                    FileName = value;
                }
            }
        }
    }
    configFile.close();
}

void CreateFileNames(const std::string& FileName, std::string& InputFileName,
                     std::string& TstampFileName, std::string& MusettFileName)
{
    // Construct input file name
    InputFileName = FileName + ".root";

    // Construct output file name with Timestamp
    std::ostringstream TStampFileNameStream;
    TStampFileNameStream << FileName << "_TSTAMP_T" << t_max << ".root";
    TstampFileName = TStampFileNameStream.str();

    // Construct output file name with MusettFileName
    std::ostringstream MusettFileNameStream;
    MusettFileNameStream << FileName << "_MUSETT_tau" << coinc_time << "_T" << t_max << "_TADC" << T_ADC << ".root";
    MusettFileName = MusettFileNameStream.str();
}


void CreateTempFile(const std::string& inputFileName)
{
  TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Failed to open input file!" << std::endl;
      return;
  }
  TTree *tree = static_cast<TTree*>(inputFile->Get("SimulatedTree"));
  if (!tree) {
      std::cerr << "Failed to get tree from input file!" << std::endl;
      inputFile->Close();
      return;
  }
  tree->SetBranchStatus("*", 0);  // Disable all branches initially
  tree->SetBranchStatus("MUSETT", 1);  // Enable the MUSETT branch only

  TMUSETTData* Mumu = nullptr;
  tree->SetBranchAddress("MUSETT", &Mumu);

  // Create a temporary file to store unsorted data
  TFile* tempFile = new TFile("temp.root", "RECREATE");
  TTree* tempTree = new TTree("EventTree", "Analysis Data");

  MusettSimuEvent* event = new MusettSimuEvent();
  tempTree->Branch("Event", &event);

  Long64_t nEntries = tree->GetEntries();
  int multX, multY;
  for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      multX = Mumu->GetDSSDXEMult();
      multY = Mumu->GetDSSDYEMult();

      for (int j = 0; j < multX; j++) {
          event->Clear();
          event->FillEvent_X(Mumu, j);
          tempTree->Fill();
      }

      for (int j = 0; j < multY; j++) {
          event->Clear();
          event->FillEvent_Y(Mumu, j);
          tempTree->Fill();
      }
  }
  tempTree->Write();
  tempFile->Close();
  inputFile->Close();
}


void ApplyTimeStamp(const std::string& inputFileName, const std::string& TstampFileName) {
    // Open the input file and get the simulated tree
    if (std::filesystem::exists(TstampFileName)) {
      std::cout << "Output file " << TstampFileName << " already exists. Skipping TStamp." << std::endl;
    return;
    }

    CreateTempFile(inputFileName);

    std::cout << "After temp file" << std::endl;

    TFile* readTempFile = new TFile("temp.root", "READ");
    TTree* tempTree = static_cast<TTree*>(readTempFile->Get("EventTree"));

    MusettSimuEvent* event = nullptr;
    tempTree->SetBranchAddress("Event", &event);

    Long64_t nEntries = tempTree->GetEntries();

    TFile* outputFile = new TFile(TstampFileName.c_str(), "RECREATE");
    TTree* outputTree = new TTree("SortedEventTree", "Events Sorted by Timestamp");
    outputFile->cd();

    MusettSimuEvent* sortedEvent = nullptr;
    outputTree->Branch("Event", &sortedEvent);


    std::vector<MusettSimuEvent*> events;
    events.reserve(nEntries);

    for (Long64_t i = 0; i < nEntries; ++i) {
        tempTree->GetEntry(i);
        events.push_back(new MusettSimuEvent(*event));
    }
    //readTempFile->Close();
    // Sort events by timestamp
    std::sort(events.begin(), events.end(), TimeComparator());

    // Open output file to store sorted data


    // Fill the tree with sorted events
    for (auto& evt : events) {
        sortedEvent = evt;
        outputTree->Fill();
    }
    outputTree->Write();
    outputFile->Close();


    // Clean up the copied events
    for (auto& evt : events) {
        delete evt;
    }

    // Remove the temporary file from disk
    readTempFile->Close();
    //std::remove("temp.root");
}


void InitializeVectors()
{
  SameTimeEvents.resize(2);  // Resize for two types: X and Y
  HitStrips.resize(2); // 0 == Y, 1 == X
  for (int i = 0; i < 2; ++i) { //Loop on X/Y
      SameTimeEvents[i].resize(NDet);
      for (int j = 0; j < NDet; ++j) {//Loop on Detectors
          SameTimeEvents[i][j].resize(NStrip);
          }
      }
}


void ClearVectors()
{
  for(int iX = 0; iX <2 ; iX++)
  {
    for (auto& pair : HitStrips[iX])
    {
      SameTimeEvents[iX][pair.first][pair.second]->Clear();
    }
  }

  HitStrips[0].clear();
  HitStrips[1].clear();



}

void TreatCoincEvent(MusettSimuEvent* event)
{
  int det = event->GetDet();
  int strip = event->GetStrip();
  std::pair<int, int> detStripPair = {det, strip};

  // Check if the detector and strip combination has not been hit
  if (HitStrips[event->IsX()].find(detStripPair) == HitStrips[event->IsX()].end())
  {
    HitStrips[event->IsX()].insert({det, strip});
    SameTimeEvents[event->IsX()][det][strip] = event;
  } else{
    double baseEnergy = SameTimeEvents[event->IsX()][det][strip]->GetEnergy();
    double baseTime = SameTimeEvents[event->IsX()][det][strip]->GetTotalTime();
    double sumEnergy = sumE2(baseEnergy,event->GetEnergy(), event->GetTotalTime()-baseTime);
    //std::cout << "*** det, strip = " << det << "," << strip ;

    SameTimeEvents[event->IsX()][det][strip]->SetEnergy(sumEnergy);
    //std::cout << ", sumEnergy =" << SameTimeEvents[event->IsX()][det][strip]->GetEnergy() << std::endl;
  }
}


void CorrelateEvents(const std::string& TstampFileName, const std::string& MusettFileName )
{
  // OOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOO
  //            Open input and output file
  // OOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOO

  TFile* inputFile = new TFile(TstampFileName.c_str(), "READ");
  TTree* inputTree = (TTree*)inputFile->Get("SortedEventTree");

  MusettSimuEvent* inputEvent = new MusettSimuEvent;

  inputTree->SetBranchAddress("Event", &inputEvent);


  TFile* outputFile = new TFile(MusettFileName.c_str(), "RECREATE");
  outputFile->cd();

  TTree* outputTree = new TTree("SimulatedTree", "Correlated Data");
  TMUSETTData* Mumu = new TMUSETTData();
  outputTree->Branch("MUSETT",&Mumu);


  // OOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOO
  //            Set Variables for the loop
  // OOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOOOOO00Ooo00ooO0OOO

  const int nEntries = (const int) inputTree->GetEntries();
  //const int nEntries = (const int) 1e6;
  inputTree->GetEntry(0);
  double lastTime = inputEvent->GetTotalTime();

  std::vector<bool> treated(nEntries, false);  // Initialize all entries to false
  std::pair<int, int> det_strip_pair;
  int entryNumber;
  int det, strip;
  double baseEnergy, baseTime, dt, totalEnergy;
  bool checkHit ;
  MusettSimuEvent* clonedEvent = nullptr;
  MusettSimuEvent* CorrEvent = nullptr;

  int Ncoinc = 0;
  int Nx, Ny;
  Nx = 0;
  Ny = 0;
  for (int i = 0; i < nEntries; i++) {
    //std::cout << std::endl;
    std::cout << "\r" // Move the cursor to the beginning of the line
          << "i/nEntries = " << i << " / " << nEntries
          << std::setw(20) << std::setfill(' ') // Fill the rest of the line if needed
          << std::flush; // Flush the stream to ensure output is updated immediately
    inputTree->GetEntry(i);
    if (treated[i])
    {
      //std::cout <<  "-> already treated, continue" << std::endl;
      continue;
    }
    //std::cout << "T = " << inputEvent->GetTotalTime() << std::endl;
    //std::cout << "∆T = " << inputEvent->GetTotalTime()-lastTime << std::endl;
    //std::cout << "(det, strip ) = (" << inputEvent->GetDet() << "," << inputEvent->GetStrip()  << "), ";
    if (inputEvent->GetTotalTime() - lastTime < coinc_time) {
      //std::cout << "-> in coinc " << std::endl;
      clonedEvent = new MusettSimuEvent(*inputEvent);
      TreatCoincEvent(clonedEvent);
      treated[i] = true;
      Ncoinc++;
    } else {
      //std::cout << "😡 Ncoinc = " << Ncoinc << std::endl;
      Ncoinc = 0;
      lastTime = inputEvent->GetTotalTime();

      for(int j = i; j < nEntries; ++j)
      {
        inputTree->GetEntry(j);
        if(inputEvent->GetTotalTime()-lastTime  > T_ADC) break;
        det = inputEvent->GetDet();
        strip = inputEvent->GetStrip();
        det_strip_pair = {det, strip};
        checkHit = HitStrips[inputEvent->IsX()].find(det_strip_pair) != HitStrips[inputEvent->IsX()].end() ;
        if(checkHit && !treated[j])
        {
          dt = inputEvent->GetTotalTime()-SameTimeEvents[inputEvent->IsX()][det][strip]->GetTotalTime();
          if(dt < T_ADC)
          {
            //std::cout << std::endl;
            if (inputEvent->IsX() == 1) {
                //std::cout << "XXXXX";
                Nx++;
            } else {
                //std::cout << "YYYYY";
                Ny++;
            }
            //std::cout << ", dt =" << dt ;

            totalEnergy = SameTimeEvents[inputEvent->IsX()][det][strip]->GetEnergy();
            //std::cout << ", E1 = " << totalEnergy << ", E2 =" << inputEvent->GetEnergy();
            totalEnergy = sumE2(totalEnergy,inputEvent->GetEnergy(),dt);

            //std::cout << ", total energy = " << totalEnergy << std::endl;
            //std::cout << std::endl;
            SameTimeEvents[inputEvent->IsX()][det][strip]->SetEnergy(totalEnergy);
            treated[j] = true;
          }
        }
        //TreatADCSum(inputEvent,treated,j);
      }
      for (auto& Xpair : HitStrips[1])
      {
        CorrEvent = SameTimeEvents[1][Xpair.first][Xpair.second];
        Mumu->SetDSSDXT(false, CorrEvent->GetDet(), CorrEvent->GetStrip(), CorrEvent->GetTime());
        Mumu->SetDSSDXE(false, CorrEvent->GetDet(), CorrEvent->GetStrip(), CorrEvent->GetEnergy());
        Mumu->SetDSSDX_Interstrip(CorrEvent->IsInterstrip());
        Mumu->SetDSSDX_Tstamp(CorrEvent->GetTimeStamp());

      }
      for (auto& Ypair : HitStrips[0])
      {
        //std::cout << "EY = " << CorrEvent->GetEnergy() << std::endl;
        CorrEvent = SameTimeEvents[0][Ypair.first][Ypair.second];
        //std::cout << "EY = " << CorrEvent->GetEnergy() << std::endl;
        Mumu->SetDSSDYT(false, CorrEvent->GetDet(), CorrEvent->GetStrip(), CorrEvent->GetTime());
        Mumu->SetDSSDYE(false, CorrEvent->GetDet(), CorrEvent->GetStrip(), CorrEvent->GetEnergy());
        Mumu->SetDSSDY_Interstrip(CorrEvent->IsInterstrip());
        Mumu->SetDSSDY_Tstamp(CorrEvent->GetTimeStamp());
      }

      outputTree->Fill();
      Mumu->Clear();
      ClearVectors();

      i -= 1; // Decrement to revisit this entry on next loop
    }
}
  outputFile->Write();
  inputFile ->Close();
  outputFile->Close();

}


void IsoldeMerging3()
{
  std::string InputFileName, TstampFileName, MusettDataFileName;
  CreateFileNames(FileName, InputFileName, TstampFileName, MusettDataFileName);
  ApplyTimeStamp(InputFileName, TstampFileName);
  InitializeVectors();
  CorrelateEvents(TstampFileName, MusettDataFileName);
}
