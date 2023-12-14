/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sandra GIRON  contact address: giron@ipno.in2p3.fr       *
 *                  Benjamin LE CROM		   lecrom@ipno.in2p3.fr              *
 * Creation Date  : march 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold exogam treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TExogamPhysics.h"
using namespace EXOGAM_LOCAL;

//	STL
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <functional>

//	NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPVDetector.h"
#include "RootInput.h"
#include "RootOutput.h"
//	ROOT
#include "TChain.h"

///////////////////////////////////////////////////////////////////////////
ClassImp(TExogamPhysics)
    ///////////////////////////////////////////////////////////////////////////
    TExogamPhysics::TExogamPhysics() {
  EventMultiplicity = 0;
  ECC_Multiplicity = 0;
  GOCCE_Multiplicity = 0;
  NumberOfHitClover = 0;
  NumberOfHitCristal = 0;
  m_Spectra = NULL;
  NumberOfClover = 0;
  m_EXO_E_RAW_Threshold = 0;
  m_EXO_E_Threshold = 0;
  m_EXO_EHG_RAW_Threshold = 0;
  m_EXO_TDC_RAW_Threshold = 0;
  m_EXO_TDC_RAW_Threshold = 0;

  m_PreTreatedData = new TExogamData;
  m_EventData = new TExogamData;
  m_EventPhysics = this;
  NumberOfClover = 0;
  CloverMult = 0;
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::PreTreat() {
  // Clearing PreTreat TExogamData
  ClearPreTreatedData();
  // Clearing local variables for pretreat
  ResetPreTreatVariable();

  //E
  m_EXO_Mult = m_EventData->GetExoMult();

  for (unsigned int i = 0; i < m_EXO_Mult; ++i) {
    
    if (m_EventData->GetExoE(i) > m_EXO_E_RAW_Threshold)
      EXO_E = fEXO_E(m_EventData, i);
    
    if (m_EventData->GetExoEHG(i) > m_EXO_EHG_RAW_Threshold)
      EXO_EHG = fEXO_EHG(m_EventData, i);
    
    if (m_EventData->GetExoTDC(i) > m_EXO_TDC_RAW_Threshold)
      EXO_TDC = fEXO_T(m_EventData, i);
    
    EXO_Outer1 = fEXO_Outer(m_EventData, i, 1);
    EXO_Outer2 = fEXO_Outer(m_EventData, i, 2);
    EXO_Outer3 = fEXO_Outer(m_EventData, i, 3);
    EXO_Outer4 = fEXO_Outer(m_EventData, i, 4);

    if(EXO_E > m_EXO_E_Threshold){
      m_PreTreatedData->SetExo(m_EventData->GetExoCrystal(i), EXO_E,
      EXO_EHG, m_EventData->GetExoTS(i), EXO_TDC, 
      m_EventData->GetExoBGO(i), m_EventData->GetExoCsI(i), EXO_Outer1,
      EXO_Outer2, EXO_Outer3, EXO_Outer4);
    } 
  }
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ResetPreTreatVariable(){
  EXO_E = -1000;
  EXO_EHG = -1000;
  EXO_TDC = -1000;
  EXO_Outer1 = -1000;
  EXO_Outer2 = -1000;
  EXO_Outer3 = -1000;
  EXO_Outer4 = -1000;
}

void TExogamPhysics::BuildPhysicalEvent() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
  //PreTreat();

  //for(unsigned int i = 0; i < m_PreTreatedData->GetExoMult(); i++){
  //  mean_free_path = ComputeMeanFreePath(m_PreTreatedData->GetExoE(i));
  //}
  std::cout << ComputeMeanFreePath(7000) << std::endl;
/*
  if(PreTreatedData -> GetECCEMult() != PreTreatedData -> GetECCTMult()) cout << PreTreatedData -> GetECCEMult() << " "
  <<  PreTreatedData -> GetECCTMult() << endl;


  for(unsigned int i = 0 ; i < PreTreatedData -> GetECCEMult(); i++) {

    // cout << i << " " << cristal_E << endl;
    // if(PreTreatedData->GetECCTTime(i) > 0)
      {
    ECC_E.push_back(PreTreatedData->GetECCEEnergy(i));
    ECC_T.push_back(PreTreatedData->GetECCTTime(i));
    ECC_CloverNumber.push_back(PreTreatedData->GetECCEClover(i));
    ECC_CristalNumber.push_back(PreTreatedData->GetECCECristal(i));

    //    cout << "BuildPhys " << PreTreatedData->GetECCEClover(i) << " " <<  PreTreatedData->GetECCECristal(i)<< " " <<
  PreTreatedData->GetECCTTime(i) << " " << endl;
      }
  }


  for(unsigned int j = 0 ; j < PreTreatedData -> GetGOCCEEMult(); j++) {
    GOCCE_E.push_back(PreTreatedData->GetGOCCEEEnergy(j));
    GOCCE_CloverNumber.push_back(PreTreatedData->GetGOCCEEClover(j));
    GOCCE_CristalNumber.push_back(PreTreatedData->GetGOCCEECristal(j));
    GOCCE_SegmentNumber.push_back(PreTreatedData->GetGOCCEESegment(j));
  }


  //int NumberOfHitClover = 0;

  int DetectorID = -1;

  for( unsigned short i = 0 ; i < PreTreatedData->GetECCEMult() ; i++ )
    {
      // cout << PreTreatedData->GetECCEClover(i) << endl;
      if( PreTreatedData->GetECCEClover(i) != DetectorID)
  {
    if(i==0)
      {
        NumberOfHitClover++;
      }
    else if(PreTreatedData->GetECCEClover(i)!= PreTreatedData->GetECCEClover(i-1) )
      {
        NumberOfHitClover++;
      }
  }
      if(NumberOfHitClover == 4) break;
      //clover_mult -> Fill(NumberOfHitClover);

    }

  //cout << "NumberOfHitClover " << NumberOfHitClover << endl;

  map<int, vector<int> > MapCristal;
  map<int, vector<int> > MapSegment;

  map<int, vector<int> > :: iterator it;    // iterator used with MapCristal
  map<int, vector<int> > :: iterator at;    // iterator used with MapSegment

  vector<int> PositionOfCristal_Buffer_ECC;
  vector<int> PositionOfSegment_Buffer_GOCCE;


  //Fill map Cristal
  for(int clo = 0; clo < NumberOfClover; clo++)
    {
      for(unsigned int k = 0; k < ECC_CloverNumber.size(); k++)
  {
    if(ECC_CloverNumber.at(k) == clo) // && ECC_CristalNumber.at(k)== cri )
      PositionOfCristal_Buffer_ECC.push_back(k);
  }
      if(PositionOfCristal_Buffer_ECC.size() != 0) MapCristal[clo] = PositionOfCristal_Buffer_ECC;

      PositionOfCristal_Buffer_ECC.clear();

    }


  //Fill map Segment
  for(int clo = 0; clo < NumberOfClover; clo++)
    {
      for(int cri = 0; cri < 4 ; cri++)
  {
    //  for(int seg = 0; seg < 4 ; seg++)
      {
        for(unsigned int m = 0; m < GOCCE_CloverNumber.size(); m++)
    {
      if(GOCCE_CloverNumber.at(m) == clo && GOCCE_CristalNumber.at(m) == cri)// && GOCCE_SegmentNumber.at(m) == seg)
        {
          // PositionOfSegment_Buffer_GOCCE.push_back(4*clo+cri);
          PositionOfSegment_Buffer_GOCCE.push_back(m);
        }
    }
      }
      if(PositionOfSegment_Buffer_GOCCE.size() != 0) MapSegment[4*clo+cri] = PositionOfSegment_Buffer_GOCCE;

      PositionOfSegment_Buffer_GOCCE.clear();
  }
    }

  // Treatment
  for(int clo = 0; clo < NumberOfClover ; clo++)
    {
      double E = 0; double T = 0;
      int mult_cristal = 0;
      int cristal = -1 , segment;

      int cristal_Emax = 0; int cristal_Emin = 0;
      int Emax = 0, Emin = 1000000;
      int Tmin = 0, Tmax = 0;

      //ADD-BACK
      it = MapCristal.find(clo);

      int cristal_cond = 0;

      if(it != MapCristal.end())
  {
    vector<int> PositionOfCristal = it -> second;

    mult_cristal = PositionOfCristal.size();
    //if(mult_cristal!=0) cristal_mult -> Fill(mult_cristal);

    // ADD-BACK
    //cout << "boucle" << endl;

    for(unsigned int k = 0; k < PositionOfCristal.size(); k++)
      {
        int indice = PositionOfCristal.at(k);

        cristal_cond += ECC_CristalNumber.at(indice);
        // cout <<  ECC_CristalNumber.at(k) << " " ECC_E.at(k) << endl;

        if(mult_cristal < 3)
    {
      E+= ECC_E.at(indice);

      if(ECC_E.at(indice) < Emin) {
        cristal_Emin = ECC_CristalNumber.at(indice);
        Emin = ECC_E.at(indice);
        Tmin = ECC_T.at(indice);
      }

      if(ECC_E.at(indice) > Emax) {
        cristal_Emax = ECC_CristalNumber.at(indice);
        Emax = ECC_E.at(indice);
        Tmax = ECC_T.at(indice);
      }
    }

        else // case of multiplicity = 3 or 4
    {
      E = -1; cristal_Emax = -1; cristal_Emin = -1; Tmax = -1; Tmin = -1;
    }

        // cout << ECC_E.at(indice) << " " << Emax << " " << cristal_Emax << " " << Emin << " " << cristal_Emin << endl;

      }

    if( (mult_cristal==1) || (mult_cristal ==2  && cristal_cond %2 == 1) )
      {
        // cout << cristal_cond << endl;

        //cristal = cristal_Emax; T = Tmax;
        //cout << Emax << " " << cristal_Emax << " " << Emin << " " << cristal_Emin << endl;

        if(E > 500) { cristal = cristal_Emax; T = Tmax; }
        else        { cristal = cristal_Emin; T = Tmin; }


        // DOPPLER CORRECTION

        at = MapSegment.find(4*clo+cristal);
        segment = -1;

        if(at != MapSegment.end())
    {
      vector<int> PositionOfSegment = at -> second;     // position of segment k in the vector

      int segment_max = -1, E_temp = -1;

      for(unsigned int m = 0; m < PositionOfSegment.size(); m++)             // loop on hit segments of cristal cri of
  clover clo
        {
          int indice = PositionOfSegment.at(m);

          if(GOCCE_E.at(indice) > 0 && GOCCE_CristalNumber.at(indice) == cristal)
      {
        if( GOCCE_E.at(indice) > E_temp )
          {
            segment_max = GOCCE_SegmentNumber.at(indice) ;
            E_temp = GOCCE_E.at(indice);
          }
      }
        }
      segment = segment_max;
    }

      }


    if(E > 0 && cristal != -1 && segment != -1)
      {
        TotalEnergy_lab.push_back(E);
        Time.push_back(T);
        CloverNumber.push_back(clo);
        CristalNumber.push_back(cristal);
        SegmentNumber.push_back(segment);

        double theta = GetSegmentAngleTheta(clo, cristal, segment);

        Theta.push_back(theta);

        double doppler_E = DopplerCorrection(E, theta);
        DopplerCorrectedEnergy.push_back(doppler_E);

        //  cout << E  << " " << clo << " " << cristal << " " << segment << " " << theta << " " << doppler_E << endl;

      }

  }  // end of condition over CristalMap

    } // loop over NumberOfClover

  CloverMult = GetClover_Mult();

  //cout << "Exogam fine" << endl;
  */
}

double TExogamPhysics::ComputeMeanFreePath(double Energy){
  auto b = Map_PhotonCS.lower_bound(Energy);
  auto a = prev(b);
  if(b == Map_PhotonCS.begin()){
    a = b;
    b++;
  }
  else if(b == Map_PhotonCS.end()){
    b--;
    a = prev(b);
  }
  std::cout << a->first << " " << b->first << " " << a->second << " " << b->second << std::endl;
  double coeff = (Energy - a->first)/(b->first - a->first);

  double PhotonCrossSection = a->second + coeff*(b->second - a->second); // mm2/g
  return 1./(GeDensity*PhotonCrossSection);
}

double TExogamPhysics::DopplerCorrection(double E, double Theta) {
  double Pi = 3.141592654;
  TString filename = "configs/beta.txt";
  ifstream file;
  // cout << filename << endl;
  file.open(filename);
  if (!file)
    cout << filename << " was not opened" << endl;

  double E_corr = 0;
  double beta = 0.;
  file >> beta;
  double gamma = 1. / sqrt(1 - beta * beta);

  E_corr = gamma * E * (1. - beta * cos(Theta * Pi / 180.));

  return (E_corr);
}

///////////////////////////////////////////////////////////////////////////

void TExogamPhysics::Clear() {
  EventMultiplicity = 0;
  ECC_Multiplicity = 0;
  GOCCE_Multiplicity = 0;
  NumberOfHitClover = 0;
  NumberOfHitCristal = 0;

  ECC_CloverNumber.clear();
  ECC_CristalNumber.clear();
  GOCCE_CloverNumber.clear();
  GOCCE_CristalNumber.clear();
  GOCCE_SegmentNumber.clear();

  // ECC
  ECC_E.clear();
  ECC_T.clear();

  // GOCCE
  GOCCE_E.clear();

  CristalNumber.clear();
  SegmentNumber.clear();
  CloverNumber.clear();

  TotalEnergy_lab.clear();
  Time.clear();
  DopplerCorrectedEnergy.clear();
  Position.clear();
  Theta.clear();
}
///////////////////////////////////////////////////////////////////////////

////	Innherited from VDetector Class	////

//	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
void TExogamPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("EXOGAM");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " EXOGAM clover found " << endl;

  // FIXME ANGLE FILE??? NOT SURE I GET IT...
  // For Doppler I guess... Something like that should be added later
  // But maybe the more stand R,THETA,PHI or X,Y,Z 
  // vector<string> token = {"ANGLE_FILE"};
  vector<string> token = {"Board, Flange, Channel0, Channel1"};
  // FIXME To be implemented in the future
  // vector<string> token = {"Board, Flange, Channel0, Channel1, R, THETA, PHI"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(token)) {
      int Board, Flange, Channel0, Channel1; // FIXME!!!! Should come from Data...
      AddClover(Board, Flange, Channel0, Channel1);
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  ReadAnalysisConfig();
}

void TExogamPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;


  // path to photon cross section
  string CSFilename = "../../Inputs/PhotonCrossSection/CoherentGe.xcom";
  string LineBuffer;

  ifstream CSFile;
  CSFile.open(CSFilename.c_str());
  
  if (!CSFile.is_open()) {
    cout << " No CS file found "
         << CSFilename << endl;
    return;
  }
  while(CSFile.good()){
    double gammaE, CrossSection;
    getline(CSFile, LineBuffer);
    istringstream ss(LineBuffer);
    ss >> gammaE >> CrossSection; // E in MeV, converted to keV, CrossSection in cm2/g  
    gammaE *= 1000.; // Convertion to keV
    CrossSection *= 100.;
    Map_PhotonCS[gammaE] = CrossSection;
  }
  // path to file
  string FileName = "./configs/ConfigExogam.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigExogam.dat found: Default parameters loaded for "
            "Analysis "
         << FileName << endl;
    return;
  }
  
  
  string DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 11, "ConfigExogam") == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {

      whatToDo = "";
      AnalysisConfigFile >> whatToDo;
      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }

      else if (whatToDo == "EXO_Threshold") {
        //AnalysisConfigFile >> DataBuffer;
        //m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        //cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
      }

    }
  }
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitSpectra() { m_Spectra = new TExogamSpectra(NumberOfClover); }

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ClearSpectra() {
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TExogamPhysics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}


//////////////////////////////////////////////////////////////////////////
void TExogamPhysics::AddClover(int Board, int Flange, int Channel0, int Channel1) {

}

// FIXME Legacy thing... Might delete later
//////////////////////////////////////////////////////////////////////////
// void TExogamPhysics::AddClover(string AngleFile) {
//   ifstream file;
//   //  TString filename = Form("posBaptiste/angles_exogam_clover%d.txt",NumberOfClover);
//   //  TString filename = Form("posz42_simu50mm/angles_exogam_clover%d.txt",NumberOfClover);
//   //  TString filename = Form("posz42_exp_stat_demiring/angles_exogam_clover%d.txt",NumberOfClover);

//   string path = "configs/";
//   TString filename = path + AngleFile;

//   cout << filename << endl;
//   file.open(filename);
//   if (!file)
//     cout << filename << " was not opened" << endl;

//   vector<double> Angles;
//   vector<vector<double>> Segment_angles;
//   vector<vector<vector<double>>> Cristal_angles;

//   Cristal_angles.clear();

//   double angle;
//   string buffer;

//   for (int i = 0; i < 4; i++) {
//     Segment_angles.clear();

//     for (int j = 0; j < 4; j++) {
//       Angles.clear();

//       for (int k = 0; k < 2; k++) {
//         file >> buffer >> angle;

//         Angles.push_back(angle); // Theta (k = 0)   Phi (k = 1)

//         // cout << angle << endl;
//         if (Angles.size() == 2)
//           cout << "Clover " << NumberOfClover << ": Theta=" << Angles[0] << " Phi=" << Angles[1] << endl;
//       }

//       Segment_angles.push_back(Angles);
//     }

//     Cristal_angles.push_back(Segment_angles);
//   }

//   Clover_Angles_Theta_Phi.push_back(Cristal_angles);

//   file.close();

//   NumberOfClover++;
// }

//	Add Parameter to the CalibrationManger
void TExogamPhysics::AddParameterToCalibrationManager() {

  CalibrationManager* Cal = CalibrationManager::getInstance();

  for (int i = 0; i < NumberOfClover; i++) {
    for (int j = 0; j < 4; j++) {
      Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Elow",
                        "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Elow");
      Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Ehigh",
                        "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Ehigh");
      Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T",
                        "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T");

      for (int k = 0; k < 4; k++) {
        Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Seg" + NPL::itoa(k) + "_E",
                          "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_Seg" + NPL::itoa(k) + "_E");
      }
    }
  }
}

//	Activated associated Branches and link it to the private member DetectorData address
//	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
void TExogamPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("EXOGAM", true);
  inputChain->SetBranchStatus("fEXO_*", true);
  inputChain->SetBranchAddress("EXOGAM", &m_EventData);

  /*
  TList* outputList = RootOutput::getInstance()->GetList();
   clover_mult = new TH1F("clover_mult","clover_mult",20,0,20);
    outputList->Add(clover_mult);
  cristal_mult = new TH1F("cristal_mult","cristal_mult",20,0,20);
  outputList->Add(cristal_mult);
  */
  }
}

/////////////////////////////////////////////////////////////////////
//   Activated associated Branches and link it to the private member DetectorPhysics address
//   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
void TExogamPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("EventMultiplicty", true);
  inputChain->SetBranchStatus("ECC_Multiplicity", true);
  inputChain->SetBranchStatus("GOCCE_Multiplicity", true);
  inputChain->SetBranchStatus("ECC_CloverNumber", true);
  inputChain->SetBranchStatus("ECC_CristalNumber", true);
  inputChain->SetBranchStatus("GOCCE_CloverNumber", true);
  inputChain->SetBranchStatus("GOCCE_CristalNumber", true);
  inputChain->SetBranchStatus("GOCCE_SegmentNumber", true);
  inputChain->SetBranchStatus("ECC_E", true);
  inputChain->SetBranchStatus("ECC_T", true);
  inputChain->SetBranchStatus("GOCCE_E", true);
  inputChain->SetBranchStatus("CristalNumber", true);
  inputChain->SetBranchStatus("SegmentNumber", true);
  inputChain->SetBranchStatus("CloverNumber", true);
  inputChain->SetBranchStatus("CloverMult", true);
  inputChain->SetBranchStatus("TotalEnergy_lab", true);
  inputChain->SetBranchStatus("Time", true);
  inputChain->SetBranchStatus("DopplerCorrectedEnergy", true);
  inputChain->SetBranchStatus("Position", true);
  inputChain->SetBranchStatus("Theta", true);
  inputChain->SetBranchAddress("EXOGAM", &m_EventPhysics);
  }
}

/////////////////////////////////////////////////////////////////////

//	Create associated branches and associated private member DetectorPhysics address
void TExogamPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("EXOGAM", "TExogamPhysics", &m_EventPhysics);

  // control histograms if needed
  /*
  TList* outputList = RootOutput::getInstance()->GetList();
  controle = new TH1F("controle","histo de controle",20,0,20);
  outputList->Add(controle);
  */
}

void TExogamPhysics::SetTreeReader(TTreeReader* TreeReader) {
   TExogamPhysicsReader::r_SetTreeReader(TreeReader);
 }

///////////////////////////////////////////////////////////////////////////
namespace EXOGAM_LOCAL {
  //	tranform an integer to a string
  double fEXO_E(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoE(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoE(i));
  }
  
  double fEXO_EHG(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoEHG(i));
    name += "_EHG";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoEHG(i));
  }
  
  double fEXO_T(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoTDC(i));
    name += "_TDC";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoTDC(i));
  }
  
  double fEXO_Outer(const TExogamData* m_EventData, const unsigned int& i, const unsigned int OuterNumber) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoE(i));
    name += "_Outer";
    name += NPL::itoa(OuterNumber);
    name += "_E";
    if(OuterNumber == 1)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter1(i));
    else if(OuterNumber == 2)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter2(i));
    else if(OuterNumber == 3)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter3(i));
    else if(OuterNumber == 4)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter4(i));
    else{
      std::cout << "WARNING: Outer number != 1-4, something is wrong\n";
      return 0;
    };
  }
  
  string itoa(int value) {
    std::ostringstream o;

    if (!(o << value))
      return "";

    return o.str();
  }
} // namespace EXOGAM_LOCAL

/////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TExogamPhysics::Construct() { return (NPL::VDetector*)new TExogamPhysics(); }

NPL::VTreeReader* TExogamPhysics::ConstructReader() { return (NPL::VTreeReader*)new TExogamPhysicsReader(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_exogam {
 public:
  proxy_exogam() {
    NPL::DetectorFactory::getInstance()->AddToken("Exogam", "Exogam");
    NPL::DetectorFactory::getInstance()->AddDetector("Exogam", TExogamPhysics::Construct);
    NPL::DetectorFactory::getInstance()->AddDetectorReader("Exogam", TExogamPhysics::ConstructReader);
  }
};

proxy_exogam p;
}
