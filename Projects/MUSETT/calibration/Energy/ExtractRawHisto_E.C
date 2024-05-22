////////////////////////////////////////////////////////////////////////////
// This macro takes a converted ROOT file from the GANIL e569 exp. and    //
// create an histogram for each strip (X and Y) of the MUST2 array filled //
// with the energy. The histograms are dumped in an output ROOT file.     //
//                                                                        //
// This is especially usefull for calibration purposes when there is no   //
// need to work directly on the TChain                                     //
////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TMust2Data.h"
#include "TMugastData.h"
#include "TMUSETTData.h"
#include "RootInput.h"


#define NBTELESCOPE 4	
#define	NBSTRIPS	128
#define NBSILI     0


void ExtractRawHistos(const char* fname = "")
{
  RootInput* Input = RootInput::getInstance("RunToTreat.txt");
  TChain* Chain = Input->GetChain();

  // connect the TfDet.c_str()Data branch  
  Chain->SetBranchStatus("MUSETT",true);
  Chain->SetBranchStatus("fMM*",true);

  TMUSETTData* rawMUSETT = new TMUSETTData();
  Chain->SetBranchAddress("MUSETT",&rawMUSETT);

  // open the output ROOT file
  TString outFileName = "./Histograms/";
  outFileName += fname;
  outFileName += "_RawMUSETTHistos.root";
  std::cout<< outFileName<< std::endl;
  TFile *outFile = new TFile(outFileName, "recreate");

  // prepare output histograms for Must2
  TH1F* hStripXEnergy[NBTELESCOPE][NBSTRIPS];
  TH1F* hStripYEnergy[NBTELESCOPE][NBSTRIPS];
  //from i=4 because I wanted only the T5,T8
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSTRIPS; j++) {
      // strips XE
      TString hnameXE     = Form("hMM%d_STRX_E%d", i, j);
      TString htitleXE    = Form("MM%d_STRX_E%d", i, j);
      hStripXEnergy[i][j] = new TH1F(hnameXE, htitleXE, 16384, 0, 16384);
      // strips YE
      TString hnameYE     = Form("hMM%d_STRY_E%d", i, j);
      TString htitleYE    = Form("MM%d_STRY_E%d", i, j);
      hStripYEnergy[i][j] = new TH1F(hnameYE, htitleYE, 16384, 0, 16384);
    }
  }
  TH1F* hSiLiEnergy[NBTELESCOPE][NBSILI];
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSILI; j++) {
      TString hnameSiLiE     = Form("hMM%d_SILI_E%d", i, j);
      TString htitleSiLiE    = Form("MM%d_SILI_E%d", i, j);
      hSiLiEnergy[i][j] = new TH1F(hnameSiLiE, htitleSiLiE, 16384, 0, 16384);
    }
  }

  // read the data
  Int_t nentries = Chain->GetEntries();
  std::cout << "TChain contains " << nentries << " events" << std::endl;
  for (Int_t i = 0; i < nentries; i++) {
    if (i%1000 == 0) std::cout << "\rEntry " << i << "\t" << i/(double)nentries * 100 << " %" << flush;
    ///////////////////////
    // read Must2 branch //
    ///////////////////////
    //brMust2->GetEntry(i);
    Chain->GetEntry(i);
    // get StripXE multiplicity
    Int_t multXE = 0; 
    multXE = rawMUSETT->GetDSSDXEMult();
    //std::cout<< "multXE" << multXE<< std::endl;
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    //
    Int_t det = 0;   
    Int_t stripX = 0;
    Int_t energy = 0;

    for (Int_t j = 0; j < multXE; j++) {
      det    = rawMUSETT->GetDSSDXEDetectorNbr(j);
      stripX = rawMUSETT->GetDSSDXEStripNbr(j);
      energy = rawMUSETT->GetDSSDXEEnergy(j);
      if ((det < NBTELESCOPE) && (stripX < NBSTRIPS)) {
        hStripXEnergy[det][stripX]->Fill(energy);
      } 
      else {
        std::cout << "Error filling histograms: X_E" << std::endl;
        rawMUSETT->Dump();
      }
    }
    // get StripYE multiplicity
    Int_t multYE=0;
    multYE = rawMUSETT->GetDSSDYEMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms

    det = 0;   
    Int_t stripY = 0;
    energy = 0;
    for (Int_t j = 0; j < multYE; j++) {
      det    = rawMUSETT->GetDSSDYEDetectorNbr(j);
      stripY = rawMUSETT->GetDSSDYEStripNbr(j);
      energy = rawMUSETT->GetDSSDYEEnergy(j);

      if ((det < NBTELESCOPE) && (stripY < NBSTRIPS)) {
        hStripYEnergy[det][stripY]->Fill(energy);
      }
      else {
        std::cout << "Error filling histograms: Y_E" << std::endl;
          rawMUSETT->Dump();
      }
    }

  }
  std::cout << '\n';

  // write on disk output file and close it
  outFile->Write();
  outFile->Close();
}
