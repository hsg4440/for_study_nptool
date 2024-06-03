////////////////////////////////////////////////////////////////////////////
// This macro takes a converted ROOT file from the GANIL e569 exp. and    //
// create an histogram for each strip (X and Y) of the MUST2 array filled //
// with the energy. The histograms are dumped in an output ROOT file.     //
//                                                                        //
// This is especially usefull for calibration purposes when there is no   //
// need to work directly on the TTree                                     //
////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TMUSETTData.h"
#include "TChain.h"



void ExtractMust2Histos(const char* fname = "run_0040", int NBTELESCOPE = 4, int NBSTRIPS = 128, int FIRSTSTRIP = 0, int FIRSTELESCOPE = 0)
{

  TString path  = "~/npTreeReader/np-tree-reader/Projects/MUSETT/RootR/";
  TChain* tree = new TChain("RD");
  tree->Add(path + fname + ".root");
  //tree->Add(path + fname + "_1.root");
  //tree->Add(path + fname + "_2.root");

  tree->SetBranchStatus("*",false);

  // connect the TMust2Data branch  
  tree->SetBranchStatus("MUSETT",true);
  tree->SetBranchStatus("fMM*",true);
  TMUSETTData *rawMUSETT; 

  std::cout << tree->GetEntries() << std::endl;

  rawMUSETT = new TMUSETTData();
  tree->SetBranchAddress("MUSETT", &rawMUSETT);



  // open the output ROOT file
  TString outFileName = "~/npTreeReader/np-tree-reader/Projects/MUSETT/calibration/Time/Histograms/";
  outFileName += fname;
  outFileName += "_RawMUSETTHistos.root";
  std::cout<< outFileName<< std::endl;
  TFile *outFile = new TFile(outFileName, "recreate");

  // prepare output histograms for Must2
  TH1F* hStripXTime[NBTELESCOPE-FIRSTELESCOPE][NBSTRIPS-FIRSTSTRIP];
  TH1F* hStripYTime[NBTELESCOPE-FIRSTELESCOPE][NBSTRIPS-FIRSTSTRIP];
  for (Int_t i = FIRSTELESCOPE; i < NBTELESCOPE; i++) {
    for (Int_t j = FIRSTSTRIP; j < NBSTRIPS; j++) {
      // strips XE
      TString hnameXT     = Form("hMM%d_STRX_T%d", i, j);
      TString htitleXT    = Form("MM%d_STRX_T%d", i, j);
      hStripXTime[i][j] = new TH1F(hnameXT, htitleXT, 16384, 0, 16384);
      // strips YE
      TString hnameYT     = Form("hMM%d_STRY_T%d", i, j);
      TString htitleYT    = Form("MM%d_STRY_T%d", i, j);
      hStripYTime[i][j] = new TH1F(hnameYT, htitleYT, 16384, 0, 16384);
    }
  }

  // read the data
  Int_t nentries = tree->GetEntries();
  //   nentries = 21000;
  std::cout << "TTree contains " << nentries << " events" << std::endl;
  for (Int_t i = 0; i < nentries; i++) {
    if (i%10000 == 0) std::cout << "Entry " << i << std::endl;
    
    tree->GetEntry(i);
    Int_t multXT = 0;
    multXT = rawMUSETT->GetDSSDXTMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    for (Int_t j = 0; j < multXT; j++) {
      Int_t det    = rawMUSETT->GetDSSDXTDetectorNbr(j);
      Int_t stripX = rawMUSETT->GetDSSDXTStripNbr(j);
      Int_t Time_ = rawMUSETT->GetDSSDXTTime(j);
      if ((det < NBTELESCOPE) && (stripX < NBSTRIPS)) {
        std::cout << det << " " << stripX << " " << Time_<< std::endl;
        hStripXTime[det][stripX]->Fill(Time_);
      }
      else {
        //std::cout << "Error filling histograms: MUSETT_X_T" << std::endl;
      }
    }
    // get StripYE multiplicity
    Int_t multYT = rawMUSETT->GetDSSDYTMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    for (Int_t j = 0; j < multYT; j++) {
      Int_t det    = rawMUSETT->GetDSSDYTDetectorNbr(j);
      Int_t stripY = rawMUSETT->GetDSSDYTStripNbr(j);
      Int_t Time_ = rawMUSETT->GetDSSDYTTime(j);
      if ((det < NBTELESCOPE) && (stripY < NBSTRIPS)) {
        hStripYTime[det][stripY]->Fill(Time_);
      }
      else {
        //std::cout << "Error filling histograms: MUSETT_Y_T" << std::endl;
      }
    }
  }

  // write on disk output file and close it
  outFile->Write();
  outFile->Close();
}
