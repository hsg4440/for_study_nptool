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

//#include "TMust2Data.h"
//#include "/home/e628/nptool/NPLib/include/TMust2Data.h"///
//#include "/home/e628/nptool/NPLib/include/RootInput.h"////

#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Detectors/MUST2/TMust2Data.h"
#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Detectors/Mugast/TMugastData.h"
#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/include/RootInput.h"

#define NBTELESCOPE 11
#define	NBSTRIPS	128

void ExtractMust2Histos(const char* fname = "run_0040", std::string fDet="Mugast")
{

  TString path  = "/home/mugast/Analysis_MUGAST_LISE_test/nptool/Projects/MUGAST/converted_runs/RootA/";
  TChain* tree = new TChain("AD");
  tree->Add(path + fname + "_000a.root");
  //tree->Add(path + fname + "_1.root");
  //tree->Add(path + fname + "_2.root");

  tree->SetBranchStatus("*",false);

  // connect the TMust2Data branch  
  tree->SetBranchStatus(fDet.c_str(),true);
  tree->SetBranchStatus("fMM*",true);
  TMust2Data *rawMust2;
  TMugastData *rawMugast; 

  cout << tree->GetEntries() << endl;

  if (fDet=="MUST2"){
    rawMust2 = new TMust2Data();
    tree->SetBranchAddress("MUST2", &rawMust2);
  }
  else if (fDet=="Mugast"){
    rawMugast = new TMugastData();
    tree->SetBranchAddress("Mugast", &rawMugast);
  }

  // TBranch *brMust2 = tree->Get Branch("MUST2");
  //brMust2->SetAddress(&rawMust2);


  // open the output ROOT file
  TString outFileName = "./Histograms/";
  outFileName += fname;
  outFileName += "_Raw" + fDet + "Histos.root";
  cout<< outFileName<< endl;
  TFile *outFile = new TFile(outFileName, "recreate");

  // prepare output histograms for Must2
  TH1F* hStripXTime[NBTELESCOPE][NBSTRIPS];
  TH1F* hStripYTime[NBTELESCOPE][NBSTRIPS];
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSTRIPS; j++) {
      // strips XE
      TString hnameXT     = Form("hMM%d_STRX_T%d", i+1, j+1);
      TString htitleXT    = Form("MM%d_STRX_T%d", i+1, j+1);
      hStripXTime[i][j] = new TH1F(hnameXT, htitleXT, 16384, 0, 16384);
      // strips YE
      TString hnameYT     = Form("hMM%d_STRY_T%d", i+1, j+1);
      TString htitleYT    = Form("MM%d_STRY_T%d", i+1, j+1);
      hStripYTime[i][j] = new TH1F(hnameYT, htitleYT, 16384, 0, 16384);
    }
  }

  // read the data
  Int_t nentries = tree->GetEntries();
  //   nentries = 21000;
  cout << "TTree contains " << nentries << " events" << endl;
  for (Int_t i = 0; i < nentries; i++) {
    if (i%10000 == 0) cout << "Entry " << i << endl;
    ///////////////////////
    // read Must2 branch //
    ///////////////////////
    //brMust2->GetEntry(i);
    tree->GetEntry(i);
    // get StripXE multiplicity
    Int_t multXT = 0;
    if (fDet=="MUST2"){
      multXT = rawMust2->GetMMStripXTMult();
    }else if (fDet=="Mugast"){
      multXT = rawMugast->GetDSSDXTMult();
    }
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    for (Int_t j = 0; j < multXT; j++) {
      Int_t det    = rawMust2->GetMMStripXTDetectorNbr(j);
      Int_t stripX = rawMust2->GetMMStripXTStripNbr(j);
      Int_t Time_ = rawMust2->GetMMStripXTTime(j);
      if ((det < NBTELESCOPE+1) && (stripX < NBSTRIPS+1)) {
        hStripXTime[det-1][stripX-1]->Fill(Time_);
      }
      else {
        cout << "Error filling histograms: Must2_X_T" << endl;
        //rawMust2->Dump();
      }
    }
    // get StripYE multiplicity
    Int_t multYT = rawMust2->GetMMStripYTMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    for (Int_t j = 0; j < multYT; j++) {
      Int_t det    = rawMust2->GetMMStripYTDetectorNbr(j);
      Int_t stripY = rawMust2->GetMMStripYTStripNbr(j);
      Int_t Time_ = rawMust2->GetMMStripYTTime(j);
      if ((det < NBTELESCOPE+1) && (stripY < NBSTRIPS+1)) {
        hStripYTime[det-1][stripY-1]->Fill(Time_);
      }
      else {
        cout << "Error filling histograms: Must2_Y_T" << endl;
        //rawMust2->Dump();
      }
    }
  }

  // write on disk output file and close it
  outFile->Write();
  outFile->Close();
}
