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


#include "../../../../NPLib/Detectors/MUST2/TMust2Data.h"
#include "../../../../NPLib/Detectors/Mugast/TMugastData.h"
#include "../../../../NPLib/include/RootInput.h"

#define NBTELESCOPE 4	
#define	NBSTRIPS	128
#define NBSILI     0


void ExtractRawHistos(const char* fname = "run_XXXX_X", std::string fDet="Mugast")
{
  RootInput* Input = RootInput::getInstance("RunToTreat.txt");
  TChain* Chain = Input->GetChain();

  // connect the TfDet.c_str()Data branch  
  Chain->SetBranchStatus(fDet.c_str(),true);
  Chain->SetBranchStatus("fMM*",true);

  TMust2Data *rawMust2;
  TMugastData *rawMugast;
  if (fDet=="MUST2"){
    rawMust2 = new TMust2Data();
    Chain->SetBranchAddress(fDet.c_str(),&rawMust2);
  }
  else if (fDet=="Mugast"){
    rawMugast = new TMugastData();
    Chain->SetBranchAddress(fDet.c_str(),&rawMugast);
  }

  // open the output ROOT file
  TString outFileName = "./Histograms/";
  outFileName += fname;
  outFileName += "_Raw" + fDet + "Histos.root";
  cout<< outFileName<< endl;
  TFile *outFile = new TFile(outFileName, "recreate");

  // prepare output histograms for Must2
  TH1F* hStripXEnergy[NBTELESCOPE][NBSTRIPS];
  TH1F* hStripYEnergy[NBTELESCOPE][NBSTRIPS];
  //from i=4 because I wanted only the T5,T8
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSTRIPS; j++) {
      // strips XE
      TString hnameXE     = Form("hMM%d_STRX_E%d", i+1, j+1);
      TString htitleXE    = Form("MM%d_STRX_E%d", i+1, j+1);
      hStripXEnergy[i][j] = new TH1F(hnameXE, htitleXE, 16384, 0, 16384);
      // strips YE
      TString hnameYE     = Form("hMM%d_STRY_E%d", i+1, j+1);
      TString htitleYE    = Form("MM%d_STRY_E%d", i+1, j+1);
      hStripYEnergy[i][j] = new TH1F(hnameYE, htitleYE, 16384, 0, 16384);
    }
  }
  TH1F* hSiLiEnergy[NBTELESCOPE][NBSILI];
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSILI; j++) {
      TString hnameSiLiE     = Form("hMM%d_SILI_E%d", i+1, j+1);
      TString htitleSiLiE    = Form("MM%d_SILI_E%d", i+1, j+1);
      hSiLiEnergy[i][j] = new TH1F(hnameSiLiE, htitleSiLiE, 16384, 0, 16384);
    }
  }

  // read the data
  Int_t nentries = Chain->GetEntries();
  cout << "TChain contains " << nentries << " events" << endl;
  for (Int_t i = 0; i < nentries; i++) {
    if (i%1000 == 0) cout << "\rEntry " << i << "\t" << i/(double)nentries * 100 << " %" << flush;
    ///////////////////////
    // read Must2 branch //
    ///////////////////////
    //brMust2->GetEntry(i);
    Chain->GetEntry(i);
    // get StripXE multiplicity
    Int_t multXE = 0; 
    if (fDet=="MUST2"){
      multXE = rawMust2->GetMMStripXEMult();
    }else if (fDet=="Mugast"){
      multXE = rawMugast->GetDSSDXEMult();
    }
    //cout<< "multXE" << multXE<< endl;
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    //
    Int_t det = 0;   
    Int_t stripX = 0;
    Int_t energy = 0;

    for (Int_t j = 0; j < multXE; j++) {
      if (fDet=="MUST2"){
        det    = rawMust2->GetMMStripXEDetectorNbr(j);
        stripX = rawMust2->GetMMStripXEStripNbr(j);
        energy = rawMust2->GetMMStripXEEnergy(j);
      } else if (fDet=="Mugast"){
        det    = rawMugast->GetDSSDXEDetectorNbr(j);
        stripX = rawMugast->GetDSSDXEStripNbr(j);
        energy = rawMugast->GetDSSDXEEnergy(j);
      }
      if ((det < NBTELESCOPE+1) && (stripX < NBSTRIPS+1)) {
        hStripXEnergy[det-1][stripX-1]->Fill(energy);
      } else {
        cout << "Error filling histograms: X_E" << endl;
        if (fDet=="MUST2"){
          rawMust2->Dump();
        } else if (fDet=="Mugast"){
          rawMugast->Dump();
        }
      }
    }
    // get StripYE multiplicity
    Int_t multYE=0;
    if (fDet=="MUST2"){
      multYE = rawMust2->GetMMStripYEMult();
    } else if (fDet=="Mugast"){
      multYE = rawMugast->GetDSSDYEMult();
    }
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms

    det = 0;   
    Int_t stripY = 0;
    energy = 0;
    for (Int_t j = 0; j < multYE; j++) {
      if (fDet=="MUST2"){
        det    = rawMust2->GetMMStripYEDetectorNbr(j);
        stripY = rawMust2->GetMMStripYEStripNbr(j);
        energy = rawMust2->GetMMStripYEEnergy(j);
      }else if (fDet=="Mugast"){
        det    = rawMugast->GetDSSDYEDetectorNbr(j);
        stripY = rawMugast->GetDSSDYEStripNbr(j);
        energy = rawMugast->GetDSSDYEEnergy(j);
      }

      if ((det < NBTELESCOPE+1) && (stripY < NBSTRIPS+1)) {
        hStripYEnergy[det-1][stripY-1]->Fill(energy);
      }
      else {
        cout << "Error filling histograms: Y_E" << endl;
        if (fDet=="MUST2"){
          rawMust2->Dump();
        }else if (fDet=="Mugast"){
          rawMugast->Dump();
        }
      }
    }
    // get SiLiE multiplicity

    //Int_t multSiLiE = rawMust2->GetMMSiLiEMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    /* for (Int_t j = 0; j < multSiLiE; j++) {
       Int_t det    = rawMust2->GetMMSiLiEDetectorNbr(j);
       Int_t pad = rawMust2->GetMMSiLiEPadNbr(j);
       Int_t energy = rawMust2->GetMMSiLiEEnergy(j);
       if ((det < NBTELESCOPE+1) && (pad < NBSTRIPS+1)) {
       hSiLiEnergy[det-1][pad-1]->Fill(energy);
       }
       else {
       cout << "Error filling histograms: Must2_SiLi_E" << endl;
       rawMust2->Dump();
       }
       }*/

  }
  cout << '\n';

  // write on disk output file and close it
  outFile->Write();
  outFile->Close();
}
