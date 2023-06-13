#include "TPISTAData.h"

int number_of_channels = 91;
int number_of_detectors = 8;
TChain* chain;
TH1F* h[8][91];

///////////////////////////////////////////////
void FillHisto()
{
  // Input file
  chain = new TChain("RD");
  chain->Add("/home/pierre/Physics/ganil/pista/analysisenv/RootR/r0181_000r.root");

  // Output file
  TFile * ofile = new TFile("histo_file_pedestal.root","recreate");

  for(int k=0; k<number_of_detectors; k++){
    for(int i=0; i<number_of_channels; i++){
      TString histo_name = Form("h_det%i_strip%i",k+1,i+1);
      h[k][i] = new TH1F(histo_name,histo_name,200,0,200);
    }
  }


  TPISTAData* pista = new TPISTAData();
  chain->SetBranchAddress("PISTA",&pista);

  int nentries = chain->GetEntries();
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);

    if(i%10000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush; 
    }

    int mult = pista->GetPISTADEMult();
    for(int j=0; j<mult; j++){
      if(mult>0){
        int det = pista->GetPISTA_DE_DetectorNbr(j);
        int strip = pista->GetPISTA_DE_StripNbr(j);
        double val = pista->GetPISTA_DE_StripEnergy(j);
        if(det>0 && strip>0 && det<9 && strip<92)
          h[det-1][strip-1]->Fill(val);
      }
    }
  }
  ofile->Write();
  ofile->Close();

}

