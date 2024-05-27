#include "convert.h"

///////////////////////////////////////////////////
void convert(int run=29){
  m_pista = new TPISTAData();
  m_fpmw = new TFPMWData();
  m_ic = new TICData();

  // input tree //
  input_tree = new TChain("AD");
  input_tree->Add(Form("/home/morfouacep/Physics/ganil/pista/analysisenv-e850-2023/RootA/r00%i*.root",run));
  InitInputTree();


  // output tree //
  TString output_path = "/home/morfouacep/Physics/ganil/pista/analysisenv-e850-2023/np_raw/";
  TString filename = output_path+Form("run_raw_%i.root",run);
  TFile* ofile = TFile::Open(filename,"recreate");
  TTree* output_tree = new TTree("RawTree","RawTree");
  //InitOutputTree();
  output_tree->Branch("PISTA","TPISTAData",&m_pista);
  output_tree->Branch("FPMW","TFPMWData",&m_fpmw);
  output_tree->Branch("IC","TICData",&m_ic);
  output_tree->Branch("T_TMW0_FPMW0",&T_TMW0_FPMW0,"T_TMW0_FPMW0/F");
  output_tree->Branch("T_TMW0_FPMW1",&T_TMW0_FPMW1,"T_TMW0_FPMW1/F");
  output_tree->Branch("T_TMW1_FPMW0",&T_TMW1_FPMW0,"T_TMW1_FPMW0/F");
  output_tree->Branch("T_TMW1_FPMW1",&T_TMW1_FPMW1,"T_TMW1_FPMW1/F");
  output_tree->Branch("FPMWPat_0RawNr",fFPMWPat_0RawNr);
  output_tree->Branch("FPMWPat_0RawM",&fFPMWPat_0RawM);
  //output_tree->Branch("TMWPat_0TS",&TMWPat_0TS);
  output_tree->Branch("fVAMOS_TS_sec",&fVAMOS_TS_sec);
  output_tree->Branch("fPISTA_TS_sec",&fPISTA_TS_sec);
  output_tree->Branch("Exo_Mult",&Exo_Mult,"Exo_Mult/I");
  output_tree->Branch("Exo_Energy",&Exo_Energy);
  output_tree->Branch("Exo_Crystal",&Exo_Crystal);


  int nentries = input_tree->GetEntries();
  for(int i=0; i<nentries; i++){
    Clear();
    input_tree->GetEntry(i);
    if(i%100000==0) cout << "\033[34m\rProcessing tree... " <<(double)i/nentries*100 << "\% done" << flush;
   

    fVAMOS_TS_sec = TMWPat_0TS*1e-8;
    fPISTA_TS_sec = PISTA_TS*1e-8;

    // IC //
    for(int i=0; i<11; i++){
      m_ic->SetIC_Section(i+1);
      m_ic->SetIC_Charge(IC[i]);
    }

    // -TMW1- //
    for(int p=0;p<TMW1_XRawM; p++){
      int strip = TMW1_XRawNr[p];
      double charge = TMW1_XRaw[p];
      m_fpmw->SetFPMW_X(1,strip,charge);
    }
    for(int p=0;p<TMW1_YRawM; p++){
      int strip = TMW1_YRawNr[p];
      double charge = TMW1_YRaw[p];
      m_fpmw->SetFPMW_Y(1,strip,charge);
    }
    // -TMW2- //
    for(int p=0;p<TMW2_XRawM; p++){
      int strip = TMW2_XRawNr[p];
      double charge = TMW2_XRaw[p];
      m_fpmw->SetFPMW_X(2,strip,charge);
    }
    for(int p=0;p<TMW2_YRawM; p++){
      int strip = TMW2_YRawNr[p];
      double charge = TMW2_YRaw[p];
      m_fpmw->SetFPMW_Y(2,strip,charge);
    }

    // -FPMW0- //
    for(int p=0;p<FPMW0_XRawM; p++){
      int strip = FPMW0_XRawNr[p];
      double charge = FPMW0_XRaw[p];
      m_fpmw->SetFPMW_X(3,strip,charge);
    }
    for(int p=0;p<FPMW0_YRawM; p++){
      int strip = FPMW0_YRawNr[p];
      double charge = FPMW0_YRaw[p];
      m_fpmw->SetFPMW_Y(3,strip,charge);
    }
    // -FPMW1- //
    for(int p=0;p<FPMW1_XRawM; p++){
      int strip = FPMW1_XRawNr[p];
      double charge = FPMW1_XRaw[p];
      m_fpmw->SetFPMW_X(4,strip,charge);
    }
    for(int p=0;p<FPMW1_YRawM; p++){
      int strip = FPMW1_YRawNr[p];
      double charge = FPMW1_YRaw[p];
      m_fpmw->SetFPMW_Y(4,strip,charge);
    }

    if(FPMW0_XRawM>0 && FPMW1_XRawM>2 ){
      for(int p=0;p<Inner6MVM; p++){
        Exo_Mult = Inner6MVM;
        Exo_Energy.push_back(Inner6MV[p]);
        Exo_Crystal.push_back(Inner6MVN[p]);;
      }
    }
 
    if(fPISTA_TS_sec>0 || fVAMOS_TS_sec>0)
      output_tree->Fill();
  }

  output_tree->Write();
  ofile->Close();

  cout << endl;

}

///////////////////////////////////////////////////
void Clear(){
  m_fpmw->Clear();
  m_ic->Clear();
  m_pista->Clear();
  Exo_Mult=0;
  Exo_Energy.clear();
  Exo_Crystal.clear();
}


///////////////////////////////////////////////////
void InitInputTree(){
  input_tree->SetBranchAddress("PISTA",&m_pista);

  // TAC
  input_tree->SetBranchStatus("T_TMW1_FPMW0_C","true");
  input_tree->SetBranchAddress("T_TMW1_FPMW0_C",&T_TMW0_FPMW0);
  input_tree->SetBranchStatus("T_TMW1_FPMW1_C","true");
  input_tree->SetBranchAddress("T_TMW1_FPMW1_C",&T_TMW0_FPMW1);
  input_tree->SetBranchStatus("T_TMW2_FPMW0_C","true");
  input_tree->SetBranchAddress("T_TMW2_FPMW0_C",&T_TMW1_FPMW0);
  input_tree->SetBranchStatus("T_TMW2_FPMW1_C","true");
  input_tree->SetBranchAddress("T_TMW2_FPMW1_C",&T_TMW1_FPMW1);

  // Pat
  input_tree->SetBranchStatus("FPMWPat_0RawNr","true");
  input_tree->SetBranchAddress("FPMWPat_0RawNr",fFPMWPat_0RawNr);
  input_tree->SetBranchStatus("FPMWPat_0RawM","true");
  input_tree->SetBranchAddress("FPMWPat_0RawM",&fFPMWPat_0RawM);
  input_tree->SetBranchStatus("TMWPat_00TS","true");
  input_tree->SetBranchAddress("TMWPat_00TS",&TMWPat_0TS);
  input_tree->SetBranchStatus("PISTA_TS","true");
  input_tree->SetBranchAddress("PISTA_TS",&PISTA_TS);

  // IC
  input_tree->SetBranchStatus("IC","true");
  input_tree->SetBranchAddress("IC",&IC);

  // FPMW0-X
  input_tree->SetBranchStatus("FPMW0_XRawM","true");
  input_tree->SetBranchAddress("FPMW0_XRawM",&FPMW0_XRawM);
  input_tree->SetBranchStatus("FPMW0_XRaw","true");
  input_tree->SetBranchAddress("FPMW0_XRaw",FPMW0_XRaw);
  input_tree->SetBranchStatus("FPMW0_XRawNr","true");
  input_tree->SetBranchAddress("FPMW0_XRawNr",FPMW0_XRawNr);
  // FPMW0-Y
  input_tree->SetBranchStatus("FPMW0_YRawM","true");
  input_tree->SetBranchAddress("FPMW0_YRawM",&FPMW0_YRawM);
  input_tree->SetBranchStatus("FPMW0_YRaw","true");
  input_tree->SetBranchAddress("FPMW0_YRaw",FPMW0_YRaw);
  input_tree->SetBranchStatus("FPMW0_YRawNr","true");
  input_tree->SetBranchAddress("FPMW0_YRawNr",FPMW0_YRawNr);
  // FPMW1-X
  input_tree->SetBranchStatus("FPMW1_XRawM","true");
  input_tree->SetBranchAddress("FPMW1_XRawM",&FPMW1_XRawM);
  input_tree->SetBranchStatus("FPMW1_XRaw","true");
  input_tree->SetBranchAddress("FPMW1_XRaw",FPMW1_XRaw);
  input_tree->SetBranchStatus("FPMW1_XRawNr","true");
  input_tree->SetBranchAddress("FPMW1_XRawNr",FPMW1_XRawNr);
  // FPMW0-Y
  input_tree->SetBranchStatus("FPMW1_YRawM","true");
  input_tree->SetBranchAddress("FPMW1_YRawM",&FPMW1_YRawM);
  input_tree->SetBranchStatus("FPMW1_YRaw","true");
  input_tree->SetBranchAddress("FPMW1_YRaw",FPMW1_YRaw);
  input_tree->SetBranchStatus("FPMW1_YRawNr","true");
  input_tree->SetBranchAddress("FPMW1_YRawNr",FPMW1_YRawNr);
  // TMW1-X
  input_tree->SetBranchStatus("TMW1_XRawM","true");
  input_tree->SetBranchAddress("TMW1_XRawM",&TMW1_XRawM);
  input_tree->SetBranchStatus("TMW1_XRaw","true");
  input_tree->SetBranchAddress("TMW1_XRaw",TMW1_XRaw);
  input_tree->SetBranchStatus("TMW1_XRawNr","true");
  input_tree->SetBranchAddress("TMW1_XRawNr",TMW1_XRawNr);
  // TMW1-Y
  input_tree->SetBranchStatus("TMW1_YRawM","true");
  input_tree->SetBranchAddress("TMW1_YRawM",&TMW1_YRawM);
  input_tree->SetBranchStatus("TMW1_YRaw","true");
  input_tree->SetBranchAddress("TMW1_YRaw",TMW1_YRaw);
  input_tree->SetBranchStatus("TMW1_YRawNr","true");
  input_tree->SetBranchAddress("TMW1_YRawNr",TMW1_YRawNr);
  // TMW2-X
  input_tree->SetBranchStatus("TMW2_XRawM","true");
  input_tree->SetBranchAddress("TMW2_XRawM",&TMW2_XRawM);
  input_tree->SetBranchStatus("TMW2_XRaw","true");
  input_tree->SetBranchAddress("TMW2_XRaw",TMW2_XRaw);
  input_tree->SetBranchStatus("TMW2_XRawNr","true");
  input_tree->SetBranchAddress("TMW2_XRawNr",TMW2_XRawNr);
  // TMW2-Y
  input_tree->SetBranchStatus("TMW2_YRawM","true");
  input_tree->SetBranchAddress("TMW2_YRawM",&TMW2_YRawM);
  input_tree->SetBranchStatus("TMW2_YRaw","true");
  input_tree->SetBranchAddress("TMW2_YRaw",TMW2_YRaw);
  input_tree->SetBranchStatus("TMW2_YRawNr","true");
  input_tree->SetBranchAddress("TMW2_YRawNr",TMW2_YRawNr);
  // Exogam
  input_tree->SetBranchStatus("Inner6MVM","true");
  input_tree->SetBranchAddress("Inner6MVM",&Inner6MVM);
  input_tree->SetBranchStatus("Inner6MV","true");
  input_tree->SetBranchAddress("Inner6MV",Inner6MV);
  input_tree->SetBranchStatus("Inner6MVN","true");
  input_tree->SetBranchAddress("Inner6MVN",Inner6MVN);
  input_tree->SetBranchStatus("Inner6MVTS","true");
  input_tree->SetBranchAddress("Inner6MVTS",Inner6MVTS);





}

///////////////////////////////////////////////////
//void InitOutputTree(){
//  output_tree->Branch("PISTA","TPISTAData",&m_pista);
//  output_tree->Branch("FPMW","TFPMWData",&m_fpmw);
//  output_tree->Branch("IC","TICData",&m_ic);
//
//  output_tree->Branch("T_TMW0_FPMW0",&T_TMW0_FPMW0,"T_TMW0_FPMW0/F");
//  output_tree->Branch("T_TMW0_FPMW1",&T_TMW0_FPMW1,"T_TMW0_FPMW1/F");
//  output_tree->Branch("T_TMW1_FPMW0",&T_TMW1_FPMW0,"T_TMW1_FPMW0/F");
//  output_tree->Branch("T_TMW1_FPMW1",&T_TMW1_FPMW1,"T_TMW1_FPMW1/F");
//}


