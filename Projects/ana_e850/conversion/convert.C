#include "convert.h"

///////////////////////////////////////////////////
void convert(int run=29){
  m_pista = new TPISTAData();
  m_fpmw = new TFPMWData();
  m_ic = new TICData();

  // input tree //
  input_tree = new TChain("AD");
  input_tree->Add(Form("/run/media/morfouacep/e850/ganil/run/e850/RootA/r00%i*",run));
  InitInputTree();


  // output tree //
  TString output_path = "/run/media/morfouacep/e850/ganil/run/e850/np_raw/";
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
  output_tree->Branch("TMWPat_0TS",&TMWPat_0TS);
  output_tree->Branch("Exo_Mult",&Exo_Mult,"Exo_Mult/I");
  output_tree->Branch("Exo_Energy",&Exo_Energy);
  output_tree->Branch("Exo_Crystal",&Exo_Crystal);


  int nentries = input_tree->GetEntries();
  for(int i=0; i<nentries; i++){
    Clear();
    input_tree->GetEntry(i);
    if(i%10000==0) cout << "\033[34m\rProcessing tree... " <<(double)i/nentries*100 << "\% done" << flush;
   


    // IC //
    for(int i=0; i<11; i++){
      m_ic->SetIC_Section(i+1);
      m_ic->SetIC_Charge(IC[i]);
    }

    // -TMW1- //
    for(int p=0;p<TMW1_XVM; p++){
      int strip = TMW1_XVN[p];
      double charge = TMW1_XV[p];
      m_fpmw->SetFPMW_X(0,strip,charge);
    }
    for(int p=0;p<TMW1_YVM; p++){
      int strip = TMW1_YVN[p];
      double charge = TMW1_YV[p];
      m_fpmw->SetFPMW_Y(0,strip,charge);
    }
    // -TMW2- //
    for(int p=0;p<TMW2_XVM; p++){
      int strip = TMW2_XVN[p];
      double charge = TMW2_XV[p];
      m_fpmw->SetFPMW_X(1,strip,charge);
    }
    for(int p=0;p<TMW2_YVM; p++){
      int strip = TMW2_YVN[p];
      double charge = TMW2_YV[p];
      m_fpmw->SetFPMW_Y(1,strip,charge);
    }

    // -FPMW0- //
    for(int p=0;p<FPMW0_XVM; p++){
      int strip = FPMW0_XVN[p];
      double charge = FPMW0_XV[p];
      m_fpmw->SetFPMW_X(2,strip,charge);
    }
    for(int p=0;p<FPMW0_YVM; p++){
      int strip = FPMW0_YVN[p];
      double charge = FPMW0_YV[p];
      m_fpmw->SetFPMW_Y(2,strip,charge);
    }
    // -FPMW1- //
    for(int p=0;p<FPMW1_XVM; p++){
      int strip = FPMW1_XVN[p];
      double charge = FPMW1_XV[p];
      m_fpmw->SetFPMW_X(3,strip,charge);
    }
    for(int p=0;p<FPMW1_YVM; p++){
      int strip = FPMW1_YVN[p];
      double charge = FPMW1_YV[p];
      m_fpmw->SetFPMW_Y(3,strip,charge);
    }

    for(int p=0;p<Inner6MVM; p++){
      Exo_Mult = Inner6MVM;
      Exo_Energy.push_back(Inner6MV[p]);
      Exo_Crystal.push_back(Inner6MVN[p]);;
    }
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



  // IC
  input_tree->SetBranchStatus("IC","true");
  input_tree->SetBranchAddress("IC",&IC);

  // FPMW0-X
  input_tree->SetBranchStatus("FPMW0_XVM","true");
  input_tree->SetBranchAddress("FPMW0_XVM",&FPMW0_XVM);
  input_tree->SetBranchStatus("FPMW0_XV","true");
  input_tree->SetBranchAddress("FPMW0_XV",FPMW0_XV);
  input_tree->SetBranchStatus("FPMW0_XVN","true");
  input_tree->SetBranchAddress("FPMW0_XVN",FPMW0_XVN);
  // FPMW0-Y
  input_tree->SetBranchStatus("FPMW0_YVM","true");
  input_tree->SetBranchAddress("FPMW0_YVM",&FPMW0_YVM);
  input_tree->SetBranchStatus("FPMW0_YV","true");
  input_tree->SetBranchAddress("FPMW0_YV",FPMW0_YV);
  input_tree->SetBranchStatus("FPMW0_YVN","true");
  input_tree->SetBranchAddress("FPMW0_YVN",FPMW0_YVN);
  // FPMW1-X
  input_tree->SetBranchStatus("FPMW1_XVM","true");
  input_tree->SetBranchAddress("FPMW1_XVM",&FPMW1_XVM);
  input_tree->SetBranchStatus("FPMW1_XV","true");
  input_tree->SetBranchAddress("FPMW1_XV",FPMW1_XV);
  input_tree->SetBranchStatus("FPMW1_XVN","true");
  input_tree->SetBranchAddress("FPMW1_XVN",FPMW1_XVN);
  // FPMW0-Y
  input_tree->SetBranchStatus("FPMW1_YVM","true");
  input_tree->SetBranchAddress("FPMW1_YVM",&FPMW1_YVM);
  input_tree->SetBranchStatus("FPMW1_YV","true");
  input_tree->SetBranchAddress("FPMW1_YV",FPMW1_YV);
  input_tree->SetBranchStatus("FPMW1_YVN","true");
  input_tree->SetBranchAddress("FPMW1_YVN",FPMW1_YVN);
  // TMW1-X
  input_tree->SetBranchStatus("TMW1_XVM","true");
  input_tree->SetBranchAddress("TMW1_XVM",&TMW1_XVM);
  input_tree->SetBranchStatus("TMW1_XV","true");
  input_tree->SetBranchAddress("TMW1_XV",TMW1_XV);
  input_tree->SetBranchStatus("TMW1_XVN","true");
  input_tree->SetBranchAddress("TMW1_XVN",TMW1_XVN);
  // TMW1-Y
  input_tree->SetBranchStatus("TMW1_YVM","true");
  input_tree->SetBranchAddress("TMW1_YVM",&TMW1_YVM);
  input_tree->SetBranchStatus("TMW1_YV","true");
  input_tree->SetBranchAddress("TMW1_YV",TMW1_YV);
  input_tree->SetBranchStatus("TMW1_YVN","true");
  input_tree->SetBranchAddress("TMW1_YVN",TMW1_YVN);
  // TMW2-X
  input_tree->SetBranchStatus("TMW2_XVM","true");
  input_tree->SetBranchAddress("TMW2_XVM",&TMW2_XVM);
  input_tree->SetBranchStatus("TMW2_XV","true");
  input_tree->SetBranchAddress("TMW2_XV",TMW2_XV);
  input_tree->SetBranchStatus("TMW2_XVN","true");
  input_tree->SetBranchAddress("TMW2_XVN",TMW2_XVN);
  // TMW2-Y
  input_tree->SetBranchStatus("TMW2_YVM","true");
  input_tree->SetBranchAddress("TMW2_YVM",&TMW2_YVM);
  input_tree->SetBranchStatus("TMW2_YV","true");
  input_tree->SetBranchAddress("TMW2_YV",TMW2_YV);
  input_tree->SetBranchStatus("TMW2_YVN","true");
  input_tree->SetBranchAddress("TMW2_YVN",TMW2_YVN);
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


