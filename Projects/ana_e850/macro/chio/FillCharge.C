TChain* chain;
int nentries;
TH1F* hQ[20];

double FF_Q;
int FPMW_Section;


/////////////////////////////////////////////////////////
void LoadRootFile(){
  chain = new TChain("PhysicsTree");
  //chain->Add("../../root/analysis/run_48.root");
  chain->Add("../../root/analysis/test.root");
}

/////////////////////////////////////////////////////////
void FillCharge(){
  LoadRootFile();
  nentries = chain->GetEntries();
  cout << "Number of entries: " << nentries << endl;

  TFile* ofile = new TFile("histo/charge.root","recreate");

  for(int i=0; i<20; i++){
    TString histo_name = Form("hQ_sec%d",i);
    hQ[i] = new TH1F(histo_name,histo_name,1000,1000,2000);
  }

  chain->SetBranchStatus("FF_Q","true");
  chain->SetBranchAddress("FF_Q",&FF_Q);

  chain->SetBranchStatus("FPMW_Section","true");
  chain->SetBranchAddress("FPMW_Section",&FPMW_Section);

  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);

    if(i%100000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
    }

    hQ[FPMW_Section]->Fill(FF_Q);
  }

  ofile->Write();
  ofile->Close();
}
