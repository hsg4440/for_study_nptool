TChain* chain;
int nentries;
TTree* otree;
double FF_AoQ;
double FF_Q;
double FF_Gamma;
int FPMW_Section;
TICPhysics* IC;

TCutG* cut[20];
TCutG* cut1;

/////////////////////////////////////////////////////////
void LoadCuts(){
  TFile* ifile = new TFile("cut_QvsAoQ.root");
  for(int i=0; i<20; i++){
    TString cut_name = Form("cut%i",i);
    cut[i] = (TCutG*)ifile->FindObjectAny(cut_name);
  }
}

/////////////////////////////////////////////////////////
void LoadRootFile(){
  chain = new TChain("PhysicsTree");
  chain->Add("../../root/analysis/run_48_init.root");

  IC = new TICPhysics();

  chain->SetBranchStatus("FF_Q13","true");
  chain->SetBranchAddress("FF_Q13",&FF_Q);

  chain->SetBranchStatus("FF_Gamma13","true");
  chain->SetBranchAddress("FF_Gamma13",&FF_Gamma);

  chain->SetBranchStatus("FF_AoQ13","true");
  chain->SetBranchAddress("FF_AoQ13",&FF_AoQ);
  
  chain->SetBranchStatus("FPMW_Section","true");
  chain->SetBranchAddress("FPMW_Section",&FPMW_Section);

  chain->SetBranchStatus("fIC","true");
  chain->SetBranchAddress("IC",&IC);
}

/////////////////////////////////////////////////////////
void InitOutputTree(){
  otree = new TTree("tree","Select Tree");

  otree->Branch("FF_Q",&FF_Q,"FF_Q/D");
  otree->Branch("FF_AoQ",&FF_AoQ,"FF_AoQ/D");
  otree->Branch("FF_Gamma",&FF_Gamma,"FF_Gamma/D");
  otree->Branch("FPMW_Section",&FPMW_Section,"FPMW_Section/I");
  otree->Branch("IC",&IC);
}

/////////////////////////////////////////////////////////
void FillSelectTree(){
  LoadRootFile();
  LoadCuts();

  nentries = chain->GetEntries();
  cout << "Number of entries: " << nentries << endl;

  TFile* ofile = new TFile("SelectTree.root","recreate");
  InitOutputTree();
 
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);

    if(i%100000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
    }
  
    if(FF_AoQ>2 && cut[FPMW_Section]->IsInside(FF_AoQ,FF_Q)){
      otree->Fill(); 
    }
 }

  otree->Write();
  delete ofile;
}
