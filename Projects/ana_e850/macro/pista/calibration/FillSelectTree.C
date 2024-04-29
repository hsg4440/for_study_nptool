TChain* chain;
int nentries;
TTree* otree;
double DeltaE;
double Eres;
double Xcalc;
double Ycalc;
double Zcalc;
int Telescope;
double Brho;
int m_2alpha;
double XTarget;
double YTarget;

TCutG* cut_12C[8];

/////////////////////////////////////////////////////////
void LoadCuts(){
  TFile* ifile = new TFile("cut/cut_12C.root");
  for(int i=0; i<8; i++){
    TString cutname = Form("cut_det%d",i+1);
    cut_12C[i] = (TCutG*)ifile->FindObjectAny(cutname);
  }
}

/////////////////////////////////////////////////////////
void LoadRootFile(){
  chain = new TChain("PhysicsTree");
  chain->Add("../../../root/analysis/run_55_Eraw.root");

  chain->SetBranchStatus("FF_Brho","true");
  chain->SetBranchAddress("FF_Brho",&Brho);

  chain->SetBranchStatus("DeltaE","true");
  chain->SetBranchAddress("DeltaE",&DeltaE);

  chain->SetBranchStatus("Eres","true");
  chain->SetBranchAddress("Eres",&Eres);

  chain->SetBranchStatus("Xcalc","true");
  chain->SetBranchAddress("Xcalc",&Xcalc);
  chain->SetBranchStatus("Ycalc","true");
  chain->SetBranchAddress("Ycalc",&Ycalc);
  chain->SetBranchStatus("Zcalc","true");
  chain->SetBranchAddress("Zcalc",&Zcalc);
 
  chain->SetBranchStatus("XTarget","true");
  chain->SetBranchAddress("XTarget",&XTarget);
  chain->SetBranchStatus("YTarget","true");
  chain->SetBranchAddress("YTarget",&YTarget);
 
  chain->SetBranchStatus("Telescope","true");
  chain->SetBranchAddress("Telescope",&Telescope);
 
  chain->SetBranchStatus("m_2alpha","true");
  chain->SetBranchAddress("m_2alpha",&m_2alpha);
 
}

/////////////////////////////////////////////////////////
void InitOutputTree(){
  otree = new TTree("tree","Select Tree");

  otree->Branch("DeltaE",&DeltaE,"DeltaE/D");
  otree->Branch("Eres",&Eres,"Eres/D");
  otree->Branch("Xcalc",&Xcalc,"Xcalc/D");
  otree->Branch("Ycalc",&Ycalc,"Ycalc/D");
  otree->Branch("Zcalc",&Zcalc,"Zcalc/D");
  otree->Branch("XTarget",&XTarget,"XTarget/D");
  otree->Branch("YTarget",&YTarget,"YTarget/D");
  otree->Branch("Telescope",&Telescope,"Telescope/I");
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

    if(Telescope>0 && Brho<0 && Xcalc != -1000 ){
      if(cut_12C[Telescope-1]->IsInside(Eres,DeltaE)){
        otree->Fill(); 
      }
    }
  }

  otree->Write();
  delete ofile;
}
