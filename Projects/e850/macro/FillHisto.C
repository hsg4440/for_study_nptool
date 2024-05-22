int number_of_channels = 57;
int number_of_detectors = 8;
TChain* chain;
TH2F* h[8][57];
TGraph* g[8][2];

///////////////////////////////////////////////
void FillHisto()
{
  // Input file
  chain = new TChain("PhysicsTree");
  chain->Add("/home/e850/workdir/Pierre/np_e850/Outputs/Analysis/run_043*");

  // Output file
  TFile * ofile = new TFile("histo_file_run43.root","recreate");

  for(int k=0; k<number_of_detectors; k++){
    for(int i=0; i<number_of_channels; i++){
      TString histo_name = Form("h_det%i_strip%i",k+1,i+1);
      h[k][i] = new TH2F(histo_name,histo_name,1000,0,4000,500,0,200);
    }
  }


  TPISTAPhysics* pista = new TPISTAPhysics();
  chain->SetBranchAddress("PISTA",&pista);

  int nentries = chain->GetEntries();
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);

    if(i%100000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush; 
    }

    int mult = pista->EventMultiplicity;
    if(mult==1){
      int det = pista->DetectorNumber[0];
      int strip = pista->E_StripNbr[0];
      double E_strip = pista->E[0];
      double E_back = pista->back_E[0];
      if(det>0 && strip>0 && det<9 && strip<58 && E_strip>300)
        h[det-1][strip-1]->Fill(E_strip,E_back);
    }
  }

  TF1* f1 = new TF1("f1","[0] + [1]*x",0,4000);
  double p0, p1;
  for(int i=0; i<number_of_detectors; i++){
    ofstream ofile1;
    string ofilename1 = "PISTA" + to_string(i+1) + "_E_run43.cal";
    ofile1.open(ofilename1.c_str());
    
    double poles0[number_of_channels], poles1[number_of_channels], strips[number_of_channels];

    for(int j=0; j<number_of_channels; j++){
      TString token = Form("PISTA_T%i_STRIP%i_E_ENERGY",i+1,j+1);

      int N = h[i][j]->GetEntries();
      if(N>0){
        h[i][j]->SetMinimum(2);
        h[i][j]->Fit("f1","q");
        p0 = f1->GetParameter(0);
        p1 = f1->GetParameter(1);
      }
      else{
        p0 = 0;
        p1 = 1;
      }
      ofile1 << token << " " << p0 << " " << p1 << endl;
      poles0[j] = p0;
      poles1[j] = p1;
      strips[j] = j;
	
    }
    ofile1.close();
    // graphs for pole0
    auto c0 = new TCanvas(Form("p0_d%d",i+1), Form("p0_d%d",i+1));
    g[i][0] = new TGraph(number_of_channels, strips, poles0);
    g[i][0]->SetName(Form("p0_detector%d",i+1));
    g[i][0]->SetTitle(Form("pole0 par strip, detecteur %d", i+1));
    g[i][0]->Draw();
    g[i][0]->Write();
    c0->Close();
    // graphs for pole1
    auto c1 = new TCanvas(Form("p1_d%d",i+1), Form("p1_d%d",i+1));
    g[i][1] = new TGraph(number_of_channels, strips, poles1);
    g[i][1]->SetName(Form("p1_detector%d",i+1));
    g[i][1]->SetTitle(Form("pole1 par strip, detecteur %d", i+1));
    g[i][1]->Draw();
    g[i][1]->Write();
    c1->Close();
  }


  ofile->Write();
  ofile->Close();


}

