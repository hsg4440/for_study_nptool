TChain* chain;
TCutG* cut10Be[8];

/////////////////////////////////////////////////////////////////
void LoadCut(string part){
  string path = "cut/"+ part + "/cut10Be.root";
  TFile* fcut = new TFile(path.c_str());
  for(int i=0; i<8; i++){
   cut10Be[i] = (TCutG*)fcut->Get(Form("cut10Be_det%i",i+1));
  }
}

/////////////////////////////////////////////////////////////////
void FillHistoMass(string part="part2"){

  LoadCut(part);

  chain = new TChain("PhysicsTree");
  if(part=="part1"){
    chain->Add("../../root/analysis/run_28.root");
    chain->Add("../../root/analysis/run_29.root");
    chain->Add("../../root/analysis/run_30.root");
    chain->Add("../../root/analysis/run_31.root");
    chain->Add("../../root/analysis/run_32.root");
    chain->Add("../../root/analysis/run_36.root");
    chain->Add("../../root/analysis/run_37.root");
    chain->Add("../../root/analysis/run_39.root");
    chain->Add("../../root/analysis/run_40.root");
    chain->Add("../../root/analysis/run_41.root");
    chain->Add("../../root/analysis/run_42.root");
    chain->Add("../../root/analysis/run_43.root");
    chain->Add("../../root/analysis/run_44.root");
    chain->Add("../../root/analysis/run_45.root");
    chain->Add("../../root/analysis/run_46.root");
    chain->Add("../../root/analysis/run_48.root");
    chain->Add("../../root/analysis/run_49.root");
    chain->Add("../../root/analysis/run_50.root");
    chain->Add("../../root/analysis/run_51.root");
    chain->Add("../../root/analysis/run_52.root");
    chain->Add("../../root/analysis/run_53.root");
  }
  else if(part=="part2"){
    chain->Add("../../root/analysis/run_55.root");
    chain->Add("../../root/analysis/run_56.root");
    chain->Add("../../root/analysis/run_57.root");
    chain->Add("../../root/analysis/run_58.root");
    chain->Add("../../root/analysis/run_59.root");
    chain->Add("../../root/analysis/run_60.root");
    chain->Add("../../root/analysis/run_61.root");
    chain->Add("../../root/analysis/run_62.root");
    chain->Add("../../root/analysis/run_63.root");
    chain->Add("../../root/analysis/run_64.root");
    chain->Add("../../root/analysis/run_65.root");
    chain->Add("../../root/analysis/run_66.root");
    chain->Add("../../root/analysis/run_67.root");
    chain->Add("../../root/analysis/run_68.root");
    chain->Add("../../root/analysis/run_69.root");
  }
  double DeltaEcorr;
  double Elab;
  double Mass;
  double Ex240Pu;
  double fTS_TMW;
  int Telescope;

  chain->SetBranchStatus("DeltaEcorr","true");
  chain->SetBranchAddress("DeltaEcorr",&DeltaEcorr);
  chain->SetBranchStatus("Elab","true");
  chain->SetBranchAddress("Elab",&Elab);
  chain->SetBranchStatus("Mass","true");
  chain->SetBranchAddress("Mass",&Mass);
  chain->SetBranchStatus("Ex240Pu","true");
  chain->SetBranchAddress("Ex240Pu",&Ex240Pu);
  chain->SetBranchStatus("fTS_TMW","true");
  chain->SetBranchAddress("fTS_TMW",&fTS_TMW);
  chain->SetBranchStatus("Telescope","true");
  chain->SetBranchAddress("Telescope",&Telescope);



  TH1F *hmass[11];
  TH2F *hPID[8];
  for(int i=0; i<11; i++){
    TString histo_name = Form("hmass_%iMeV",i+6);
    hmass[i] = new TH1F(histo_name, histo_name,1000,80,160); 
  }
  for(int i=0; i<8; i++){
    TString histo_name = Form("PID_det%i",i+1);
    hPID[i] = new TH2F(histo_name,histo_name,1000,0,200,1000,0,60);
  }
  //TH2F* hPID = new TH2F("PID","PID",1000,0,200,1000,0,60);
  TH1F* hEx = new TH1F("Ex","Ex",500,0,20);
  TH1F* hMassAll = new TH1F("hmass_all","hmass_all",1000,80,160);


  int nentries = chain->GetEntries();
  cout << "**** Number of entries to be treated: " << nentries << " ****" << endl; 
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    if(i%1000000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
    }

    if(fTS_TMW>0){

      if(Telescope>0) hPID[Telescope-1]->Fill(Elab,DeltaEcorr);

      if(cut10Be[Telescope-1]->IsInside(Elab,DeltaEcorr)){
        hEx->Fill(Ex240Pu);
        hMassAll->Fill(Mass);
        if(Ex240Pu>5.5 && Ex240Pu<6.5) hmass[0]->Fill(Mass);
        else if(Ex240Pu>6.5 && Ex240Pu<7.5) hmass[1]->Fill(Mass);
        else if(Ex240Pu>7.5 && Ex240Pu<8.5) hmass[2]->Fill(Mass);
        else if(Ex240Pu>8.5 && Ex240Pu<9.5) hmass[3]->Fill(Mass);
        else if(Ex240Pu>9.5 && Ex240Pu<10.5) hmass[4]->Fill(Mass);
        else if(Ex240Pu>10.5 && Ex240Pu<11.5) hmass[5]->Fill(Mass);
        else if(Ex240Pu>11.5 && Ex240Pu<12.5) hmass[6]->Fill(Mass);
        else if(Ex240Pu>12.5 && Ex240Pu<13.5) hmass[7]->Fill(Mass);
        else if(Ex240Pu>13.5 && Ex240Pu<14.5) hmass[8]->Fill(Mass);
        else if(Ex240Pu>14.5 && Ex240Pu<15.5) hmass[9]->Fill(Mass);
        else if(Ex240Pu>15.5 && Ex240Pu<16.5) hmass[10]->Fill(Mass);
      }
    }
  }

  TString filename = Form("histo_mass_%s.root",part.c_str());
  TFile *ofile = new TFile(filename,"recreate");
  for(int i=0; i<11; i++){
    hmass[i]->Write();
  }
  hEx->Write();
  hMassAll->Write();
  for(int i=0; i<8; i++)
    hPID[i]->Write();
  ofile->Close();



}
