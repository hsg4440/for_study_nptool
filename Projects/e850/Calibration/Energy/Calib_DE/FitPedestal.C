int NumberOfDetectors=8;
int NumberOfStrips=91;
//////////////////////////////////////////////
void FitPedestal(){
  TFile* input_file = new TFile("histo_file_pedestal.root");


  for(int i=0; i<NumberOfDetectors; i++){
    ofstream ofile;
    string ofilename = "PISTA" + to_string(i+1) + "_DE.ped";
    ofile.open(ofilename.c_str());
    for(int j=0; j<NumberOfStrips; j++){
      TString histo_name = Form("h_det%i_strip%i",i+1,j+1);
      TH1F* h = (TH1F*)gDirectory->FindObjectAny(histo_name);


      TF1* f1 = new TF1("f1","gaus");
      h->Fit(f1,"q");
      TString token = Form("PISTA_T%i_STRIP%i_DE_PEDESTAL",i+1,j+1);
      ofile << token << " " << f1->GetParameter(1) << endl;
    } 
    ofile.close();
  }

}
