TH1F* h1[20];
TCanvas* c1[20];
void Fit(TH1F* histo,int section);

ofstream ofile;

//////////////////////////////////////////////
void LoadHisto(){
  TFile* ifile = new TFile("histo/charge.root","read");
  for(int i=0; i<20; i++){
    TString histo_name = Form("hQ_sec%i",i);
    h1[i] = (TH1F*) ifile->FindObjectAny(histo_name);
  }
}

//////////////////////////////////////////////
void FitCharge(){
  LoadHisto();

  string filename = "charge.cal";
  ofile.open(filename.c_str());
  for(int i=0; i<20; i++){
    Fit(h1[i],i);
  }

  ofile.close();
}

//////////////////////////////////////////////
void Fit(TH1F* histo, int section){
  TString canvas_name = Form("c1_sec%i",section);
  c1[section] = new TCanvas(canvas_name,canvas_name,1200,600);
  c1[section]->Divide(2,1);
  c1[section]->cd(1);
  histo->Draw();

  TGraphErrors* gerr = new TGraphErrors();

  int NumberOfQ = 9;
  TF1* g[NumberOfQ];
  Double_t para[3*NumberOfQ];
  double mean = 1435;
  for(int i=0; i<NumberOfQ; i++){
    g[i] = new TF1(Form("g%i",i),"gaus",mean-20,mean+20);
    g[i]->SetParLimits(0,0.1,5000);
    g[i]->SetParLimits(1,mean-20,mean+20);
    g[i]->SetParLimits(2,5,11);

    histo->Fit(g[i],"qR");
    g[i]->Draw("lsame");
    g[i]->GetParameters(&para[3*i]);

    mean+=41;
  }

  TString total_func = "gaus(0)";
  for(int i=0; i<NumberOfQ; i++){
    TString gaus_func = Form("gaus(%i)",3*i);
    total_func += "+" + gaus_func;
  }
  TF1* total = new TF1("total",total_func,1400,1800);
  total->SetParameters(para);
  for(int i=0; i<NumberOfQ; i++){
    total->SetParLimits(3*i+2,5,11);
  }
  total->SetLineColor(4);
  histo->Fit(total,"Rq");

  int Q[20] = {34,34,34,34,34,34,34,35,34,36,36,36,36,36,36,36,35,35,35,35};
  int Qi = 33;
  for(int i=0; i<NumberOfQ; i++){
    double mean = total->GetParameter(3*i+1);
    double mean_err = total->GetParError(3*i+1);

    gerr->SetPoint(i,mean,Q[section]-1);
    gerr->SetPointError(i,mean_err,0);
    Q[section]++;
    Qi++;
  }

  gerr->SetMarkerStyle(8);
  gerr->SetMarkerSize(1.5);
  c1[section]->cd(2);
  gerr->Draw("ap");

  TF1* my_func = new TF1("f1","pol1",1000,2000);
  gerr->Fit("pol1");


  double p0 = gerr->GetFunction("pol1")->GetParameter(0);
  double p1 = gerr->GetFunction("pol1")->GetParameter(1);

  ofile << "SECTION" << section << " " << p0 << " " << p1 << endl;


}
