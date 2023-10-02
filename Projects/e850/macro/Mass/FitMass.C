void Fit(TH1F* hmass, string Energy);
TFile* ofile;

//////////////////////////////////////////////////////
void FitMass(){

  TFile* ifile = new TFile("histo_mass_all.root","read");

  ofile = new TFile("histo_mass_fitted.root","recreate");

  for(int i=0; i<10; i++){
    int Ex = i+6;
    string Energy = to_string(Ex) + "MeV";
    TString histo_name = Form("hmass_%iMeV",i+6);
    TH1F* h1 = (TH1F*)ifile->FindObjectAny(histo_name);
    Fit(h1,Energy);
  }

  ofile->Close();
}

//////////////////////////////////////////////////////
void Fit(TH1F* hmass,string Energy){
  hmass->Rebin(2);
  hmass->Draw();

  TGraph* gyield = new TGraph();
  int Amin = 82;
  int Amax = 153;

  int A = Amin;

  int NumberOfA = Amax - Amin;
  Double_t para[3*NumberOfA];

  TF1* g[NumberOfA];
  double Integral[NumberOfA];
  double Integral_err[NumberOfA];
  double total_integral = 0;
  double Yield[NumberOfA];

  for(int i=0; i<NumberOfA; i++){
    g[i] = new TF1(Form("g%i",i),"gaus",A-0.4,A+0.4);   
    g[i]->SetParLimits(0,0.1,20000);
    g[i]->SetParLimits(1,A-0.1,A+0.1);
    g[i]->SetParLimits(2,0.1,0.45);
  
    hmass->Fit(g[i],"qR");

    g[i]->Draw("lsame");
    g[i]->GetParameters(&para[3*i]);

    A++;

  }


  TString total_func = "gaus(0)";
  for(int i=1; i<NumberOfA; i++){
    TString gaus_func = Form("gaus(%i)",3*i);
    total_func += "+" + gaus_func;
  }

  TF1* total = new TF1("total",total_func,Amin-0.5, Amin+NumberOfA-0.5);

  A = Amin;
  for(int i=0; i<NumberOfA; i++){
    total->SetParLimits(3*i+1,A-0.05,A+0.05);
    total->SetParLimits(3*i+2,0.25,0.35);
    A++;
  }
  total->SetParameters(para);
  total->SetLineColor(4);

  hmass->Fit(total,"Rq");

  A = Amin;
  for(int i=0; i<NumberOfA; i++){
      
    double Amplitude = total->GetParameter(3*i);
    double mean = total->GetParameter(3*i+1);
    double sigma = total->GetParameter(3*i+2);
    double Amplitude_err = total->GetParError(3*i);
    double sigma_err = total->GetParError(3*i+2);
    
    /*cout << "A=  " << A << endl;
    cout << "Amplitude= " << Amplitude << endl;
    cout << "Mean= " << mean << endl;
    cout << "Sigma= " << sigma << endl;*/

    Integral[i] = Amplitude*sigma*sqrt(2*TMath::Pi());
    Integral_err[i] = sqrt(2*TMath::Pi()*(pow(sigma*Amplitude_err,2) + pow(Amplitude*sigma_err,2)));
    //Integral_err[i] = sqrt(Integral[i]);
    

    total_integral += Integral[i];
    A++;
  }

  hmass->Write();
  total->Write();

  string output_filename = "Mass_Yield_" + Energy + ".dat";
  ofstream ofile;
  ofile.open(output_filename.c_str());

  // Yield calculation //
  A = Amin;
  for(int i=0; i<NumberOfA; i++){
    Yield[i] = Integral[i]/total_integral*200;
    gyield->SetPoint(i,A,Yield[i]);
    double yield_err = Integral_err[i]/total_integral*200;
    ofile << A << " " << Yield[i] << " " << yield_err << endl;
    if(A==99) cout << "A=99 -> Y= " << Yield[i] << endl;
    if(A==147) cout << "A=147 -> Y= " << Yield[i] << endl;
    A++;
  }
  ofile.close();

  gyield->SetMarkerStyle(8);
  gyield->SetMarkerSize(1);

}
