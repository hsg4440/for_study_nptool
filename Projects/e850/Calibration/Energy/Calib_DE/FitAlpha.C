#include "Calibrator.h"

int NumberOfDetectors = 8;
int NumberOfStrips = 91;

double pedestal[91];

///////////////////////////////////////////////
void FitAlpha()
{
  TFile* input_file = new TFile("histo_file_alpha.root");

  DefineSource();
 
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TCanvas* c2 = new TCanvas("c2","c2",600,600);


  for(int i=0; i<NumberOfDetectors; i++){
    /*ifstream ifile;
    string ifilename = "PISTA" + to_string(i+1) + "_DE.ped";
    ifile.open(ifilename.c_str());
    string token;
    for(int s=0; s<NumberOfStrips; s++){
      ifile >> token >> pedestal[s];

      cout << token << " " << pedestal[s] << endl;
    }
    ifile.close();*/

    ofstream ofile;
    string ofilename = "PISTA" + to_string(i+1) + "_DE.cal";
    ofile.open(ofilename.c_str());
    for(int j=0; j<NumberOfStrips; j++){
      c1->cd();
      TString histo_name = Form("h_det%i_strip%i",i+1,j+1);
      TH1F* h = (TH1F*)input_file->FindObjectAny(histo_name);
      h->Draw();
      
      TString token = Form("PISTA_T%i_STRIP%i_DE_ENERGY",i+1,j+1);
      if(Finder(h,mean,sigma)){
        c2->cd();
        Calibration(h,mean,sigma,a,b,0);
        cout << i+1 << " " << j+1 << " " << 0 << endl;
        ofile << token << " "  << b << " " << a << endl;
      }
      else 
        ofile << token << " " << 0 << " " << 1 << endl;
    }
    ofile.close();
  }
}

/////////////////////////////////
bool Finder(TH1F* h, Double_t *mean, Double_t *sigma)
{

  /////////////////////////////////////////////////
  //                                             //
  //	           ALPHA  FINDER                   //
  //                                             //
  /////////////////////////////////////////////////

  for(int k=0; k<m_NumberOfPeak; k++)
  {
    mean[k]=0;
    sigma[k]=0;
  }

  Double_t resolsig=3;
  Double_t  resolsigTSpec=1;
  Double_t seuil=0.2;

  //////// Peak finder

  TSpectrum *s = new TSpectrum(m_NumberOfPeak,resolsigTSpec);

  Int_t nfound = s->Search(h,resolsig,"new",seuil);
  Double_t *xpeaks = s->GetPositionX();

  /// Sort in growing order the array

  if(nfound>1)
  {
    for(Int_t p=0;p<nfound;p++)
    {
      for(Int_t i=0;i<nfound-1;i++)
      {
        if(xpeaks[i]>xpeaks[i+1])
        {
          Double_t varia=xpeaks[i];
          xpeaks[i]=xpeaks[i+1];
          xpeaks[i+1]=varia;
        }
      }
    }
  }

  Double_t linf=0, lsup=0;

  if(nfound == m_NumberOfPeak)
  {
    cout << "Number of peaks found matches number of expected peaks" << endl;
    for (Int_t p=0;p<nfound;p++)
    {

      linf = xpeaks[p]-5;
      lsup = xpeaks[p]+5;

      TF1 *gauss = new TF1("gauss","gaus",linf,lsup);
      h->Fit(gauss,"RQ");
      mean[p] = gauss->GetParameter(1);
      error_par[p] = gauss->GetParError(1);
      sigma[p]= gauss->GetParameter(2);

      //sigma_fit->Fill(gauss->GetParameter(2));
    }
  }

  if(nfound!=m_NumberOfPeak)
  {
    for (Int_t p=0;p<m_NumberOfPeak;p++)
    {
      cout << "Warning, Number of peak different of 5!!! " << h->GetName() << " nfound = " << nfound << endl;
      mean[p]=-1;
      sigma[p]=-1;
      return false ;
    }
  }

  return true ;
}


/////////////////////////////////
Double_t Calibration(TH1F* hist, Double_t* mean, Double_t* sigma, Double_t &a , Double_t &b, double pedestal)
{
  TGraphErrors *gr = new TGraphErrors();

  for (unsigned int p = 0; p < m_NumberOfPeak; p++) {
    cout << p << " " << mean[p] << " " << pedestal << " " << Source_E[p] << endl;
    errors[p] = 0.0001;
    gr->SetPoint(p, mean[p]-pedestal, Source_E[p]);
    //gr->SetPoint(p, mean[p], Source_E[p]);
    gr->SetPointError(p, error_par[p], errors[p]);
    gr->SetPointError(p, sigma[p], errors[p]);

  }

  //gr->SetPoint(3,4095,60);
  //gr->SetPointError(3,5,1);

  gr->SetMarkerStyle(8);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);
  gr->SetMarkerSize(1);
  gr->Fit("pol1", "q");

  chi2=gr->GetFunction("pol1")->GetChisquare()/gr->GetFunction("pol1")->GetNDF();

  a = gr->GetFunction("pol1")->GetParameter(1);
  b =  gr->GetFunction("pol1")->GetParameter(0);
  err_a = gr->GetFunction("pol1")->GetParError(1);
  err_b = gr->GetFunction("pol1")->GetParError(0);

  cout << "p0= " << b << endl;
  cout << "p1= " << a << endl;
  cout << "Energy value at channel=4095 -> " << b + a*4095 << " MeV" << endl;

  /*TF1 *f1 = new TF1("f2","[0]+[1]*x",0,16384);
    f1->SetParameter(0,gr->GetFunction("pol1")->GetParameter(0));
    f1->SetParameter(0,gr->GetFunction("pol1")->GetParameter(1));
    f1->SetLineColor(4);
    f1->SetLineWidth(2);*/

  gr->Draw("ap");
  //f1->Draw("lsame");

  for (Int_t p = 0; p < m_NumberOfPeak; p++) {
    residue[p] = Source_E[p] - (a*mean[p] + b);
    //error_res[p] = sqrt( pow(a,2)*pow(error_par[p],2) );
    error_res[p] = sqrt( pow(a,2)*pow(error_par[p],2) + pow(err_a,2)*pow(mean[p],2));
    Resolution = a*sigma[p];
  }


  TString token;

  //cout << token << " " << b << " " << a << endl;
  // look at the dispersion around Pedestals
  dispersion = Pedestal + b/a ;


  //delete f1;
  return a ;

}

/////////////////////////////
void DefineSource()
{
  /// Information used in the summary
  Source_Number_Peak = m_NumberOfPeak;
  Source_isotope = new TString[Source_Number_Peak] ;
  Source_E = new Double_t[Source_Number_Peak] ;
  Source_Sig = new Double_t[Source_Number_Peak] ;

  // Subtracting eloss from 0.5 um of Al
  Source_isotope[0]="ThSource"; Source_E[0]   = 5.15-0.0819 ; Source_Sig[0] = 0.00014 ;
  Source_isotope[1]="ThSource"; Source_E[1]   = 5.48-0.0784 ; Source_Sig[1] = 0.00014 ;
  Source_isotope[2]="ThSource"; Source_E[2]   = 5.80-0.0749  ; Source_Sig[2] = 0.00014 ;

  // Subtracting eloss from 3.5 um of Al
  //Source_isotope[0]="ThSource"; Source_E[0]   = 5.15-0.586 ; Source_Sig[0] = 0.00014 ;
  //Source_isotope[1]="ThSource"; Source_E[1]   = 5.48-0.565 ; Source_Sig[1] = 0.00014 ;
  //Source_isotope[2]="ThSource"; Source_E[2]   = 5.80-0.542  ; Source_Sig[2] = 0.00014 ;


  cout << "/**** Source characteristics ****/" << endl;
  for(int i=0; i<m_NumberOfPeak; i++){
    cout << Source_E[i] << " MeV" << endl;
  }

  return;
}
