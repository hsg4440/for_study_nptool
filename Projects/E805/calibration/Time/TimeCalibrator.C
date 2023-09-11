#include<stdlib.h>
#include<stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TString.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;
// #define TELESCOPE_START 1	
// #define TELESCOPE_END 4
#define STRIP_START	1
#define STRIP_END 	128


Int_t CurrentTelescope = 0 ;
Int_t CurrentStrip     = 0 ;


void TimeGene(TH1F*, TString, Int_t, string, int);
Int_t Finder(TH1F*, Double_t*, Double_t*, Int_t);
TF1* Calib(TGraphErrors*,Double_t*, Double_t*, Int_t);
Double_t pol2 (Double_t *,Double_t *);
		
ofstream peaks_file, calib_file, dispersion_file,calib_online_file; 
	
void TimeCalibrator(const char* frun = "run_XXXX", TString xy = "X",Int_t npeak= 32,string fDet="MUST2", int TELESCOPE_START=1, int TELESCOPE_END=4)//49
{
	// open the ROOT file to process
	TString path  = "./Histograms/";
	TString inFileName = frun;
  TFile *inFile;

  if (fDet=="MUST2"){
    inFile = new TFile(path + inFileName +"_RawMUST2Histos.root");
  }
  else if (fDet=="Mugast"){
    inFile = new TFile(path + inFileName +"_RawMugastHistos.root");
  }  


  TString str;
  TString str1;
  TString str2;
  TString strbuff;
  TString strbuff2;
  TString fname;
  TString fname2;

	TString hname;
	
	TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->cd();
	//c1->Divide(16,8);

	for (int det = TELESCOPE_START-1; det < TELESCOPE_END ; det++){
			ostringstream number ;
			number << det+1  ;
			CurrentTelescope = det+1 ;
		if (xy == "X"){ 
			//////// Input Files ///////////
			str = "hMM"+number.str()+"_STRX_T";
			/////// Output Files ///////////
			if (fDet=="MUST2"){
        str1 = "Cal_Str_X_T_MM"+number.str()+"_"+frun;
      }
      else if (fDet=="Mugast"){
        str1 = "Cal_Str_X_T_MG"+number.str()+"_"+frun;
      }
    
    }	 
		
    else if (xy == "Y"){ 	
			//////// Input Files ///////////
			str = "hMM"+number.str()+"_STRY_T";
			/////// Output Files ///////////
      if (fDet=="MUST2"){
        str1 = "Cal_Str_Y_T_MM"+number.str()+"_"+frun;
      }
      else if (fDet=="Mugast"){
        str1 = "Cal_Str_Y_T_MG"+number.str()+"_"+frun;
      }	
		}	 
		else {cout<< "ERROR FOR X or Y PARAMETER"<< endl;}

		fname =  "./Coeff/" + str1 + ".peak";
		peaks_file.open(( (string)fname ).c_str());

		fname2 = "./Coeff/"+ str1 + ".cal";
		calib_file.open( ( (string)fname2 ).c_str() );
		
		for (Int_t j = STRIP_START-1; j < STRIP_END; j++)
		{
			///// Get the histogram of det i and strip j /////
			CurrentStrip=j+1;
			number.seekp(0);
			number << j+1;
			hname = str+number.str();
			TH1F *hist = (TH1F*) inFile->Get(((string)hname).c_str());
		
			TimeGene(	hist,
								xy,
								npeak, fDet, det);
		}
	peaks_file.close();
	calib_file.close();
	dispersion_file.close();
	
	}
	
}



void TimeGene(TH1F *hist, TString xy, Int_t npeak, string fDet, int DetNbr)
	{
	
  // hist->GetXaxis()->SetRangeUser(0,16384);
  //hist->GetXaxis()->SetRangeUser(250,8300);
  TH1F* hist_clone = (TH1F*)hist->Clone();

 // int Nbins1 = hist_clone->GetNbinsX();
 // cout<<"Nbins before "<<Nbins1<<endl;
 
  hist_clone->Rebin(16);
  int Nbins = hist_clone->GetNbinsX();
 // cout<<"Nbins after "<<Nbins<<endl;
  for(int i =0; i<Nbins; i++){
      if(hist_clone->GetBinCenter(i) > 8170 ) hist_clone->SetBinContent(i,0);
     // if (DetNbr!=10){//range change for the anular telescope 11
      if(hist_clone->GetBinCenter(i) < 3800 ) hist_clone->SetBinContent(i,0);
     // }
     // else {
      //   if(hist_clone->GetBinCenter(i) < 100 ) hist_clone->SetBinContent(i,0);
     // }  
 } 
//  hist_clone->Smooth(10);

  
  double Ecut = 0;
  int maxbin = hist_clone->GetMaximumBin(); 
  while(hist_clone->GetBinContent(maxbin) != 0){
      maxbin++;
  }
  double maxE =  hist_clone->GetBinCenter(maxbin);
  //cout << "maxbin : "<< maxbin << endl;
  //cout << "maxE : "<< maxE << endl;
/*
  Nbins = hist->GetNbinsX();
  for(int i =0; i<Nbins; i++){
      if(hist->GetBinCenter(i) < maxE ) hist->SetBinContent(i,0);
      if(hist->GetBinCenter(i) > 8180 ) hist->SetBinContent(i,0);
  } 
  */
	Double_t *mean = new Double_t[npeak];
  Double_t *sigma = new Double_t[npeak];

	Int_t nfound = Finder(hist_clone, mean, sigma, npeak);

  if (fDet=="MUST2"){	
    if (xy == "X"){
      peaks_file << "MUST2_T" << CurrentTelescope << "_Si_X" << CurrentStrip << "_T " << endl;
      for(int i=0; i<nfound; i++) {
        peaks_file << i<< "\t" << mean[i] <<  endl;
      }
    }
    else if (xy == "Y"){
      peaks_file << "MUST2_T" << CurrentTelescope << "_Si_Y" << CurrentStrip << "_T " << endl;
      for(int i=0; i<nfound; i++) {
        peaks_file << i<< "\t" << mean[i] <<  endl;
      }
    }
  }
  else if (fDet=="Mugast"){	
    if (xy == "X"){
      peaks_file << "Mugast_T" << CurrentTelescope << "_DSSD_X" << CurrentStrip << "_T " << endl;
      for(int i=0; i<nfound; i++) {
        peaks_file << i<< "\t" << mean[i] <<  endl;
      }
    }
    else if (xy == "Y"){
      peaks_file << fDet <<"_T" << CurrentTelescope << "_DSSD_Y" << CurrentStrip << "_T " << endl;
      for(int i=0; i<nfound; i++) {
        peaks_file << i<< "\t" << mean[i] <<  endl;
      }
    }
  }

	if(nfound<npeak) cout << "On "<< hist->GetName() << endl ;
	
	// Fit npeak
	TGraphErrors* gr_MM= new TGraphErrors(nfound);
	TF1* f = Calib(gr_MM, mean, sigma, nfound);
	
	Double_t a = f->GetParameter(0);
	Double_t b = f->GetParameter(1);
	Double_t c = f->GetParameter(2);
  
  if(fDet=="MUST2"){
    if (xy == "X")
      calib_file << fDet << "_T" << CurrentTelescope << "_Si_X" << CurrentStrip << "_T " << a << " " << b << "	" << c <<endl ;

    else if (xy == "Y")
      calib_file << fDet << "_T" << CurrentTelescope << "_Si_Y" << CurrentStrip << "_T " << a << " " << b << " " << c << endl ;
  }
  else if(fDet=="Mugast"){
    if (xy == "X")
      calib_file << fDet << "_T" << CurrentTelescope << "_DSSD_X" << CurrentStrip << "_T " << a << " " << b << "	" << c <<endl ;

    else if (xy == "Y")
      calib_file << fDet << "_T" << CurrentTelescope << "_DSSD_Y" << CurrentStrip << "_T " << a << " " << b << " " << c << endl ;
  }



}

Int_t Finder(TH1F *h, Double_t *mean, Double_t *sigma, Int_t npeak)
{

	/////////////////////////////////////////////////
	//																						 //
	//	           	PEAK  FINDER									 //
	//																					 	 //
	/////////////////////////////////////////////////


	for(int k=0; k<npeak; k++)
			{
		    	mean[k]=0;
		        sigma[k]=0;
			}

	//Double_t resolsig=5;
	Double_t resolsig=2;
	Float_t  resolsigTSpec=1;
	Double_t seuil=0.05;
	Int_t nmax=npeak+1;   // nombre de pics maximum a prendre

	//////// Peak finder

	TSpectrum *s = new TSpectrum(nmax,resolsigTSpec);

	Int_t nfound = s->Search(h,resolsig,"new",seuil);
 	
  Double_t *xpeaks = s->GetPositionX();
  cout << "NUMBER OF PEAKS" << " " << nfound << endl;
	/// remise dans l ordre du tableau ordre croissant

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
					cout <<"XPeak: " << xpeaks[i]<<endl;
				}	  
			}
		}
	}

	
	Float_t linf=0, lsup=0; 
	
	// if(nfound == npeak || nfound > npeak)
	// {
  	//struct FitResult{
	//	double mean;
	//	double sigma;
	//};
	//map<int, FitResult> results;

		for (Int_t p=0;p<nfound;p++)
		{   
		linf = xpeaks[p]-5;
		lsup = xpeaks[p]+5;

		TF1 *gauss = new TF1("gauss","gaus",linf,lsup);
		h->Fit(gauss,"RQ");
		mean[p] = gauss->GetParameter(1);
		sigma[p]= gauss->GetParameter(2);
		//if (gauss->GetParameter(2)<10) {
		//	results[p].mean = gauss->GetParameter(1); 
		//	results[p].sigma = gauss->GetParameter(2);
		//}
    // std::cout <<  gauss->GetParameter(1) << std::endl;
    // std::cout <<  gauss->GetParameter(2)<< std::endl;
		 }
	// }
	//nfound=results.size();
	//int cnt = 0;
	//for (const auto & result: results){
	//	mean[cnt]=result.second.mean;
	//	sigma[cnt]=result.second.sigma;
	//	cnt++;	
	//	cout <<"Mean: " << mean[cnt]<<endl;
	//}
	
	


	if(nfound<npeak)
	{
	 cout << "attention, nombre de pics ("<<nfound<<") different de "<<npeak<<" !!!" ;
		/*
		for (Int_t p=0;p<npeak;p++)
		{
			mean[p]=-1;
			sigma[p]=-1;
			//return false ;
		}
		*/
	 
	}
	
	return nfound ;
}

TF1* Calib(TGraphErrors *gr, Double_t *mean, Double_t *sigma, Int_t npeak)
{  
	Double_t energy[npeak];
  //Double_t errors[npeak];
  
	
	for(int i = 0 ; i < npeak ; i ++)
	{
	energy[i] = (npeak-i)*20;
	}
			
	
	for (Int_t p = 0; p < npeak; p++) {
		gr->SetPoint(p, mean[p], energy[p]);
		gr->SetPointError(p, 0, 0);
		//gr->SetPointError(p, sigma[p], errors[p]);    
	}
	
	TF1 *f1 = new TF1("f1",pol2,mean[0],mean[npeak-1], npeak);
	//TF1 *f1 = new TF1("f1",Form("[0]*(x-%f)",Pedestal));
	gr->Fit("f1", "Q" );
	// TCanvas*c2 = new TCanvas();
	// c2->cd();
	gr->Draw("AP");
  // gPad->WaitPrimitive();

	//Double_t a = f1 -> GetParameter(0);
	
	return f1 ;
}

Double_t pol2 (Double_t *x,Double_t *par){

	Double_t result = 0;
	result = par[0]+ x[0]*par[1]+ x[0]*x[0]*par[2];
	return result;
}


	
