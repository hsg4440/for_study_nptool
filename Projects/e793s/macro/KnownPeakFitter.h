#include "TMath.h"
#include "math.h"
#include <cmath>
#include "stdlib.h"

vector<double> means = { 0.000,
                           0.143,
                           0.279,
			   0.728,
			   0.967,
			   1.409,
			   1.978,
			   2.407,
			   2.909,
			   3.254,
			   3.601,
			   3.830,
			   //4.16,
			   //4.29,
			   //4.35
			   4.363
                           };
const int numPeaks = means.size(); 


vector<double> means_dt = { 0.000,
	                   0.587,
			   0.691,
			   //0.886,
			   //1.738,
                           1.944,
			   2.233, 
			   2.732,
			   3.377,//DOUBLET
                           //3.344,
                           //3.410,
			   4.297
                           };

const int numPeaks_dt = means_dt.size(); 


array<double,27> knowngammas = { 0.143,
					0.279,
					0.449,
					0.968,
					1.130,
					1.410,
					1.267,
					//0.575,
					1.013,
					1.838,
					1.981,
					1.000,
					2.412,
					2.767,
					2.518,
					3.325,
					2.878,
					3.605,
					2.839,
					2.734,
					3.522,
					3.076,
					3.875,
					0.834,
					3.325,
					3.77,
					4.037,
					4.364					
					};

/*
Double_t f_bg(Double_t *x, Double_t *par){
  // Flat bg [0] + semicircle [1]*sqrt(6.183^2 - (x-10.829)^2) 
  Float_t xx = x[0];
  Double_t f;
  Double_t a = TMath::Power(6.183,2);
  Double_t b = TMath::Power(xx-10.829,2);
  if(a > b){ f = par[0] + (par[1]*TMath::Sqrt(a-b)); }
  else{ f = par[0]; }
  return f;
}

Double_t f_peak(Double_t *x, Double_t *par){
  float xx = x[0];
  double f = (par[2]/(par[0]*sqrt(2*pi)))*exp(-0.5*pow((xx-par[1])/par[0],2));
  return f;
}
*/

string f_full_string(){
  string result = "gaus(0)";
  for(int pk=1; pk<numPeaks; pk++){
    result += "+gaus(" + to_string(3*pk) + ")"; 
  }
  return result;
}

Double_t f_full(Double_t *x, Double_t *par) {
  float xx = x[0];
  double result, norm;
  result = 0;

  // Add N peaks
  for(int pk=0; pk<numPeaks; pk++){
    result += (par[1+(pk*3)]*exp(-0.5*pow(((xx-par[2+(pk*3)]-par[0])/par[3+(pk*3)]),2))); //REGULAR GAUS
  }
  return result;
}

Double_t f_full_dt(Double_t *x, Double_t *par) {
  float xx = x[0];
  double result, norm;
  // Flat background
//  result = par[0];
  result = 0;
  // Add N peaks
  for(int pk=0; pk<numPeaks_dt; pk++){
    result += (par[3+(pk*3)]/(par[1+(pk*3)]*sqrt(2*pi)))
	      //* exp(-0.5*pow((xx-par[2+(pk*3)])/par[1+(pk*3)],2));
	      * exp(-0.5*pow((xx-par[2+(pk*3)]-par[0])/par[1+(pk*3)],2)); //added par 0 as shift in energy
  }
  return result;
}

vector<vector<double>> FitKnownPeaks_RtrnArry(TH1F* hist, double slideshift){
  double minFit=-1.0, maxFit=4.644; 
  double binWidth = hist->GetXaxis()->GetBinWidth(3);
  //double sigma = 0.14;

  double sigma = 0.0;// = 0.14;
  int numPeaks_temp = 0;
  vector<double> means_temp; 

  //Gradient on mean value
  double gradient;
  if(reactionName=="47K(d,p)"){
    sigma = 0.14;
    gradient = 0.989;
    numPeaks_temp = numPeaks;
    means_temp = means;
  } else if (reactionName=="47K(d,t)"){
    sigma = 0.18;
    gradient = 0.986;
    numPeaks_temp = numPeaks_dt;
    means_temp = means_dt;
  }

  cout << BOLDGREEN 
       << "FITTING WITH A REAL:OBSERVED Ex GRADIENT GRADIENT OF " 
       << gradient 
       << RESET 
       << endl;

  //hist->Sumw2();

  //Build individual peak fit functions
  string nameBase = "Peak ";
  string function = "gaus";//"[0]*exp(-0.5*pow(((x-[1])/[2]),2))";//REGULAR GAUS
  TF1 **allPeaks = new TF1*[numPeaks_temp];
  for(int i=0; i<numPeaks_temp; i++) {
    string nameHere = nameBase;
    nameHere +=to_string(i);
    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), minFit, maxFit);
    allPeaks[i]->SetLineColor(kBlack);  
    allPeaks[i]->SetLineStyle(7);  
    allPeaks[i]->SetNpx(500);
  } 

  //Subtract flat background equal to smallest bin in range
  hist->GetXaxis()->SetRange(hist->FindBin(-0.8), hist->FindBin(-0.0));
  double bgmin = hist->GetBinContent(hist->GetMinimumBin());
  hist->GetXaxis()->UnZoom();
  cout << "Subtracting background of " << bgmin << endl;
  for(int b=1; b<hist->GetNbinsX() ; b++){
      hist->SetBinContent(b,hist->GetBinContent(b)-bgmin);
  }

  //Build background function
  TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
  bg->SetLineColor(kGreen);
  bg->SetLineStyle(9);
  bg->SetParNames("Background");

  //Build total function
  string s_full = f_full_string();
  TF1 *full = new TF1("fitAllPeaks", s_full.c_str(), minFit, maxFit);
  full->SetLineColor(kRed);
  const int numParams = (numPeaks_temp*3);
  for(int i=0; i<numPeaks_temp; i++) {
    full->FixParameter((i*3)+1,means_temp.at(i)*gradient);
    full->FixParameter((i*3)+2,sigma);
    full->SetParameter((i*3)+0,1e1);
    full->SetParLimits((i*3)+0,0.0,1e5);
  }
  full->SetNpx(500);

  //Set and fix specific parameters
  for (int i=0; i<numPeaks; i++){
    if(means[i]>4.0){
      full->SetParameter((i*3)+2,0.23);
    }
    if(abs(means[i]-0.279)<0.05){
      full->FixParameter((i*3)+0,0.0);
    }
  }

  //Fit full function to histogram
  //hist->Fit(full, "RWQB", "", minFit, maxFit);
  hist->Fit(full, "RWQBL", "", minFit, maxFit); cout << GREEN << "FITTING WITH MAX LIKELYHOOD METHOD" << RESET << endl;

  hist->Draw();

  //Extract fitted variables, assign them to individual fits, and draw them
  const Double_t* finalPar = full->GetParameters();
  const Double_t* finalErr = full->GetParErrors();
  for (int i=0; i<numPeaks_temp; i++){
    allPeaks[i]->SetParameters(finalPar[0+(i*3)], finalPar[1+(i*3)], finalPar[2+(i*3)]);
  }
  bg->SetParameter(0,finalPar[0]);
  bg->Draw("SAME");
  full->Draw("SAME");

  //Make irresolvable doublet clear
  if(reactionName=="47K(d,p)"){
    for (int i=0; i<numPeaks_temp; i++){
      if(abs(means[i]-3.83)<0.1){
        cout << BOLDCYAN 
             << "PEAK #" << i 
             << " IS FITTING IRRESOLVABLE DOUBLET!" 
             << RESET 
             << endl;
        allPeaks[i]->SetLineColor(kCyan+1);
      }
      if(means[i]>4.0){
        cout << BOLDMAGENTA 
             << "PEAK #" << i 
             << " IS FITTING UNKNOWN ZONE WITH SIGMA = "
             << full->GetParameter((i*3)+2)
             << RESET 
             << endl;
        allPeaks[i]->SetLineColor(kMagenta);
      }
    }
  } else if(reactionName=="47K(d,t)"){
    for (int i=0; i<numPeaks_temp; i++){
      if(abs(means[i]-3.37)<0.1){
        cout << BOLDCYAN 
             << "PEAK #" << i 
             << " IS FITTING IRRESOLVABLE DOUBLET!" 
             << RESET 
             << endl;
        allPeaks[i]->SetLineColor(kCyan+1);
      }
    }
  }

  //Draw all peaks
  for (int i=0; i<numPeaks_temp; i++){
    allPeaks[i]->Draw("SAME");
  }

  //Write to screen
  cout << "===========================" << endl;
  cout << "== PEAK =========== AREA ==" << endl;
  
  //Loop over every peak
  vector<vector<double>> allpeaks;
  for(int i=0; i<numPeaks_temp; i++){
    //For AREA = HEIGHT * SIGMA * SQRT(2*PI)
    double A = (finalPar[0+(i*3)] * finalPar[2+(i*3)] * sqrt(2*pi)) / binWidth;
//    cout << "AREA = " << finalPar[0+(i*3)] << " * "
//	              << finalPar[2+(i*3)] << " * "
//		      << sqrt(2*pi)        << " = "
//		      << A                 << endl;

    double deltaA = A * sqrt( pow( finalErr[0+(i*3)]/finalPar[0+(i*3)] ,2) + pow( finalErr[2+(i*3)]/finalPar[2+(i*3)] ,2) );
    cout << "DELTAAREA = " 
	 << finalErr[0+(i*3)] << " / " 
	 << finalPar[0+(i*3)] << " and "
	 << finalErr[2+(i*3)] << " / " 
	 << finalPar[2+(i*3)] << " = " 
	 << deltaA << endl;

    cout << fixed << setprecision(3) 
	 << " #" << i << "  " 
	 << finalPar[(i*3)+1] << "\t" << setprecision(0)
	 << A << "\t+- " 
	 << deltaA << setprecision(3);
    cout << "    SQRT: " << sqrt(A) << endl;

    vector<double> onepeak; //energy, area and error for one peak
    onepeak.push_back(finalPar[(i*3)+1]);
    onepeak.push_back(A);
    onepeak.push_back(deltaA);
    allpeaks.push_back(onepeak);
    cout << "-------------" << endl;
  }

  return allpeaks;
}

////vector<vector<double>> FitKnownPeaks_dt_RtrnArry(TH1F* hist, double slideshift){
////  double minFit=-1.0, maxFit=8.0; 
////  double binWidth = hist->GetXaxis()->GetBinWidth(3);
////  double sigma = 0.18;
////
////  hist->Sumw2();
////
////  /* Construct flat BG to subtract */ 
////  /**
////  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
////  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
////  cout << " REMOVING FLAT BG OF 36 COUNTS!!!!" << endl;
////  double ConstBG = 36.0; double ErrBG = 1.0;
////  int xbins = hist->GetXaxis()->GetNbins();
////  double xmin = hist->GetXaxis()->GetXmin();
////  double xmax = hist->GetXaxis()->GetXmax();
////  TH1F *FlatBG = new TH1F("FlatBG","FlatBG", xbins, xmin, xmax);
////  for(int i=0; i<xbins;i++){
////    FlatBG->SetBinContent(i,ConstBG);
////    FlatBG->SetBinError(i,ErrBG);
////  }
////  hist->Add(FlatBG,-1);
////  **/
////
////  //Build individual peak fit functions
////  string nameBase = "Peak ";
////  string function = "([2]/([0]*sqrt(2*pi)))*exp(-0.5*pow((x-[1])/[0],2))";
////  TF1 **allPeaks = new TF1*[numPeaks_dt];
////  for(int i=0; i<numPeaks_dt; i++) {
////    string nameHere = nameBase;
////    nameHere +=to_string(i);
////
////    allPeaks[i] = new TF1(nameHere.c_str(), function.c_str(), minFit, maxFit);
////    //allPeaks[i] = new TF1(nameHere.c_str(), f_peak, -1, 5);
////    allPeaks[i]->SetLineColor(kBlack);  
////    allPeaks[i]->SetLineStyle(7);  
////    allPeaks[i]->SetParNames("Sigma", "Mean", "Area*BinWidth");
////  } 
////
////  //Subtract flat background equal to smallest bin in range
////  hist->GetXaxis()->SetRange(hist->FindBin(-0.8), hist->FindBin(-0.0));
////  double bgmin = hist->GetBinContent(hist->GetMinimumBin());
////  hist->GetXaxis()->UnZoom();
////  cout << "Subtracting background of " << bgmin << endl;
////  for(int b=1; b<hist->GetNbinsX() ; b++){
////      hist->SetBinContent(b,hist->GetBinContent(b)-bgmin);
////  }
////  
////  ////Subtract flat background equal to average bin in range
////  //double bgavg = 0;
////  //int bgcount = 0;
////  //for(int i = hist->FindBin(-0.8); i<hist->FindBin(0.1); i++){
////  //  bgavg += hist->GetBinContent(i);
////  //  bgcount++;
////  //}
////  //bgavg = bgavg/(double)bgcount;
////  //hist->GetXaxis()->UnZoom();
////  //cout << "Subtracting background of " << bgavg << endl;
////  //for(int b=1; b<hist->GetNbinsX() ; b++){
////  //    hist->SetBinContent(b,hist->GetBinContent(b)-bgavg);
////  //}
//// 
////
////  //Build background function
////  TF1 *bg = new TF1 ("bg","[0]",minFit, maxFit);
////  bg->SetLineColor(kGreen);
////  bg->SetLineStyle(9);
////  bg->SetParNames("Background");
////
////  //Build IMPROVED total function
////  TF1 *full = new TF1("fitAllPeaks", f_full_dt, minFit, maxFit, (int) 1+(3*numPeaks_dt));
////  full->SetLineColor(kRed);
////  const int numParams = (numPeaks_dt*3)+1;
////  for(int i=0; i<numPeaks_dt; i++) {
////    //full->FixParameter((i*3)+1,sigma);
////    //full->FixParameter((i*3)+2,means_dt.at(i));
////    //full->SetParameter((i*3)+3,1e1);
////    //full->SetParLimits((i*3)+3,0.0,1e5);
////    ////full->SetParLimits((i*3)+3,10.0,1e5);
////    //
////    //
////    //
////    //
////    full->FixParameter((i*3)+1,means.at(i)); //GRADIENT REMOVED FOR NOW
////    full->FixParameter((i*3)+2,sigma);
////    full->SetParameter((i*3)+0,1e1);
////    full->SetParLimits((i*3)+0,0.0,1e5);
////
////
////  }
////  //full->SetParLimits(0,0.,40.); /* FOR TOTAL SPECTRUM FITTING */
////  //full->SetParLimits(0,0.,10.); /* FOR ANGLE GATED FITTING */
////  //full->FixParameter(0,0.); /* FOR ANGLE GATED FITTING WITH BG SUBTRACTED */
////  
////  //Fit full function to histogram
////  //hist->Fit(full, "RWQB", "", minFit, maxFit);
////  hist->Fit(full, "RWQBL", "", minFit, maxFit); cout << GREEN << "FITTING WITH MAX LIKELYHOOD METHOD" << RESET << endl;
////  hist->Draw();
////
////  //Extract fitted variables, assign them to individual fits, and draw them
////  const Double_t* finalPar = full->GetParameters();
////  const Double_t* finalErr = full->GetParErrors();
////  for (int i=0; i<numPeaks_dt; i++){
////    allPeaks[i]->SetParameters(sigma, means_dt.at(i), finalPar[3+(i*3)]);
////  }
////  bg->SetParameter(0,finalPar[0]);
////  bg->Draw("SAME");
////  full->Draw("SAME");
////
////  for (int i=0; i<numPeaks_dt; i++){
////    allPeaks[i]->Draw("SAME");
////  }
////
//// /* Error propogation:
////  * (Abin) +- deltaAbin, B+-0 (no uncertainty)
////  * A = Abin/B
////  * deltaA/A = deltaAbin/Abin
////  * deltaA = A x deltaAbin/Abin
////  */
////
////  //Write to screen
////  cout << "===========================" << endl;
////  cout << "== PEAK =========== AREA ==" << endl;
////  
////  vector<vector<double>> allpeaks;
////  for(int i=0; i<numPeaks_dt; i++){
////    double A = finalPar[(i*3)+3]/binWidth;
////    double deltaA = A *  (finalErr[(i*3)+3]/finalPar[(i*3)+3]);
////
////    cout << fixed << setprecision(3) 
////	 << " #" << i << "  " 
////	 << finalPar[(i*3)+2] << "\t" << setprecision(0)
////	 << A << "\t+- " 
////	 << deltaA << setprecision(3)
////	 << endl;
////
////    vector<double> onepeak; //energy, area and error for one peak
////    onepeak.push_back(finalPar[(i*3)+2]);
////    onepeak.push_back(A);
////    onepeak.push_back(deltaA);
////    allpeaks.push_back(onepeak);
////  }
////  cout << " BG  " << full->GetParameter(0) 
////       << " +- " << full->GetParError(0) << endl;
////
////  return allpeaks;
////}

void FitKnownPeaks(TH1F* hist){
  //Shell function to call Rtrn_Arry without writing vector<vector<double>> to screen
  vector<vector<double>> shell = FitKnownPeaks_RtrnArry(hist,0.0);
}

void FitKnownPeaks_dt(TH1F* hist){
  //Shell function to call Rtrn_Arry without writing vector<vector<double>> to screen
  //vector<vector<double>> shell = FitKnownPeaks_dt_RtrnArry(hist,0.0);
  vector<vector<double>> shell = FitKnownPeaks_RtrnArry(hist,0.0);
}

