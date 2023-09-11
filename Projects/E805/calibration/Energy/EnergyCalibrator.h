// STL
#include<stdlib.h>
#include<stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
using namespace std;

// Root
#include "TString.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLine.h"

//NPTool
//#include "/home/mugast/analysis/nptool/NPLib/Physics/NPEnergyLoss.h"
//#include "/home/mugast/analysis/nptool/NPLib/Core/NPSystemOfUnits.h"

#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Physics/NPEnergyLoss.h"
#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Core/NPSystemOfUnits.h"

using namespace NPUNITS;
using namespace NPL;

/// DEFINING GLOBAL VARIABLE

double mean_extrapolation=100;

/// Parameter used in the macro
// Micrometer
Double_t AlThickness;
Double_t SiThickness;
EnergyLoss EL_Al("./EnergyLossTable/alpha_Al.G4table" , "G4Table", 100) ;
EnergyLoss EL_Si("./EnergyLossTable/alpha_Si.G4table" , "G4Table", 100) ;
// Information about the calibration condition (use Latex marks-up)

const TString xy                  = "Y";
const TString Experiment          = "MUGAST_LISE23";
const TString Run_Period          = "16/02/23"; 
const TString Operator            = "Valerian and Hugo";
const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
const TString Comment             = "MUST2";
const char* frun = "run_19";

// const TString Experiment          = "e793s";
// const TString Run_Period          = "08/03/21";
// const TString Operator            = "Franco";
// const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
// const TString Comment             = "Mugast";

// const TString Experiment          = "Test";
// const TString Run_Period          = "27/02/19";
// const TString Operator            = "VALERIAN and CYRIL";
// const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
// const TString Comment             = "Mugast";

/* const TString Experiment          = "E744"; */
/* const TString Run_Period          = "Mai 2018"; */
/* const TString Operator            = "Valerian"; */
/* const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm"; */
/* const TString Comment             = "T5 with 100V on DSSD"; */
/* const char* frun = "run_0196"; */
/*
   const TString Experiment          = "E744";
   const TString Run_Period          = "June 2018";
   const TString Operator            = "Valerian";
   const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
   const TString Comment             = "Y T1-T5";
   const char* frun = "run_3001";
   */

//const TString Experiment          = "RIBF57";
//const TString Run_Period          = "April 2010, Riken BigRIPS, Run 3";
//const TString Operator            = "Adrien MATTA";
//const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const TString Comment             = "Source at 0$^{\\circ}$";
//const char* frun = "RIBF57_runx";

//const TString Experiment          = "e569s";
//const TString Run_Period          = "july 2009, Ganil VAMOS, Run 19-21";
//const TString Operator            = "Adrien MATTA";
//const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const TString Comment             = "Source at 0$^{\\circ}$ facing Telescope 1,2,3,4";
//const char*   frun                = "e569s_run_19_21";

//const TString Experiment          = "e530";
//const TString Run_Period          = "march 2009, Ganil LISE, Run 3";
//const TString Operator            = "Adrien MATTA";
//const TString Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const TString Comment             = "Source at 0$^{\\circ}$ facing Telescope 1,2,3,4";
//const char* frun = "e530_run_0003";

int Telescope_Number=0;
const int Strip_Start=1;
const int Strip_End=128;

// choosing a method for the fit
const TString method = "ZeroExtrapolation" ;
//const TString method = "ZeroForce";
const bool RefitWithSatellite = false;
const bool Pedestals_Aligned = true;   
Int_t CurrentTelescope = 0;
Int_t CurrentStrip     = 0;
TString folder;
TString main_name;
TCanvas* Tsummary;
TCanvas* Buffer;
double sigma_fit_centroid = 0.;
double sigma_fit_sigma = 0.;

map<int,string> BadStrip;

// Defining the array used after (order needs to be diffent for X and Y )
Int_t NumberOfIsotope;

// Source original value
Int_t Source_Number_Peak;
TString*  Source_isotope;
Double_t* Source_branching_ratio;
Double_t* Source_E;
Double_t* Source_Sig;

// Source corrected value
Double_t* energyX;
Double_t* errorsX;
Double_t* energyY;
Double_t* errorsY;

// Calibration Coefficient
Double_t a ;
Double_t b ;

Double_t* mean      = new Double_t[3]; 
Double_t* sigma     = new Double_t[3];
Double_t* error_par = new Double_t[3];

TGraph* ZeroDispersion ;
ofstream peaks_file, calib_file, dispersion_file , calib_online_file, latex_file;; 

TH1F* sigma_fit  ;
TH1F* Dispersion ;
TGraph* coeff_a ;
TGraph* coeff_b ;

TFile *inFile;

/// Function Header
void AutoCalibration(int first_telescope ,int last_telescope, std::string fDet);
void EnergyCalibrator(std::string);
Double_t Pedestals(TH1F *);
void Alpha(TH1F*, TString, Double_t,std::string);
bool Finder(TH1F*, TString , Double_t*, Double_t*);
Double_t Calib_ZeroForceMethod(string ,TGraphErrors*,float, Double_t*, Double_t*,std::string);
Double_t Calib_ZeroExtrapolationMethod(TH1F* hist ,string ,TGraphErrors*,float, Double_t*, Double_t*, Double_t &a , Double_t &b,std::string);
void LatexSummaryHeader(TString xy,std::string);
void LatexSummaryEnder();
void LatexSummaryTelescope();
void DefineSource(TString sourceName="3 alphas");


void Find_Satellites(TH1F *h);
Double_t source_Pu(Double_t *x, Double_t *par);
Double_t source_Am(Double_t *x, Double_t *par);
Double_t source_Cm(Double_t *x, Double_t *par);

