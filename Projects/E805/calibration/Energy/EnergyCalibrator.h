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

#include "../../../../NPLib/Physics/NPEnergyLoss.h"
#include "../../../../NPLib/Core/NPSystemOfUnits.h"

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

const std::string xy                  = "X";
const std::string Experiment          = "MUGAST_LISE23";
const std::string Run_Period          = "16/02/23"; 
const std::string Operator            = "Valerian and Hugo";
const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
const std::string Comment             = "MUST2";
const char* frun = "";

// const std::string Experiment          = "e793s";
// const std::string Run_Period          = "08/03/21";
// const std::string Operator            = "Franco";
// const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
// const std::string Comment             = "Mugast";

// const std::string Experiment          = "Test";
// const std::string Run_Period          = "27/02/19";
// const std::string Operator            = "VALERIAN and CYRIL";
// const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
// const std::string Comment             = "Mugast";

/* const std::string Experiment          = "E744"; */
/* const std::string Run_Period          = "Mai 2018"; */
/* const std::string Operator            = "Valerian"; */
/* const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm"; */
/* const std::string Comment             = "T5 with 100V on DSSD"; */
/* const char* frun = "run_0196"; */
/*
   const std::string Experiment          = "E744";
   const std::string Run_Period          = "June 2018";
   const std::string Operator            = "Valerian";
   const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
   const std::string Comment             = "Y T1-T5";
   const char* frun = "run_3001";
   */

//const std::string Experiment          = "RIBF57";
//const std::string Run_Period          = "April 2010, Riken BigRIPS, Run 3";
//const std::string Operator            = "Adrien MATTA";
//const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const std::string Comment             = "Source at 0$^{\\circ}$";
//const char* frun = "RIBF57_runx";

//const std::string Experiment          = "e569s";
//const std::string Run_Period          = "july 2009, Ganil VAMOS, Run 19-21";
//const std::string Operator            = "Adrien MATTA";
//const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const std::string Comment             = "Source at 0$^{\\circ}$ facing Telescope 1,2,3,4";
//const char*   frun                = "e569s_run_19_21";

//const std::string Experiment          = "e530";
//const std::string Run_Period          = "march 2009, Ganil LISE, Run 3";
//const std::string Operator            = "Adrien MATTA";
//const std::string Source              = "3 alpha peaks $^{239}$Pu, $^{241}$Am, $^{244}$Cm";
//const std::string Comment             = "Source at 0$^{\\circ}$ facing Telescope 1,2,3,4";
//const char* frun = "e530_run_0003";

int Telescope_Number=0;
const int Strip_Start=1;
const int Strip_End=128;

// choosing a method for the fit
const std::string method = "ZeroExtrapolation" ;
//const std::string method = "ZeroForce";
const bool RefitWithSatellite = false;
const bool Pedestals_Aligned = true;   
Int_t CurrentTelescope = 0;
Int_t CurrentStrip     = 0;
std::string folder;
std::string main_name;
TCanvas* Tsummary;
TCanvas* Buffer;
double sigma_fit_centroid = 0.;
double sigma_fit_sigma = 0.;

map<int,string> BadStrip;

// Defining the array used after (order needs to be diffent for X and Y )
Int_t NumberOfIsotope;

// Source original value
Int_t Source_Number_Peak;
std::string*  Source_isotope;
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
void Alpha(TH1F*, std::string, Double_t,std::string);
bool Finder(TH1F*, std::string , Double_t*, Double_t*);
Double_t Calib_ZeroForceMethod(string ,TGraphErrors*,float, Double_t*, Double_t*,std::string);
Double_t Calib_ZeroExtrapolationMethod(TH1F* hist ,string ,TGraphErrors*,float, Double_t*, Double_t*, Double_t &a , Double_t &b,std::string);
void LatexSummaryHeader(std::string xy,std::string);
void LatexSummaryEnder();
void LatexSummaryTelescope();
void DefineSource(std::string sourceName="3 alphas");


void Find_Satellites(TH1F *h);
Double_t source_Pu(Double_t *x, Double_t *par);
Double_t source_Am(Double_t *x, Double_t *par);
Double_t source_Cm(Double_t *x, Double_t *par);

