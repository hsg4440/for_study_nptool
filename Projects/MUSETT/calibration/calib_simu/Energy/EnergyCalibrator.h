// STL
#include<stdlib.h>
#include<stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

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
#include "NPEnergyLoss.h"
#include "NPSystemOfUnits.h"


using namespace NPUNITS;
using namespace NPL;

/// DEFINING GLOBAL VARIABLE

double mean_extrapolation=100;

/// Parameter used in the macro
// Micrometer
Double_t AlThickness;
Double_t SiThickness;
Double_t CThickness;
EnergyLoss EL_Al("./EnergyLossTable/alpha_Al.G4table" , "G4Table", 100) ;
EnergyLoss EL_Si("./EnergyLossTable/alpha_Si.G4table" , "G4Table", 100) ;
EnergyLoss EL_C("./EnergyLossTable/alpha_C.G4table" , "G4Table", 100) ;
// Information about the calibration condition (use Latex marks-up)

const TString xy                  = "Y";
const TString Experiment          = "double_alpha";
const TString Run_Period          = "21/06/2023"; 
const TString Operator            = "Hugo";
const TString Source              = "3 alpha peaks $^{222}$Ra, $^{218}$Rn, $^{214}$Po";
const TString Comment             = "MUSETT";
const char* frun = "run_0274";




int Telescope_Number=0;
const int Strip_Start=0;
const int Strip_End=127;

const int Xmin = 8700;
const int Xmax = 9500;

const int Ymin = 6800;
const int Ymax = 8000;

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
void EnergyCalibrator();
Double_t Pedestals(TH1F *);
void Alpha(TH1F*, TString, Double_t);
bool Finder(TH1F*, TString , Double_t*, Double_t*);
Double_t Calib_ZeroForceMethod(string ,TGraphErrors*,float, Double_t*, Double_t*);
Double_t Calib_ZeroExtrapolationMethod(TH1F* hist ,string ,TGraphErrors*,float, Double_t*, Double_t*, Double_t &a , Double_t &b);
void LatexSummaryHeader(TString xy);
void LatexSummaryEnder();
void LatexSummaryTelescope();
void DefineSource(TString sourceName="3 alphas");


void Find_Satellites(TH1F *h);
Double_t source_Pu(Double_t *x, Double_t *par);
Double_t source_Am(Double_t *x, Double_t *par);
Double_t source_Cm(Double_t *x, Double_t *par);

