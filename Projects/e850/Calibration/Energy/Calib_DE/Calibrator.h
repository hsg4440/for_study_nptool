TFile* InputFile;


unsigned int m_NumberOfPeak=3;
unsigned int m_NumberOfDetectors=20;


Int_t Source_Number_Peak;
TString *Source_isotope;
Double_t *Source_E;
Double_t *Source_Sig;

vector<TH1F*> h;

bool Pedestals_Aligned = false;

// Calibration Coefficient
ofstream output_file;
ofstream output_file_ped;
Double_t Pedestal;
Double_t Pedestal_err;
Double_t a ;
Double_t b ;
Double_t err_a ;
Double_t err_b ;
TH2F* h_chi2;
Double_t dispersion;
Double_t chi2;
Double_t chi2_min = 100.;

Double_t* mean      = new Double_t[m_NumberOfPeak];
Double_t* sigma     = new Double_t[m_NumberOfPeak];
Double_t* error_par = new Double_t[m_NumberOfPeak];
Double_t* errors    = new Double_t[m_NumberOfPeak];
Double_t* residue   = new Double_t[m_NumberOfPeak];
Double_t* error_res = new Double_t[m_NumberOfPeak];

Double_t Resolution;

void OpenRootFile();
void DefineSource();
bool Finder(TH1F *h, Double_t *mean, Double_t *sigma);
Double_t Calibration(TH1F* hist, Double_t* mean, Double_t* sigma, Double_t &a , Double_t &b, double pedestal=0);
Double_t Pedestals(TH1F *hist);

