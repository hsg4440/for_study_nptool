TPISTAData* m_pista;
TFPMWData* m_fpmw;
TICData* m_ic;
//TFile* ofile;
TChain* input_tree;
//TTree* output_tree;

// TAC //
float T_TMW0_FPMW0;
float T_TMW0_FPMW1;
float T_TMW1_FPMW0;
float T_TMW1_FPMW1;

// Exogam //
int Exo_Mult;
vector<float> Exo_Energy;
vector<int> Exo_Crystal;

// Pat
UShort_t fFPMWPat_0RawNr[20];
Int_t fFPMWPat_0RawM;
ULong64_t TMWPat_0TS;
ULong64_t fVAMOS_TS_sec;
ULong64_t PISTA_TS;
ULong64_t fPISTA_TS_sec;

// IC //
float IC[11];

// TW1 //
int TMW1_XVM;
float TMW1_XV[92];
UShort_t TMW1_XVN[92];
int TMW1_YVM;
float TMW1_YV[92];
UShort_t TMW1_YVN[92];

// TMW2 //
int TMW2_XVM;
float TMW2_XV[92];
UShort_t TMW2_XVN[92];
int TMW2_YVM;
float TMW2_YV[92];
UShort_t TMW2_YVN[92];

// FPMW0 //
int FPMW0_XVM;
float FPMW0_XV[992];
UShort_t FPMW0_XVN[992];
int FPMW0_YVM;
float FPMW0_YV[160];
UShort_t FPMW0_YVN[160];

// FPMW1 //
int FPMW1_XVM;
float FPMW1_XV[992];
UShort_t FPMW1_XVN[992];
int FPMW1_YVM;
float FPMW1_YV[160];
UShort_t FPMW1_YVN[160];

// Exogam //
Int_t Inner6MVM;
Float_t Inner6MV[12];
UShort_t Inner6MVN[12];
ULong64_t Inner6MVTS[12];

void InitInputTree();
void InitOutputTree();
void Clear();
