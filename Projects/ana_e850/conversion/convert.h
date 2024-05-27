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
Int_t TMW1_XRawM;
UShort_t TMW1_XRawNr[92];
UShort_t TMW1_XRaw[92];
Int_t TMW1_YRawM;
UShort_t TMW1_YRawNr[92];
UShort_t TMW1_YRaw[92];

// TW2 //
Int_t TMW2_XRawM;
UShort_t TMW2_XRawNr[92];
UShort_t TMW2_XRaw[92];
Int_t TMW2_YRawM;
UShort_t TMW2_YRawNr[92];
UShort_t TMW2_YRaw[92];

// FPMW0 //
Int_t FPMW0_XRawM;
UShort_t FPMW0_XRawNr[992];
UShort_t FPMW0_XRaw[992];
Int_t FPMW0_YRawM;
UShort_t FPMW0_YRawNr[160];
UShort_t FPMW0_YRaw[160];

// FPMW1 //
Int_t FPMW1_XRawM;
UShort_t FPMW1_XRawNr[992];
UShort_t FPMW1_XRaw[992];
Int_t FPMW1_YRawM;
UShort_t FPMW1_YRawNr[160];
UShort_t FPMW1_YRaw[160];

// Exogam //
Int_t Inner6MVM;
Float_t Inner6MV[12];
UShort_t Inner6MVN[12];
ULong64_t Inner6MVTS[12];

void InitInputTree();
void InitOutputTree();
void Clear();
