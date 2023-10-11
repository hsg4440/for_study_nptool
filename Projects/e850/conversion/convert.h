TPISTAData* m_pista;
TFPMWData* m_fpmw;
TICData* m_ic;
TChain* input_tree;
TTree* output_tree;

// TAC //
float T_TMW0_FPMW0;
float T_TMW0_FPMW1;
float T_TMW1_FPMW0;
float T_TMW1_FPMW1;

// IC //
float IC[11];

// TW1 //
int TMW1_XVM;
float TMW1_XV[92];
short TMW1_XVN[92];
int TMW1_YVM;
float TMW1_YV[92];
short TMW1_YVN[92];

// TMW2 //
int TMW2_XVM;
float TMW2_XV[92];
short TMW2_XVN[92];
int TMW2_YVM;
float TMW2_YV[92];
short TMW2_YVN[92];

// FPMW0 //
int FPMW0_XVM;
float FPMW0_XV[992];
short FPMW0_XVN[992];
int FPMW0_YVM;
float FPMW0_YV[160];
short FPMW0_YVN[160];

// FPMW1 //
int FPMW1_XVM;
float FPMW1_XV[992];
short FPMW1_XVN[992];
int FPMW1_YVM;
float FPMW1_YV[160];
short FPMW1_YVN[160];

void InitInputTree();
void InitOutputTree();
void Clear();
