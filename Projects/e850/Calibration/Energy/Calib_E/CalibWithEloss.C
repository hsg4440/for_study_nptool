int number_of_channels = 57;
int number_of_detectors = 8;
TChain* chain;
TH2F* h[8];
TCutG *cut[8];


///////////////////////////////////////////////
void CalibWithEloss()
{
  // Energy loss table 
  NPL::EnergyLoss C12Al = NPL::EnergyLoss("EnergyLossTable/C12_Al.G4table","G4Table",100);
  NPL::EnergyLoss C12Si = NPL::EnergyLoss("EnergyLossTable/C12_Si.G4table","G4Table",100);

  // Input file
  chain = new TChain("PhysicsTree");
  chain->Add("../../../../../Outputs/Analysis/run_57_nocalib.root");

  // cuts et histo
  TFile *fcut=new TFile("cut/cut_file.root","read");
  
  // Output file
  TFile * ofile = new TFile("backE_histo_file_run57_eloss.root","recreate");
  
  for(int k=0; k<number_of_detectors; k++){
    TString histo_name = Form("backE_h_det%i",k+1);
    h[k] = new TH2F(histo_name,histo_name,1000,0,60000,500,0,200);
    cut[k]=(TCutG*)fcut->Get(Form("cutC_det%i",k+1));
  }

  // reaction
  NPL::Reaction *r1 = new NPL::Reaction("238U(12C,12C)238U@1417");
  TGraph *g1 = r1->GetKinematicLine3();

  // TPISTAPhysics
  TPISTAPhysics* pista = new TPISTAPhysics();
  chain->SetBranchAddress("PISTA",&pista);
  
  // branche vamos
  ULong64_t fTS_TMW;
  chain->SetBranchStatus("fTS_TMW",true);
  chain->SetBranchAddress("fTS_TMW",&fTS_TMW);
 
  double thickness[8] = {106,108,108,105,106,108,108,106};

  int nentries = chain->GetEntries();
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    

    if(i%10000==0){
      cout << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush; 
    }

    int mult = pista->EventMultiplicity;
    if(mult==1){
      int det = pista->DetectorNumber[0];
      int stripE = pista->E_StripNbr[0];
      int stripDE = pista->DE_StripNbr[0];
      double E_strip = pista->E[0];
      double DE_strip = pista->DE[0];
      double E_back = pista->back_E[0];
      TVector3 PosPista = TVector3(pista->PosX[0],pista->PosY[0],pista->PosZ[0]);
      TVector3 PosTarget = TVector3(0,0,0);
      TVector3 pos = PosPista - PosTarget;
      double theta  = pos.Theta()/TMath::Pi()*180;

      double Elab_th = C12Si.EvaluateEnergyFromDeltaE(DE_strip,thickness[det-1]*micrometer,0,10,200,0.05);
      // E residual evolution //
      // -1- Elab - Eloss in entrance dead layer of DE
      double Eres = C12Al.Slow(Elab_th,0.5*micrometer,0);
      // Subtraction from measured DE
      Eres = Eres - DE_strip;
      // -2- Elab - Eloss in back dead layer of DE
      Eres = C12Al.Slow(Eres,0.5*micrometer,0);
      // -3- Elab - Eloss in entrance dead layer of E
      Eres = C12Al.Slow(Eres,0.5*micrometer,0);
      
      int condition = cut[det-1]->IsInside(E_back,DE_strip);
      
      if(det>0 && det<9 && stripE>0 && stripE<58 && stripDE>0 && stripDE<92 && E_back>18000 && fTS_TMW==0 && condition){
        h[det-1]->Fill(E_back,Eres);
      }
    }
  }



  TF1* f1 = new TF1("f1","[0] + [1]*x",0,60000);
  //f1-> SetParLimits(0,0,10);
  //f1-> SetParLimits(1,0.002,0.003);
  double p0, p1;
  
  ofstream ofile1;
  string ofilename1 = "PISTA_BACK_E_eloss.cal";
  ofile1.open(ofilename1.c_str());
  
  for(int i=0; i<number_of_detectors; i++){
    

    TString token = Form("PISTA_T%i_BACK_E",i+1); //PISTA_T%i_BACK_DE  PISTA_T%i_backE_ENERGY

    int N = h[i]->GetEntries();
    if(N>0){
      //h[i]->SetMinimum(2);
      h[i]->Fit("f1","qr");
      p0 = f1->GetParameter(0);
      p1 = f1->GetParameter(1);
    }
    else{
      p0 = 0;
      p1 = 1;
    }

    ofile1 << token << " " << p0 << " " << p1 << endl;
  }
  ofile1.close();

  ofile->Write();
  ofile->Close();


}

