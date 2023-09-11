#include "TCutG.h"
#include "TFile.h"

using namespace std;

TCutG *cut[34];
TSpline3 *gSpline3[34];

int IdentifyZ(double,double);
int FindClosestZ(double, double);

void FillHistoPISTA(int Run)
{

  //  const int Nrun=10;
  const int Nrun=1;
  //  int RunNumber[Nrun]={29,41,42,43,44,45,49,50,51,52};
  int RunNumber[Nrun]={Run};
  //  int RunNumber[Nrun]={45,49,50,51};

  double Reference[34];
  for(int i=0;i<34;i++) Reference[i]=10000-i/34.*(10000-6539);

  TChain *chain=new TChain("PhysicsTree");
  for(int i=0;i<Nrun;i++) 
  {
    chain->Add(Form("../Outputs/Analysis/run_0%d.root",RunNumber[i]));
    cout << "RunNumber=" << RunNumber[i] << endl;
  }

  ///////// Cuts

  TFile *fcut=new TFile("cut_PISTA_PID.root");
  TCutG *cutBe10=(TCutG*)fcut->Get("cutBe10");
  TCutG *cutC12=(TCutG*)fcut->Get("cutC12");
  TCutG *cutC14=(TCutG*)fcut->Get("cutC14");
  TCutG *cutChio1=(TCutG*)fcut->Get("cutChio2");

  TSpline3 *fspline=(TSpline3*)fcut->Get("hChio1_DE_E_pfx");
  
  TFile *fcutChio=new TFile("cut_Chio_Indu.root");
  for(int i=0;i<34;i++) cut[i]=(TCutG*)fcutChio->Get(Form("Z%d",i+1));

  TFile *Fspline=new TFile("spline_Chio.root");
  for(int i=0;i<34;i++) gSpline3[i]=(TSpline3*)Fspline->Get(Form("fspline_%d",i+1));

  //////// Variables and branches

  double DeltaEcorr, Eres;
  ULong64_t fTS_TMW;
  double Ex240Pu, Ex236U, Ex238U, Ex240PuNoCorr;
  double Chio_DE, Chio_E;
  Float_t fTAC_MW1_PISTA;

  chain->SetBranchStatus("DeltaEcorr",true);
  chain->SetBranchAddress("DeltaEcorr",&DeltaEcorr);
  chain->SetBranchStatus("Eres",true);
  chain->SetBranchAddress("Eres",&Eres);
  chain->SetBranchStatus("fTS_TMW",true);
  chain->SetBranchAddress("fTS_TMW",&fTS_TMW);
  chain->SetBranchStatus("fTAC_MW1_PISTA",true);
  chain->SetBranchAddress("fTAC_MW1_PISTA",&fTAC_MW1_PISTA);
  chain->SetBranchStatus("Ex240Pu",true);
  chain->SetBranchAddress("Ex240Pu",&Ex240Pu);
  chain->SetBranchStatus("Ex240PuNoCorr",true);
  chain->SetBranchAddress("Ex240PuNoCorr",&Ex240PuNoCorr);
  
  chain->SetBranchStatus("Ex236U",true);
  chain->SetBranchAddress("Ex236U",&Ex236U);
  chain->SetBranchStatus("Ex238U",true);
  chain->SetBranchAddress("Ex238U",&Ex238U);
  chain->SetBranchStatus("Chio_DE",true);
  chain->SetBranchAddress("Chio_DE",&Chio_DE);
  chain->SetBranchStatus("Chio_E",true);
  chain->SetBranchAddress("Chio_E",&Chio_E);

  /////// Histograms

  TH2F* DE_E=new TH2F("DE_E","DE_E",1000,0,200,1000,0,60);
  TH1F* hEx240Pu=new TH1F("hEx240Pu","hEx240Pu",500,-5,20);
  TH1F* hEx240PuNoCorr=new TH1F("hEx240PuNoCorr","hEx240PuNoCorr",500,-5,20);
  TH1F* hEx240Pu_div=new TH1F("hEx240Pu_div","hEx240Pu_div",500,-5,20);
  TH1F* hEx240Pu_chio=new TH1F("hEx240Pu_chio","hEx240Pu_chio",500,-5,20);
  TH1F* hEx236U=new TH1F("hEx236U","hEx236U",500,-5,20);
  TH1F* hEx238U=new TH1F("hEx238U","hEx238U",500,-5,20);
  TH1F* hEx238U_div=new TH1F("hEx238U_div","hEx238U_div",500,-5,20);
  TH1F* hEx238U_chio=new TH1F("hEx238U_chio","hEx238U_chio",500,-5,20);

  TH2F* hChio_DE_E=new TH2F("hChio_DE_E","hChio_DE_E",1000,0,20000,1000,0,30000);
  TH2F* hChio_DE_E_corr=new TH2F("hChio_DE_E_corr","hChio_DE_E_corr",1000,0,20000,1000,0,30000);
  TH2F* hChio1_DE_E=new TH2F("hChio1_DE_E","hChio1_DE_E",1000,0,20000,1000,0,20000);
  TH2F* hChioDE_Ex240Pu=new TH2F("hChioDE_Ex240Pu","hChioDE_Ex240Pu",500,0,20,900,0,30000);
  TH2F* hChioDEcorr_Ex240Pu=new TH2F("hChioDEcorr_Ex240Pu","hChioDEcorr_Ex240Pu",500,0,20,1000,0,30000);
  TH2F* hChioDEcorr_Ex238U=new TH2F("hChioDEcorr_Ex238U","hChioDEcorr_Ex238U",500,0,20,1000,0,30000);
  TH2F* hChioDEcorr_Ex236U=new TH2F("hChioDEcorr_Ex236U","hChioDEcorr_Ex236U",500,0,20,1000,0,30000);

  /////// Filling of histograms

  int Nentries=chain->GetEntries();

  for(int e=0;e<Nentries;++e)
  {
    if(e%1000000==0) cout << e*100./Nentries << "% \r" << flush;

    chain->GetEntry(e);


    if(fTS_TMW==0){
      if(cutBe10->IsInside(Eres,DeltaEcorr)) {
        hEx240Pu_div->Fill(Ex240Pu);
      }
      else if(cutC12->IsInside(Eres,DeltaEcorr)) {
        hEx238U_div->Fill(Ex238U);
      }
    }
    
    else if(fTS_TMW>0 && fTAC_MW1_PISTA>0)
    {
      DE_E->Fill(Eres,DeltaEcorr);
      hChio_DE_E->Fill(Chio_E,Chio_DE);

      if(Chio_E>4000 && Chio_E<15000) 
	{
	  int Zid=FindClosestZ(Chio_E,Chio_DE); //IdentifyZ(Chio_E,Chio_DE);
	  //	  if(Zid==-1) {Zid=FindClosestZ(Chio_E,Chio_DE); }//cout << "Zid=" << Zid << " Chio_E=" << Chio_E << " Chio_DE=" << Chio_DE << endl;}
	  if(Zid==-1) continue;
	  double reference=Reference[Zid];
	  double delta=gSpline3[Zid]->Eval(Chio_E) - gSpline3[Zid]->Eval(reference);
	  hChio_DE_E_corr->Fill(Chio_E,Chio_DE-delta);
	}

      if(cutBe10->IsInside(Eres,DeltaEcorr)) 
      {
        hEx240Pu->Fill(Ex240Pu);
        if(Chio_E>5000 && Chio_E<10000) hChioDE_Ex240Pu->Fill(Ex240Pu,Chio_DE);
        if(Chio_E>1000 && Chio_E<15000){
	  int Zid=FindClosestZ(Chio_E,Chio_DE); //IdentifyZ(Chio_E,Chio_DE);
	  //	  if(Zid==-1) Zid=FindClosestZ(Chio_E,Chio_DE); 
	  if(Zid==-1) continue;
	  double reference=Reference[Zid];
	  double delta=gSpline3[Zid]->Eval(Chio_E) - gSpline3[Zid]->Eval(reference);
	  hChioDEcorr_Ex240Pu->Fill(Ex240Pu,Chio_DE-delta);
	  //          hChioDEcorr_Ex240Pu->Fill(Ex240Pu,Chio_DE*10000./fspline->Eval(Chio_E));
          hEx240Pu_chio->Fill(Ex240Pu);
          hEx240PuNoCorr->Fill(Ex240PuNoCorr);
        }
      }
      else if(cutC12->IsInside(Eres,DeltaEcorr)) 
	{
	  hEx238U->Fill(Ex238U);
	  if(Chio_E>1000 && Chio_E<15000){ 
            hChioDEcorr_Ex238U->Fill(Ex238U,Chio_DE*10000./fspline->Eval(Chio_E));
            hEx238U_chio->Fill(Ex238U);
          }
	}
      else if(cutC14->IsInside(Eres,DeltaEcorr)) 
	{
	  hEx236U->Fill(Ex236U);
          if(Chio_E>1000 && Chio_E<15000){
            int Zid=FindClosestZ(Chio_E,Chio_DE); //IdentifyZ(Chio_E,Chio_DE);
            //	  if(Zid==-1) Zid=FindClosestZ(Chio_E,Chio_DE); 
            if(Zid==-1) continue;
	    double reference=Reference[Zid];
            double delta=gSpline3[Zid]->Eval(Chio_E) - gSpline3[Zid]->Eval(10000);

            //hChioDEcorr_Ex236U->Fill(Ex236U,Chio_DE*10000./fspline->Eval(Chio_E));
            hChioDEcorr_Ex236U->Fill(Ex236U,Chio_DE-delta);
          }
        }

      if(cutChio1->IsInside(Chio_E,Chio_DE)) hChio1_DE_E->Fill(Chio_E,Chio_DE);

  }
}

TH1F* hProfile=(TH1F*)hChio1_DE_E->ProfileX();
// TSpline3();

TFile *OutFile;
if(Nrun>1) OutFile=new TFile(Form("hist/HistoPISTA_run%d-%d.root",RunNumber[0],RunNumber[Nrun-1]), "recreate");
else if(Nrun==1) OutFile=new TFile(Form("hist/HistoPISTA_run%d.root",RunNumber[0]), "recreate");

DE_E->Write();
hChio_DE_E->Write();
hChio1_DE_E->Write();
hChio_DE_E_corr->Write();

hEx240Pu->Write();
hEx240PuNoCorr->Write();
hEx240Pu_div->Write();
hEx240Pu_chio->Write();
hEx238U->Write();
hEx238U_div->Write();
hEx238U_chio->Write();
hEx236U->Write();

hChioDE_Ex240Pu->Write();
hChioDEcorr_Ex240Pu->Write();
hChioDEcorr_Ex238U->Write();
hChioDEcorr_Ex236U->Write();

OutFile->Close();
}


int IdentifyZ(double Chio_E, double Chio_DE)
{

  if(cut[0]->IsInside(Chio_E,Chio_DE)) return 0;
  else if(cut[1]->IsInside(Chio_E,Chio_DE)) return 1;
  else if(cut[2]->IsInside(Chio_E,Chio_DE)) return 2;
  else if(cut[3]->IsInside(Chio_E,Chio_DE)) return 3;
  else if(cut[4]->IsInside(Chio_E,Chio_DE)) return 4;
  else if(cut[5]->IsInside(Chio_E,Chio_DE)) return 5;
  else if(cut[6]->IsInside(Chio_E,Chio_DE)) return 6;
  else if(cut[7]->IsInside(Chio_E,Chio_DE)) return 7;
  else if(cut[8]->IsInside(Chio_E,Chio_DE)) return 8;
  else if(cut[9]->IsInside(Chio_E,Chio_DE)) return 9;
  else if(cut[10]->IsInside(Chio_E,Chio_DE)) return 10;
  else if(cut[11]->IsInside(Chio_E,Chio_DE)) return 11;
  else if(cut[12]->IsInside(Chio_E,Chio_DE)) return 12;
  else if(cut[13]->IsInside(Chio_E,Chio_DE)) return 13;
  else if(cut[14]->IsInside(Chio_E,Chio_DE)) return 14;
  else if(cut[15]->IsInside(Chio_E,Chio_DE)) return 15;
  else if(cut[16]->IsInside(Chio_E,Chio_DE)) return 16;
  else if(cut[17]->IsInside(Chio_E,Chio_DE)) return 17;
  else if(cut[18]->IsInside(Chio_E,Chio_DE)) return 18;
  else if(cut[19]->IsInside(Chio_E,Chio_DE)) return 19;
  else if(cut[20]->IsInside(Chio_E,Chio_DE)) return 20;
  else if(cut[21]->IsInside(Chio_E,Chio_DE)) return 21;
  else if(cut[22]->IsInside(Chio_E,Chio_DE)) return 22;
  else if(cut[23]->IsInside(Chio_E,Chio_DE)) return 23;
  else if(cut[24]->IsInside(Chio_E,Chio_DE)) return 24;
  else if(cut[25]->IsInside(Chio_E,Chio_DE)) return 25;
  else if(cut[26]->IsInside(Chio_E,Chio_DE)) return 26;
  else if(cut[27]->IsInside(Chio_E,Chio_DE)) return 27;
  else if(cut[28]->IsInside(Chio_E,Chio_DE)) return 28;
  else if(cut[29]->IsInside(Chio_E,Chio_DE)) return 29;
  else if(cut[30]->IsInside(Chio_E,Chio_DE)) return 30;
  else if(cut[31]->IsInside(Chio_E,Chio_DE)) return 31;
  else if(cut[32]->IsInside(Chio_E,Chio_DE)) return 32;
  else if(cut[33]->IsInside(Chio_E,Chio_DE)) return 33;
  else return -1;
}

int FindClosestZ(double Chio_E, double Chio_DE)
{
  double DEmin=100000;
  int Zmin=-1;

  for(int i=2;i<34;i++)
  {
    //cout << "i= " << i << endl;
    double test_DE=gSpline3[i]->Eval(Chio_E);
    //cout << test_DE << endl;
    if(fabs(test_DE-Chio_DE)<DEmin && test_DE>1000)
    {
      //cout << "DEmin=" << DEmin << " / test_DE-Chio_DE=" << fabs(test_DE-Chio_DE) << " / test_DE=" << test_DE << " /Chio_DE=" << Chio_DE << endl;
      DEmin=fabs(test_DE-Chio_DE);
      Zmin=i;
    }
  }
  return Zmin;
}
