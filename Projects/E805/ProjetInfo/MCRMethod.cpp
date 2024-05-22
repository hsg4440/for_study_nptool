#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<sstream>


#include<string.h>
#include"MCRMethod.h"
#include"mini.h"
using namespace std;




MCRMethod::MCRMethod()
{
  //Init();
  // MetroHast();
  Nappe();
  // NappeBis();
  //MinimizationFunction();
}


MCRMethod::MCRMethod(double lambd, double ini, double final, int n1)
{
  tini= ini;
  tfinale=final;
  lambda=lambd;
  n=n1;
}

MCRMethod::~MCRMethod()
{
}

//void MCRMethod::MetroHast() ///////////
/*{
  srand(time(NULL));
  double a,b,c;
  int i1,i2;
  fstream fich;
  double min;
  int*li=(int*)malloc(n*sizeof(int));
  int*lican=(int*)malloc(n*sizeof(int));
  int*limin=(int*)malloc(n*sizeof(int));
  double delt;
  double t;
  dist obj(n,li);
  dist objcan(n,lican);
  dist objmin(n,limin);
  obj.minimisation::initialise(recuit::read());
  objcan.minimisation::initialise(recuit::read());
  objmin.minimisation::initialise(recuit::read());
  obj.set();
  objcan.set();
  objmin.set();
  for(int k=0;k<n;k++)
    li[k]=k;
  for(int k=0;k<n;k++)
    limin[k]=k;
  objmin.initialise();
  min=objmin.read();
  t=tini;
  fich.open("MetroHast3.res",ios::out);
  while (t>tfinale)
    {
      i1= rand() % n;
      i2= rand() % n;
      for(int k=0;k<n;k++)
	  lican[k]=li[k];
      lican[i1]=li[i2];
      lican[i2]=li[i1];
      obj.initialise();
      objcan.initialise();
      objmin.initialise();
      a=obj.read();
      b=objcan.read();
      c=objmin.read();
      delt= b-a;
      if (exp(-delt/t)>=(double(rand())/double(RAND_MAX))){
	
	for(int k=0;k<n;k++){
	  li[k]=lican[k];
	}
	fich << t << " " << b << endl;
      }
      if (c > b){
	for(int k=0;k<n;k++){
	  limin[k]=lican[k];
	}
	min = b;
      }
      t=t*lambda;
      //cout << x << "  " << min << endl;
    }
  reponse=min;
  cout << endl << "L'ordre de parcours des " << n << " villes dans la configuration optimale est: " << endl;
  for(int k=0;k<n;k++)
    cout << limin[k] << " ";
  cout << endl << endl;
  fich.close();
  
    
}*/

void MCRMethod::MetroHast() ///////////
{
  int itmax = TMath::Log(TemperatureFinale/TemperatureInitiale)/TMath::Log(Lambda);
  std::cout << "itmax " << itmax << std::endl;
  srand(time(NULL));
  clock_t start = clock(), current, end;
  while (Temperature>TemperatureFinale){
        current = clock();
        Double_t Frac = 1.0*it/itmax;
        Double_t Time = ((double) (current-start)/CLOCKS_PER_SEC);
        Double_t TimeLeft = Time*(1/Frac - 1.);

        std::cout << "\rIteration : " << it
             << "/" << itmax
             << " --- "
             << Form("%.2f",100.*it/itmax) <<" %"
             << " --- "
             <<Form("%.00f it/sec",it/Time)
            <<" --- "
           << " Time Left : "<< Form("%d min ",(int)TimeLeft/60)
           << Form("%01.00d sec",(int)TimeLeft%60)
           << std::flush;
    RandomStep();
    double currdist = MinimizationFunction();
    dist.push_back(currdist);
    if(it > 0)
      std::cout << "\n///////////// " << exp((dist[it-1] - dist[it])/Temperature) << "\n";
    if (it == 0 || exp((dist[it-1] - dist[it])/Temperature)>=(double(rand())/double(RAND_MAX))){
      for(unsigned int i = 0; i < 4; i++){
        CATSPosXY[i] = CATSPosXYCan[i];
      }
      if(dist[itmin] > dist[it]){
        itmin = it;
      }
    }
    else
    {
      dist[it] = dist[it-1];
    }
    pos[0].push_back(CATSPosXY[0]);
    pos[1].push_back(CATSPosXY[1]);
    pos[2].push_back(CATSPosXY[2]);
    pos[3].push_back(CATSPosXY[3]);
    std::cout << "Temperature : " << Temperature << std::endl;
    std::cout << "Dist : " << dist[it] << std::endl;
    std::cout << "Pos : " << pos[0][it] << " " << pos[1][it] << " " << pos[2][it] << " " << pos[3][it] << " " << std::endl;
    Temperature *= Lambda;
    it++;  
    Graph->SetPoint(it,it,dist[it]);
  }
    
    std::cout << "\n\nFIN DE L ALGO :\n";
    std::cout << "Temperature : " << Temperature << std::endl;
    std::cout << "Dist : " << dist[itmin] << std::endl;
    std::cout << "Pos : " << pos[0][itmin] << " " << pos[1][itmin] << " " << pos[2][itmin] << " " << pos[3][itmin] << " " << std::endl;
  //std::cout << "Min " << std::min_element(dist.begin(), dist.end());
  Graph->Draw();
}

void MCRMethod::RandomStep(){
  for(unsigned int i = 0; i < 4; i++){
    CATSPosXYCan[i] = CATSPosXY[i];
  }
  unsigned int i = rand() % 4;
  if(CATSPosXY[i] <= CATSPosXYLim[2*i])
    CATSPosXYCan[i] = CATSPosXY[i]+Step;
  else if(CATSPosXY[i] >= CATSPosXYLim[2*i +1])
    CATSPosXYCan[i] = CATSPosXY[i]-Step;
  else
    CATSPosXYCan[i] = CATSPosXY[i] + 2*((rand() % 2) -0.5)*Step;
}

void MCRMethod::Init(){
  Graph = new TGraph();
  it = 0;
  itmin = 0;
  TemperatureInitiale = 10;
  Temperature = TemperatureInitiale;
  TemperatureFinale = 0.1;
  Step = 0.1;
  runmask[0] = "r0367_mask1";
  runmask[1] = "r0368_mask1";
  for(unsigned int i = 0; i < 2; i++){
    Chain[i] = new TChain("PhysicsTree");
    Chain[i]->Add(path+runmask[i]+".root"); //CATS1
    if (!(Chain[i]->GetEntries())){
      printf("ERROR in the Chain !! \n");
      return;
    }
  TreeReader[i] = new TTreeReader(Chain[i]);
  CATSPhysics_[i]= new TTreeReaderValue<TCATSPhysics>(*TreeReader[i],"CATS");
  DistRatio[i] = -(CATSPosZ[1] - MASKPosZ[i])/(CATSPosZ[1] - CATSPosZ[0]);
  std::cout << "test " << DistRatio[i] << std::endl;
  }
} 

double MCRMethod::MinimizationFunction(){

    for(unsigned int i = 0; i < 2; i++){
      auto HistMask = new TH2F(Form("HistMask%i",i+1),Form("HistMask%i",i+1),500,-50,50,500,-50,50);
      auto SpecMask = new TSpectrum2(50,2);
      TreeReader[i]->Restart();
      while (TreeReader[i]->Next())
      {
      CATSPhysics = &**(CATSPhysics_[i]);
      if(CATSPhysics->GetCATSMult()== 2){
        double X1Eff, Y1Eff, X2Eff, Y2Eff;
        for(unsigned int i = 0; i < CATSPhysics->GetCATSMult(); i++){
          UShort_t Det = CATSPhysics->GetCATSDet(i);
          if(Det == 1){
          X1Eff = CATSPhysics->GetCATSPosX(i)+CATSPosXYCan[0];
          Y1Eff = CATSPhysics->GetCATSPosY(i)+CATSPosXYCan[1];
          }
          else if(Det == 2){
          X2Eff = CATSPhysics->GetCATSPosX(i)+CATSPosXYCan[2];
          Y2Eff = CATSPhysics->GetCATSPosY(i)+CATSPosXYCan[3];
          }
          // std::cout << "test " << X2Eff - (X2Eff - X1Eff)*DistRatio[i] << " " << Y2Eff - (Y2Eff - Y1Eff)*DistRatio[i] << std::endl;
        }
          HistMask->Fill(X2Eff + (X2Eff - X1Eff)*DistRatio[i],Y2Eff + (Y2Eff - Y1Eff)*DistRatio[i]);
      }

  }
    //new TCanvas();
   // HistMask->Draw();
    Int_t nfound = SpecMask->Search(HistMask,2,"col",0.15);
    Double_t* PeaksPosX = SpecMask->GetPositionX();
    Double_t* PeaksPosY = SpecMask->GetPositionY();
    Double_t FirstX;
    Double_t FirstY;
    if(nfound > 0){
      FirstX = PeaksPosX[0];
      FirstY = PeaksPosY[0];
      BotLeft[i][0] = FirstX;
      BotLeft[i][1] = FirstY;
      BotRight[i][0] = FirstX;
      BotRight[i][1] = FirstY;
      TopLeft[i][0] = FirstX;
      TopLeft[i][1] = FirstY;
  }
    else{
      std::cout << "ERROR, No peak was found" << std::endl;
      return 0;
  }

    for(unsigned int j = 0; j < nfound; j++){
      if(-BotLeft[i][0]-BotLeft[i][1] < -PeaksPosX[j]-PeaksPosY[j]){
        BotLeft[i][0] = PeaksPosX[j];
        BotLeft[i][1] = PeaksPosY[j];
      }
      else if(BotRight[i][0]-BotRight[i][1] < PeaksPosX[j]-PeaksPosY[j]){
        BotRight[i][0] = PeaksPosX[j];
        BotRight[i][1] = PeaksPosY[j];
      }
      else if(-TopLeft[i][0]+TopLeft[i][1] < -PeaksPosX[j]+PeaksPosY[j]){
        TopLeft[i][0] = PeaksPosX[j];
        TopLeft[i][1] = PeaksPosY[j];
      }
  }
    std::cout << "\nBotLeft " << BotLeft[i][0] << " " << BotLeft[i][1] << std::endl; 
    std::cout << "BotRight " << BotRight[i][0] << " " << BotRight[i][1] << std::endl; 
    std::cout << "TopLeft " << TopLeft[i][0] << " " << TopLeft[i][1] << std::endl; 
    HistMask->Delete();
    } 

  double distdiff = 0;
  for(unsigned int i = 0; i < 2; i++){
    distdiff+= TMath::Sqrt(pow(BotLeft[i][0] - BotLeftGeo[2*i],2) + pow(BotLeft[i][1] - BotLeftGeo[2*i+1],2));
    distdiff+= TMath::Sqrt(pow(TopLeft[i][0] - TopLeftGeo[2*i],2) + pow(TopLeft[i][1] - TopLeftGeo[2*i+1],2));
    distdiff+= TMath::Sqrt(pow(BotRight[i][0] - BotRightGeo[2*i],2) + pow(BotRight[i][1] - BotRightGeo[2*i+1],2));
  }

  std::cout << distdiff << std::endl;
  return distdiff;
}

void MCRMethod::Nappe(){
  int X1 = 40;
  int X2 = 40;
  double stepx1 = 0.1;
  double stepx2 = 0.1;
  string Side = "X";


  TH1Map = new std::map<int,std::map<string,TH1F*>>;
  double norm[2*X1+1][2*X2+1];
  for(unsigned int i = 0; i < 2; i++){
    for(int i1 = -X1; i1 <= X1; i1++){
      for(int i2 = -X2; i2 <= X2; i2++){
        (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)] = new TH1F(Form("h_m%i_x1%i_x2%i",i, i1, i2),Form("h_m%i_x1%i_x2%i",i, i1, i2),200,-10,10);
      }
    }
  }
  // std::cout << "test2" << std::endl;
  runmask[0] = "r314_mask1";
  runmask[1] = "r315_mask2";
  for(unsigned int i = 0; i < 2; i++){
    Chain[i] = new TChain("PhysicsTree");
    Chain[i]->Add(path+runmask[i]+".root"); //CATS1
    if (!(Chain[i]->GetEntries())){
      printf("ERROR in the Chain !! \n");
      return;
    }
  // std::cout << "test3" << std::endl;
  TreeReader[i] = new TTreeReader(Chain[i]);
  CATSPhysics_[i]= new TTreeReaderValue<TCATSPhysics>(*TreeReader[i],"CATS");

  }
  
  
  TString CName = "CUTm_1";
  CFile[0] = new TFile(CName+".root");
  CUT[0] = (TCutG*)CFile[0]->FindObjectAny(CName);
  std::cout << CUT[0] << std::endl;
  // CName = "CUTm_2";
  CName = "CUTm_2";
  CFile[1] = new TFile(CName+".root");
  CUT[1] = (TCutG*)CFile[1]->FindObjectAny(CName);
  std::cout << CUT[1] << std::endl;
  for(unsigned int i = 0; i < 2; i++){
    while (TreeReader[i]->Next()){
      CATSPhysics = &**(CATSPhysics_[i]);
      if(CATSPhysics->PositionX.size() == 2){
  // std::cout << "test5 " << CATSPhysics->GetCATSMult() << " " << CATSPhysics->DetNumber[0] << " " << (CATSPhysics->PositionX).size() << " " << CATSPhysics->PositionY[i] << " " << CUT[i]->IsInside(CATSPhysics->PositionX[i],CATSPhysics->PositionY[i]) << std::endl;
        if(CATSPhysics->DetNumber[0] == 1 && CUT[i]->IsInside(CATSPhysics->PositionX[i],CATSPhysics->PositionY[i])){
          for(int i1 = -X1; i1 <= X1; i1++){
            for(int i2 = -X2; i2 <= X2; i2++){
              // std::cout << "test " <<  i1 << " " << i2 << std::endl;
  // std::cout << "test4" << std::endl;
              double x1, x2;
              if(Side == "X"){
                x1 = CATSPhysics->PositionX[0] +i1*stepx1;
                x2 = CATSPhysics->PositionX[1] +i2*stepx2;
              }
              else if(Side == "Y"){
                x1 = CATSPhysics->PositionY[0] +i1*stepx1;
                x2 = CATSPhysics->PositionY[1] +i2*stepx2;
              }
              (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)]->Fill(ProjectOnCats(i,x1,x2));
              norm[X1+i1][X2+i2] = 0;

            }
          }
        }
      }
    }
  }
  
  (*TH1Map)[0][Form("h_m%i_x1%i_x2%i",0, 0, 0)]->Draw();
  (*TH1Map)[0][Form("h_m%i_x1%i_x2%i",0, 1, 1)]->Draw("same");
  (*TH1Map)[0][Form("h_m%i_x1%i_x2%i",0, 2, 2)]->Draw("same");
  auto c1 = new TCanvas;
  TH2F* HeatMap[2];
  auto fitfunc = new TF1("fitfunc","gausn",-10,10);
  auto SumHeatMap = new TH2F("Sumheatmap","Sumheatmap",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  HeatMap[0] = new TH2F("heatmap_0","heatmap_0",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  HeatMap[1] = new TH2F("heatmap_1","heatmap_1",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  TSpectrum* Spec = new TSpectrum(1,1.);
  for(unsigned int i = 0; i < 2; i++){
    for(int i1 = -X1; i1 <= X1; i1++){
      for(int i2 = -X2; i2 <= X2; i2++){
        Int_t nfound = Spec->Search((*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)],2,"",0.15);
        Double_t* PeaksPosX = Spec->GetPositionX();
        fitfunc->SetParameters(200,PeaksPosX[0],1);
        fitfunc->SetParLimits(1,-10,10);
        (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)]->Fit(fitfunc,"QL","",-10,10);
        // std::cout << i1 << " " << i2 << " " << PeaksPosX[0] << std::endl;
        double value;
        if(Side == "X")
          value = abs(fitfunc->GetParameter(1) - MaskPosX[i]);
        else if(Side == "Y")
          value = abs(fitfunc->GetParameter(1) - MaskPosY[i]);
        HeatMap[i]->SetBinContent(X1+i1+1,X2+i2+1,value);
        //HeatMap[i]->SetBinContent(i1+1,i2+1,value);
        std::cout << norm[X1+i1][X2+i2] << std::endl;
        norm[X1+i1][X2+i2] += value*value;
        if(i ==1){
          norm[X1+i1][X2+i2] = sqrt(norm[X1+i1][X2+i2]);
     
         }
      }
    }
  }
  for(int i1 = -X1; i1 <= X1; i1++){
    for(int i2 = -X2; i2 <= X2; i2++){
      SumHeatMap->SetBinContent(X1+i1+1,X2+i2+1,norm[X1+i1][X2+i2]);
  
    }
  }
  // SumHeatMap->Add(HeatMap[0],HeatMap[1],1,1);
  for(unsigned int i = 0; i < 2; i++){
  auto c = new TCanvas();
  HeatMap[i]->Draw("colz");
  }
  auto c = new TCanvas();
  SumHeatMap->Draw("colz");
}


float MCRMethod::ProjectOnCats(unsigned int i, double x1, double x2){
  double tmask = (CATSPosZ[i] - MASKPosZ[i])/(CATSPosZ[1] - CATSPosZ[0]); 
  // double tmaskt = (-CATSPosZ[0])/(CATSPosZ[1] - CATSPosZ[0]); 
  if(i == 0)
    return x1 -(x2-x1)*tmask;
  else if(i == 1)
    return x2 - (x2-x1)*tmask;
  // else if(i == 1)
    // return x1 + (x2-x1)*tmask;
  else
    exit(1);
};

void MCRMethod::NappeBis(){
  int X1 = 10;
  int X2 = 10;
  double stepx1 = 0.2;
  double stepx2 = 0.2;
  string Side = "X";


  TH1Map = new std::map<int,std::map<string,TH1F*>>;
  double norm[2*X1+1][2*X2+1];
  for(unsigned int i = 0; i < 2; i++){
    for(int i1 = -X1; i1 <= X1; i1++){
      for(int i2 = -X2; i2 <= X2; i2++){
          (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)] = new TH1F(Form("h_m%i_x1%i_x2%i",i, i1, i2),Form("h_m%i_x1%i_x2%i",i, i1, i2),300,-30,30);
      }
    }
  }
  runmask[0] = "NPA_360";// called runmask, but here runmask is the target run
  runmask[1] = "NPA_364";// this one is the normal run to analyze pd
  for(unsigned int i = 0; i < 2; i++){
    Chain[i] = new TChain("PhysicsTree");
    Chain[i]->Add(path+runmask[i]+".root"); //CATS1
    if (!(Chain[i]->GetEntries())){
      printf("ERROR in the Chain !! \n");
      return;
    }
  // std::cout << "test3" << std::endl;
  TreeReader[i] = new TTreeReader(Chain[i]);
  CATSPhysics_[i]= new TTreeReaderValue<TCATSPhysics>(*TreeReader[i],"CATS");

  }
  MUST2Physics_= new TTreeReaderValue<TMust2Physics>(*TreeReader[1],"MUST2");
  ZDDPhysics_= new TTreeReaderValue<TZDDPhysics>(*TreeReader[1],"ZDD");
  TACPhysics_= new TTreeReaderValue<TTACPhysics>(*TreeReader[1],"TAC");
  M2_CsI_E_p_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_CsI_E_p");
  M2_CsI_E_d_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_CsI_E_d");
  M2_CsI_E_t_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_CsI_E_t");
  M2_CsI_E_a_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_CsI_E_a");
  M2_ELab_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_ELab");
  M2_ThetaLab_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_ThetaLab");
  M2_ThetaCM_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_ThetaCM");
  M2_X_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_X");
  M2_Y_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_Y");
  M2_Z_ = new TTreeReaderValue<std::vector<double>>(*TreeReader[1],"M2_Z");
  
  TString CName = "CUT_360";
  CFile[0] = new TFile(CName+".root");
  CUT[0] = (TCutG*)CFile[0]->FindObjectAny(CName);
  std::cout << CUT[0] << std::endl;
  // CName = "CUTm_2";
  CName = "CUT_364";
  CFile[1] = new TFile(CName+".root");
  CUT[1] = (TCutG*)CFile[1]->FindObjectAny(CName);
  
  ULong64_t GetEntries = Chain[1]->GetEntries();
  std::cout << CUT[1] << std::endl;
  for(unsigned int i = 0; i < 2; i++){
  clock_t start = clock(), current, end;
    while (TreeReader[i]->Next()){
        Long64_t cEntry = TreeReader[i]->GetCurrentEntry();;
        if (cEntry%100000 == 0){
            current = clock();
            Double_t Frac = 1.0*cEntry/GetEntries;
            Double_t Time = ((double) (current-start)/CLOCKS_PER_SEC);
            Double_t TimeLeft = Time*(1/Frac - 1.);

            std::cout << "\rEntry : " << cEntry
                 << "/" << GetEntries
                 << " --- "
                 << Form("%.2f",100.*cEntry/GetEntries) <<" %"
                 << " --- "
                 <<Form("%.00f RunEvt/sec",cEntry/Time)
                <<" --- "
               << " Time Left : "<< Form("%d min ",(int)TimeLeft/60)
               << Form("%01.00d sec",(int)TimeLeft%60)
               << std::flush;
        }
      CATSPhysics = &**(CATSPhysics_[i]);
      if(CATSPhysics->PositionX.size() == 2 && CATSPhysics->PositionY.size()==2){
  // std::cout << "test5 " << CATSPhysics->GetCATSMult() << " " << CATSPhysics->DetNumber[0] << " " << (CATSPhysics->PositionX).size() << " " << CATSPhysics->PositionY[i] << " " << CUT[i]->IsInside(CATSPhysics->PositionX[i],CATSPhysics->PositionY[i]) << std::endl;
        if(condition(i)){
          for(int i1 = -X1; i1 <= X1; i1++){
            for(int i2 = -X2; i2 <= X2; i2++){
              // std::cout << "test " <<  i1 << " " << i2 << std::endl;
  // std::cout << "test4" << std::endl;
              double x1, x2;
              if(Side == "X"){
                x1 = CATSPhysics->PositionX[0] +i1*stepx1;
                x2 = CATSPhysics->PositionX[1] +i2*stepx2;
              }
              else if(Side == "Y"){
                x1 = CATSPhysics->PositionY[0] +i1*stepx1;
                x2 = CATSPhysics->PositionY[1] +i2*stepx2;
              }
              (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)]->Fill(ProjectOnTarget(i,x1,x2));
              norm[X1+i1][X2+i2] = -1000;

            }
          }
        }
      }
    }
  }
  (*TH1Map)[0][Form("h_m%i_x1%i_x2%i",0, 0, 0)]->Draw();
  new TCanvas;
  (*TH1Map)[1][Form("h_m%i_x1%i_x2%i",1, 1, 1)]->Draw("");
  (*TH1Map)[1][Form("h_m%i_x1%i_x2%i",1, 2, 2)]->Draw("same");

  auto c1 = new TCanvas;
  TH2F* HeatMap[2];
  TF1 *fitfunc[2];
  fitfunc[0] = new TF1("fitfunc","gausn",-10,10);
  // fitfunc[1] = new TF1("fitfunc","gausn(0) + gausn(3)",-10,10);
  auto SumHeatMap = new TH2F("Sumheatmap","Sumheatmap",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  HeatMap[0] = new TH2F("heatmap_0","heatmap_0",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  HeatMap[1] = new TH2F("heatmap_1","heatmap_1",2*X1+1, -X1*stepx1, (X1+1)*stepx1, 2*X2+1, -X2*stepx2, (X2+1)*stepx2);
  for(unsigned int i = 0; i < 2; i++){
    TSpectrum* Spec = new TSpectrum(i+1,1.);
    for(int i1 = -X1; i1 <= X1; i1++){
      for(int i2 = -X2; i2 <= X2; i2++){
        
        Int_t nfound = Spec->Search((*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)],2,"",0.15);
        double value;
        if(nfound == 1){
          Double_t* PeaksPosX = Spec->GetPositionX();
        
          fitfunc[0]->SetParameters(200,PeaksPosX[0],1);
          (*TH1Map)[i][Form("h_m%i_x1%i_x2%i",i, i1, i2)]->Fit(fitfunc[0],"QL","",-30,30);

          value = fitfunc[0]->GetParameter(1);
        }
        else
          value = 1e5;
        HeatMap[i]->SetBinContent(X1+i1+1,X2+i2+1,value);

        if(i== 0)
          norm[X1+i1][X2+i2] = value;
        else if(i== 1)
          norm[X1+i1][X2+i2] = abs(norm[X1+i1][X2+i2] - value);
             
      }
    }
  }
  for(int i1 = -X1; i1 <= X1; i1++){
    for(int i2 = -X2; i2 <= X2; i2++){
      SumHeatMap->SetBinContent(X1+i1+1,X2+i2+1,norm[X1+i1][X2+i2]);
  
    }
  }
  // SumHeatMap->Add(HeatMap[0],HeatMap[1],1,1);
  for(unsigned int i = 0; i < 2; i++){
  auto c = new TCanvas();
  HeatMap[i]->Draw("colz");
  }
  auto c = new TCanvas();
  SumHeatMap->Draw("colz");
}

bool MCRMethod::condition(unsigned int i){
    return CUT[i]->IsInside(CATSPhysics->PositionOnTargetX,CATSPhysics->PositionOnTargetY);
}

float MCRMethod::ProjectOnTarget(unsigned int i, double x1, double x2){
  double tmaskt = (-CATSPosZ[0])/(CATSPosZ[1] - CATSPosZ[0]); 
  return x1 + (x2-x1)*tmaskt;
};


void MCRMethod::Unallocate(){
  M2_CsI_E_p     = **M2_CsI_E_p_ ;
  M2_CsI_E_d     = **M2_CsI_E_d_ ;
  M2_CsI_E_t     = **M2_CsI_E_t_ ;
  M2_CsI_E_a     = **M2_CsI_E_a_ ;
  M2_ELab       = **M2_ELab_ ;
  M2_ThetaLab   = **M2_ThetaLab_ ;
  M2_ThetaCM     = **M2_ThetaCM_ ;
  M2_X           = **M2_X_ ;
  M2_Y           = **M2_Y_ ;
  M2_Z           = **M2_Z_ ;
  // M2_dE          = **M2_dE_ ;
  Must2Physics   = **MUST2Physics_ ;

}




double MCRMethod::read()
{
  return reponse;
}
