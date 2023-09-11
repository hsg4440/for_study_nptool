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
  Init();
  MetroHast();
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
  runmask[0] = "r0314_algo";
  runmask[1] = "r0315_algo";
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









double MCRMethod::read()
{
  return reponse;
}
