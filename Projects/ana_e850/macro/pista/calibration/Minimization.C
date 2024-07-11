TChain *chain;
double* parameter;
const int nb_parameter = 27;
int iteration = 0;
int m_telescope;

int NumericalMinimization(const char* minName = "Minuit", const char* algoName = "");
double ConstantFactor(const double* parameter);

NPL::Reaction* elastic;
double TargetThickness = 0.44*micrometer;
double BeamEnergy = 5.955*238;

NPL::EnergyLoss C12C;
NPL::EnergyLoss C12Si;
NPL::EnergyLoss C12Al;
NPL::EnergyLoss U238C;
TGraph* geloss_C12C;
TGraph* geloss_C12Al;

TH2F* hEx;
TH2F* hEx_Tel;
TH2F* hEx_Theta;
TH1D* hExmean;

///////////////////////////////////////////////////
void Minimization(int telescope=3){
 
  C12C  = EnergyLoss("../../../EnergyLossTable/C12_C.G4table","G4Table",100);
  C12Si = EnergyLoss("../../../EnergyLossTable/C12_C.G4table","G4Table",100);
  C12Al = EnergyLoss("../../../EnergyLossTable/C12_C.G4table","G4Table",100);
  U238C = EnergyLoss("../../../EnergyLossTable/U238_C.G4table","G4Table",100);
  geloss_C12C = new TGraph("../../../EnergyLossTable/C12_C.dat");
  geloss_C12Al = new TGraph("../../../EnergyLossTable/C12_Al.dat");
 
  elastic = new NPL::Reaction("238U(12C,12C)238U@1412");
  BeamEnergy = U238C.Slow(BeamEnergy,0.5*TargetThickness,0);
  elastic->SetBeamEnergy(BeamEnergy);

  hEx = new TH2F("hEx","hEx",7000,0,7000,500,-10,10);
  hEx_Tel = new TH2F("hEx_Tel","hEx_Tel",8,1,9,500,-10,10);
  hEx_Theta = new TH2F("hEx_Theta","hEx_Theta",500,30,60,500,-10,10);
 
  m_telescope = telescope;
  chain = new TChain("tree"); 
  chain->Add("SelectTree.root");

  double buffer[nb_parameter];
  for(int i=0; i<8; i++){
    buffer[3*i] = 0;
    buffer[3*i+1] = 0.0029;
    buffer[3*i+2] = 0.;
  }
  buffer[24] = 0;
  buffer[25] = 0;
  buffer[26] = 0;
  //buffer[27] = 0;
  //buffer[28] = 0;

  parameter = new double[nb_parameter];
  parameter = buffer;

  ConstantFactor(parameter);
  NumericalMinimization("Minuit","Migrad");

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  hEx->Draw("colz");

  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  hEx_Theta->Draw("colz");
  c2->cd(2);
  hEx_Tel->Draw("colz");
}



////////////////////////////////////////////////////
int NumericalMinimization(const char* minName, const char* algoName){

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(100000);
  min->SetTolerance(0.001);
  min->SetPrecision(0.0001);
  min->SetPrintLevel(1);

  ROOT::Math::Functor f(&ConstantFactor,nb_parameter);

  min->SetFunction(f);

  for(int i=0; i<8; i++){
    int indice = 3*i;
    string par0_name = "p0_E" + to_string(i+1);
    string par1_name = "p1_E" + to_string(i+1);
    string par2_name = "p2_E" + to_string(i+1);

    min->SetLimitedVariable(indice,par0_name,0,0.01,-2,2);
    min->SetLimitedVariable(indice+1,par1_name,0.0028,0.0001,0.97*0.0029,1.03*0.0029);
    min->SetLimitedVariable(indice+2,par2_name,0.0,0.0001,0,1e-8);
  }
  min->SetLimitedVariable(24,"Xt",0,0.01,-2,2);
  min->SetLimitedVariable(25,"Yt",0,0.01,-2,2);
  min->SetLimitedVariable(26,"Zt",0,0.01,-2,2);
  //min->SetLimitedVariable(27,"ThetaX",0,0.01,-0.5,0.5);
  //min->SetLimitedVariable(28,"ThetaY",0,0.01,-0.5,0.5);

  min->Minimize();

  const double* xs = min->X();

  ofstream ofile;
  string filename = "PISTA_BACK_E_min.cal";
  ofile.open(filename.c_str());

  for(int i=0; i<8; i++){
    string token = "PISTA_T" + to_string(i+1) + "_BACK_E";

    ofile << token << " " << xs[3*i] << " " << xs[3*i+1] << " " << xs[3*i+2] << endl;
  }
  ofile.close();
 
  ofstream ofile2;
  string filename2 = "AnalysisConfig_min.dat";
  ofile2.open(filename2.c_str());
  ofile2 <<  "AnalysisConfig" << endl;
  ofile2 <<  " VAMOS_ANGLE 20 " << endl;
  ofile2 <<  " BRHO_REF 1.1 " << endl;
  ofile2 <<  " XTARGET_OFFSET " << xs[24] << endl;
  ofile2 <<  " YTARGET_OFFSET " << xs[25] << endl;
  ofile2 <<  " ZTARGET_OFFSET " << xs[26] << endl;
  //ofile2 <<  " BEAM_THETAX " << xs[27] << endl;
  //ofile2 <<  " BEAM_THETAY " << xs[28] << endl;
  ofile2.close();

  cout << "MinValue = " << min->MinValue() << endl;
  ConstantFactor(xs);

  return 0;

}

////////////////////////////////////////////////////
double ConstantFactor(const double* parameter){

  iteration++;

  double DeltaE;
  double Eres;
  double Xcalc;
  double Ycalc;
  double Zcalc;
  int Telescope;
  double XTarget;
  double YTarget;


  chain->SetBranchStatus("DeltaE","true");
  chain->SetBranchAddress("DeltaE",&DeltaE);

  chain->SetBranchStatus("Eres","true");
  chain->SetBranchAddress("Eres",&Eres);

  chain->SetBranchStatus("Xcalc","true");
  chain->SetBranchAddress("Xcalc",&Xcalc);

  chain->SetBranchStatus("Ycalc","true");
  chain->SetBranchAddress("Ycalc",&Ycalc);

  chain->SetBranchStatus("Zcalc","true");
  chain->SetBranchAddress("Zcalc",&Zcalc);

  chain->SetBranchStatus("XTarget","true");
  chain->SetBranchAddress("XTarget",&XTarget);

  chain->SetBranchStatus("YTarget","true");
  chain->SetBranchAddress("YTarget",&YTarget);

  chain->SetBranchStatus("Telescope","true");
  chain->SetBranchAddress("Telescope",&Telescope);



  int nentries = 1e5;//chain->GetEntries();

  hEx_Tel->Reset();
  hEx_Theta->Reset();
  double Elab = 0;
  double ThetaLab = 0;
  double Ecal = 0;
  double distance=100;
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    int indice = 3*(Telescope-1);
    Ecal = parameter[indice] + parameter[indice+1]*Eres + parameter[indice+2]*Eres*Eres;
   
    //Elab = DeltaE + Ecal;
    Elab = geloss_C12Al->Eval(Ecal);
    Elab = geloss_C12Al->Eval(Elab);
    Elab += DeltaE;
    Elab = geloss_C12Al->Eval(Elab);
    Elab = geloss_C12C->Eval(Elab);


    if(Xcalc!=-1000 && Ycalc!=-1000 && Zcalc!=-1000){
      TVector3 PositionOnTarget = TVector3(XTarget+parameter[24],YTarget+parameter[25],parameter[26]);
      TVector3 HitPosition = TVector3(Xcalc,Ycalc,Zcalc);

      TVector3 BeamDirection = TVector3(0,0,1);
      //BeamDirection.RotateX(parameter[27]*3.1415/180);
      //BeamDirection.RotateY(parameter[28]*3.1415/180);
      TVector3 HitDirection = HitPosition - PositionOnTarget;
      ThetaLab = HitDirection.Angle(BeamDirection);

      double Ex = elastic->ReconstructRelativistic(Elab,ThetaLab);
      if(abs(Ex)<5){
        hEx->Fill(iteration,Ex);
        hEx_Theta->Fill(ThetaLab*180./3.1415,Ex);
        hEx_Tel->Fill(Telescope,Ex);
      }

    }
  }

  double Emean[8];
  double sigma[8];
  double AvEmean = 0;
  for(int i=0; i<8; i++){
    hExmean = hEx_Tel->ProjectionY(Form("hExmean_%d",i+1),i+1,i+1);
    Emean[i] = hExmean->GetMean();
    double Emin = Emean[i]-3;
    double Emax = Emean[i]+3;
    hExmean->GetXaxis()->SetRangeUser(Emin,Emax);

    Emean[i] = hExmean->GetMean();
    sigma[i] = hExmean->GetRMS();

    distance += pow(Emean[i],2) + pow(sigma[i]-0.6, 2);
    AvEmean += Emean[i]; 
  }

  AvEmean/=8;
 
  if(iteration%1==0){
    cout << "Number of iterations: " << iteration  << endl;
    for(int k=0; k<nb_parameter; k++){
      cout << parameter[k] << endl;
    }
    cout << "Emean= " << AvEmean << endl;
    cout << "distance= " << distance << endl;
  }

  //delete hExmean;
  return distance;
}
