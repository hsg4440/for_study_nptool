TChain *chain;
double* parameter;
const int nb_parameter = 21;
int iteration = 0;
int m_telescope;

int NumericalMinimization(const char* minName = "Minuit", const char* algoName = "");
double ConstantFactor(const double* parameter);

NPL::Reaction* elastic;

TH2F* hEx;

///////////////////////////////////////////////////
void Minimization(int telescope=3){

  elastic = new NPL::Reaction("238U(12C,12C)238U@1417");
  hEx = new TH2F("hEx","hEx",4000,0,4000,500,-10,10);

  m_telescope = telescope;
  chain = new TChain("tree"); 
  chain->Add("SelectTree.root");

  double buffer[nb_parameter];
  for(int i=0; i<8; i++){
    buffer[2*i] = 0;
    buffer[2*i+1] = 0.0029;
  }
  buffer[16] = 0;
  buffer[17] = 0;
  buffer[18] = 0;
  buffer[19] = 0;
  buffer[20] = 0;

  parameter = new double[nb_parameter];
  parameter = buffer;

  ConstantFactor(parameter);
  NumericalMinimization("Minuit","Migrad");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  hEx->Draw("colz");
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
    int indice = 2*i;
    string par0_name = "p0_E" + to_string(i+1);
    string par1_name = "p1_E" + to_string(i+1);

    min->SetLimitedVariable(indice,par0_name,0,0.01,-5,5);
    min->SetLimitedVariable(indice+1,par1_name,0.0029,0.0001,0.95*0.0029,1.05*0.0029);
  }
  min->SetLimitedVariable(16,"Xt",0,0.01,-5,5);
  min->SetLimitedVariable(17,"Yt",0,0.01,-5,5);
  min->SetLimitedVariable(18,"Zt",0,0.01,-1,1);
  min->SetLimitedVariable(19,"ThetaX",0,0.01,-1,1);
  min->SetLimitedVariable(20,"ThetaY",0,0.01,-1,1);

  min->Minimize();

  const double* xs = min->X();

  ofstream ofile;
  string filename = "PISTA_BACK_E_min.cal";
  ofile.open(filename.c_str());
  cout << "**********************" << endl;
  cout << "Minimum : " << endl;
  cout <<  "p0= " << xs[0] << endl;
  cout <<  "p1= " << xs[1] << endl;
  cout <<  "Xt= " << xs[16] << endl;
  cout <<  "Yt= " << xs[17] << endl;
  cout <<  "Zt= " << xs[18] << endl;
  cout <<  "ThetaX= " << xs[19] << endl;
  cout <<  "ThetaY= " << xs[20] << endl;

  for(int i=0; i<8; i++){
    string token = "PISTA_T" + to_string(i+1) + "_BACK_E";

    ofile << token << " " << xs[2*i] << " " << xs[2*i+1] << endl;
  }
  ofile.close();
 
  ofstream ofile2;
  string filename2 = "AnalysisConfig_min.dat";
  ofile2.open(filename2.c_str());
  ofile2 <<  "AnalysisConfig" << endl;
  ofile2 <<  " XTARGET_OFFSET " << xs[16] << endl;
  ofile2 <<  " YTARGET_OFFSET " << xs[17] << endl;
  ofile2 <<  " ZTARGET_OFFSET " << xs[18] << endl;
  ofile2 <<  " BEAM_THETAX " << xs[19] << endl;
  ofile2 <<  " BEAM_THETAY " << xs[20] << endl;
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

  chain->SetBranchStatus("Telescope","true");
  chain->SetBranchAddress("Telescope",&Telescope);



  int nentries = 1e6;//chain->GetEntries();

  TH1F* hExmean = new TH1F("hExmean","hExmean",500,-10,10);
  double Elab = 0;
  double ThetaLab = 0;
  double Ecal = 0;
  double distance=100;
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    int indice = 2*(Telescope-1);
    Ecal = parameter[indice] + parameter[indice+1]*Eres;
    Elab = DeltaE + Ecal;

    if(Xcalc!=-1000 && Ycalc!=-1000 && Zcalc!=-1000 && Eres>25e3){
      TVector3 PositionOnTarget = TVector3(parameter[16],parameter[17],parameter[18]);
      TVector3 HitPosition = TVector3(Xcalc,Ycalc,Zcalc);

      TVector3 BeamDirection = TVector3(0,0,1);
      BeamDirection.RotateX(parameter[19]*3.1415/180);
      BeamDirection.RotateY(parameter[20]*3.1415/180);
      TVector3 HitDirection = HitPosition - PositionOnTarget;
      ThetaLab = HitDirection.Angle(BeamDirection);

      double Ex = elastic->ReconstructRelativistic(Elab,ThetaLab);

      if(abs(Ex)<10){
        hEx->Fill(iteration,Ex);
        hExmean->Fill(Ex);
      }

    }
  }

  double Emean = hExmean->GetMean();
  double Emin = Emean - 3;
  double Emax = Emean + 3;
  hExmean->GetXaxis()->SetRangeUser(Emin,Emax);
  TF1* f1 = new TF1("gaus","gaus",Emin, Emax);
  Emean = hExmean->GetMean();
  double sigma = hExmean->GetRMS();
 
  distance = pow(Emean,2) + pow(sigma-0.8,2);
  if(iteration%1==0){
    cout << "Number of iterations: " << iteration  << endl;
    cout << "Emean= " << Emean << endl;
    cout << "distance= " << distance << endl;
  }

  hExmean->Draw();
  //delete hExmean;
  return distance;
}
