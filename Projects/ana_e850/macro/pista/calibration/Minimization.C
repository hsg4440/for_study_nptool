TChain *chain;
double* parameter;
const int nb_parameter = 29;
int iteration = 0;
int m_telescope;

int NumericalMinimization(const char* minName = "Minuit", const char* algoName = "");
double ConstantFactor(const double* parameter);

NPL::Reaction* elastic;

TH2F* hEx;

///////////////////////////////////////////////////
void Minimization(int telescope=3){

  elastic = new NPL::Reaction("238U(12C,12C)238U@1417");
  hEx = new TH2F("hEx","hEx",7000,0,7000,500,-10,10);

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
  buffer[27] = 0;
  buffer[28] = 0;

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
    int indice = 3*i;
    string par0_name = "p0_E" + to_string(i+1);
    string par1_name = "p1_E" + to_string(i+1);
    string par2_name = "p2_E" + to_string(i+1);

    min->SetLimitedVariable(indice,par0_name,0,0.01,-3,3);
    min->SetLimitedVariable(indice+1,par1_name,0.0029,0.0001,0.97*0.00285,1.03*0.00285);
    min->SetLimitedVariable(indice+2,par2_name,0.0,0.0001,0,1e-8);
  }
  min->SetLimitedVariable(24,"Xt",0,0.01,-5,5);
  min->SetLimitedVariable(25,"Yt",0,0.01,-5,5);
  min->SetLimitedVariable(26,"Zt",0,0.01,-1,1);
  min->SetLimitedVariable(27,"ThetaX",0,0.01,-0.5,0.5);
  min->SetLimitedVariable(28,"ThetaY",0,0.01,-0.5,0.5);

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
  ofile2 <<  " XTARGET_OFFSET " << xs[24] << endl;
  ofile2 <<  " YTARGET_OFFSET " << xs[25] << endl;
  ofile2 <<  " ZTARGET_OFFSET " << xs[26] << endl;
  ofile2 <<  " BEAM_THETAX " << xs[27] << endl;
  ofile2 <<  " BEAM_THETAY " << xs[28] << endl;
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



  int nentries = 1e5;//chain->GetEntries();

  TH1F* hExmean = new TH1F("hExmean","hExmean",500,-10,10);
  TH2F* hEx_Theta = new TH2F("hEx_Theta","hEx_Theta",500,30,60,500,-10,10);
  double Elab = 0;
  double ThetaLab = 0;
  double Ecal = 0;
  double distance=100;
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    int indice = 3*(Telescope-1);
    Ecal = parameter[indice] + parameter[indice+1]*Eres + parameter[indice+2]*Eres*Eres;
    Elab = DeltaE + Ecal;

    if(Xcalc!=-1000 && Ycalc!=-1000 && Zcalc!=-1000){
      TVector3 PositionOnTarget = TVector3(parameter[24],parameter[25],parameter[26]);
      TVector3 HitPosition = TVector3(Xcalc,Ycalc,Zcalc);

      TVector3 BeamDirection = TVector3(0,0,1);
      BeamDirection.RotateX(parameter[27]*3.1415/180);
      BeamDirection.RotateY(parameter[28]*3.1415/180);
      TVector3 HitDirection = HitPosition - PositionOnTarget;
      ThetaLab = HitDirection.Angle(BeamDirection);

      double Ex = elastic->ReconstructRelativistic(Elab,ThetaLab);

      if(abs(Ex)<10){
        hEx->Fill(iteration,Ex);
        hExmean->Fill(Ex);
        hEx_Theta->Fill(ThetaLab*180./3.1415,Ex);
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
 
  distance = pow(Emean,2) + pow(sigma-0.3,2);
  if(iteration%1==0){
    cout << "Number of iterations: " << iteration  << endl;
    cout << "Emean= " << Emean << endl;
    cout << "distance= " << distance << endl;
  }

  hEx_Theta->Draw("colz");
  //delete hExmean;
  return distance;
}
