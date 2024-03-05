TChain *chain;
double* parameter;
const int nb_parameter = 10;
int iteration = 0;
int m_section = 9;

int NumericalMinimization(const char* minName = "Minuit", const char* algoName = "");
double ConstantFactor(const double* parameter);

///////////////////////////////////////////////////
void Minimization(int section=9){
  m_section = section;
  chain = new TChain("tree"); 
  chain->Add("SelectTree.root");

  double buffer[nb_parameter] = {1,1,1,1,1,1,1,1,1,1};
  parameter = new double[nb_parameter];
  parameter = buffer;

  ConstantFactor(parameter);
  NumericalMinimization("Minuit","Migrad");
}



////////////////////////////////////////////////////
int NumericalMinimization(const char* minName, const char* algoName){

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(10000);
  min->SetTolerance(0.01);
  min->SetPrecision(0.01);
  min->SetPrintLevel(1);

  ROOT::Math::Functor f(&ConstantFactor,nb_parameter);

  min->SetFunction(f);

  for(int i=0; i<nb_parameter; i++){
    string seg_name = "seg"+to_string(i+1);
    min->SetLimitedVariable(i,seg_name,parameter[i],0.01,parameter[i]-0.5,parameter[i]+0.5);
  }
  //min->SetLimitedVariable(nb_parameter-1,"scalor",parameter[nb_parameter-1],0.01,0,2);

  min->Minimize();

  const double* xs = min->X();

  ofstream ofile;
  string filename = "Chio_E_sec" + to_string(m_section) + ".cal";
  ofile.open(filename.c_str());
  double init_par[nb_parameter] = {0.8686,0.7199,0.6233,0.4697,0.9787,0.9892,2.1038,1.9429,1.754,2.5};
  cout << "**********************" << endl;
  cout << "Minimum : " << endl;
  for(int i=0; i<nb_parameter; i++){ 
    cout << "par" << i+1 << " = " << xs[i] << " -> " << xs[i]*init_par[i] << endl;
    string parname;
    parname = "IC_SEC"+to_string(m_section)+"_SEG"+to_string(i+1)+"_ALIGN";
    ofile << parname << " " << xs[i]*init_par[i] << endl;
  }
  ofile.close();
  cout << "MinValue = " << min->MinValue() << endl;
  ConstantFactor(xs);

  return 0;

}

////////////////////////////////////////////////////
double ConstantFactor(const double* parameter){

  double distance=0;
  iteration++;

  double FF_Gamma;
  double FF_AoQ;
  double FF_Q;
  double FF_Etot;
  int FPMW_Section;


  TICPhysics* IC = new TICPhysics();

  chain->SetBranchStatus("FF_Gamma","true");
  chain->SetBranchAddress("FF_Gamma",&FF_Gamma);

  chain->SetBranchStatus("FF_AoQ","true");
  chain->SetBranchAddress("FF_AoQ",&FF_AoQ);

  chain->SetBranchStatus("FF_Etot","true");
  chain->SetBranchAddress("FF_Etot",&FF_Etot);

  chain->SetBranchStatus("FF_Q","true");
  chain->SetBranchAddress("FF_Q",&FF_Q);

  chain->SetBranchStatus("FPMW_Section","true");
  chain->SetBranchAddress("FPMW_Section",&FPMW_Section);


  chain->SetBranchStatus("fIC","true");
  chain->SetBranchAddress("IC",&IC);


  int nentries = chain->GetEntries();

  double Etot=0;
  double M1;
  double Q;
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);
    Etot = 0;
    M1 = 0;
    Q = 0;
    if(FPMW_Section==m_section && FF_AoQ>2){
      for(int j=0; j<nb_parameter-1; j++){
        Etot += parameter[j]*IC->fIC[j];

      }
      Etot = parameter[nb_parameter-1]*Etot;
      M1 = Etot/931.5016/(FF_Gamma-1);
      Q = M1/FF_AoQ;

      distance += pow((Q - 1540),2);
      /*if(cut1->IsInside(FF_AoQ,FF_Q)){
        if(Q>0) distance += pow((Q - 1410),2);
      }
      else if(cut2->IsInside(FF_AoQ,FF_Q)){
        if(Q>0) distance += pow((Q - 1605),2);
      }
      else if(cut3->IsInside(FF_AoQ,FF_Q)){
        if(Q>0) distance += pow((Q - 1800),2);
      }*/

      //cout << distance << " " << Etot << " " << M1 << " " << Q << endl;
    }
  }
  if(iteration%1==0){
    cout << "Number of iterations: " << iteration  << endl;
    cout << "distance= " << distance << endl;
  }

  return distance;
}
