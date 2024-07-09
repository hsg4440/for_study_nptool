double GetAsymmetry(string, int Zb);
double GetRMS(string, int Zb);
double GetAsymmetryMoller(string, int Zb);
double GetAsymmetrySmooth(string, int Zb);


TH2I* h1;
TH2I* h2;
TH2I* h3;
TGraph* gStableNuclei;
TH2I* hStableNuclei;

vector<int> vZ;
vector<int> vA;
vector<int> vN;
vector<string> yield_file;

int NumberOfZ=0;
int NumberOfN=0;

////////////////////////////////////////////
void MapAsymmetry()
{
  gROOT->SetStyle("pierre_style");
  
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  
  NPL::Particle* npparticle = new NPL::Particle(1,1);
  gStableNuclei = npparticle->GetStableNuclei();
  gStableNuclei->SetMarkerStyle(21);
  gStableNuclei->SetFillColor(1);
  
  hStableNuclei = new TH2I("hstable","hstable",150,0,150,100,0,100);
  int Npoint = gStableNuclei->GetN();
  double gN, gZ;
  for(int i=0; i<Npoint; i++){
     gStableNuclei->GetPoint(i,gN,gZ);
     hStableNuclei->Fill((int)gN, (int)gZ);
  }
  hStableNuclei->SetMarkerStyle(1);
  hStableNuclei->SetFillColor(1);

  h2 = new TH2I("h2","h2",150,0,150,95,0,95);
  h2->GetXaxis()->SetTitle("Neutron number N");
  h2->GetYaxis()->SetTitle("Proton number Z");
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->SetMarkerStyle(1);
  h2->SetMarkerColor(0);
  
  h3 = new TH2I("h3","h3",150,0,150,100,0,100);
  h3->GetXaxis()->SetTitle("Neutrons, N");
  h3->GetYaxis()->SetTitle("Protons, Z");
  h3->GetZaxis()->SetTitle("Asymmetry");
  h3->GetXaxis()->CenterTitle();
  h3->GetYaxis()->CenterTitle();
  h3->GetZaxis()->CenterTitle();
  
  h1 = new TH2I("h1","h1",150,0,150,100,0,100);
  h1->GetXaxis()->SetTitle("N");
  h1->GetYaxis()->SetTitle("Z");
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();


  string filename = "list_all.dat";
  //string filename = "list_s455.dat";
  ifstream ifile;
  ifile.open(filename.c_str());

  ofstream ofile;
  string output_filename = "asymmetry_file.dat";
  ofile.open(output_filename.c_str());

  int Z=0;
  int Zprev=0;
  int A;
  int N;
  NPL::Particle* iso;
  while(ifile>>Z>>A){
    N = A-Z;
    //ifile >> Z >> A;
    iso = new NPL::Particle(Z,A);
    string isoname = iso->GetName();
    if(Z==92)
      isoname.resize(4);
    else
      isoname.resize(5);
    string yield_filename = "dat/yield_" + isoname + ".dat";
    yield_file.push_back(yield_filename);
    //cout << yield_filename << endl;
    vZ.push_back(Z);
    vA.push_back(A);
    vN.push_back(A-Z);
  
    double asy = GetAsymmetry(yield_filename, Z)*100;
    double asy_moller = GetAsymmetryMoller(yield_filename, Z)*100;
    double asy_smooth = GetAsymmetrySmooth(yield_filename, Z);
    //double asy = GetRMS(yield_filename, Z)/Z*100;
  
    ofile << Z << " " << A-Z << " " << 102 - asy_smooth*100 << endl;

    h1->Fill(N,Z,asy);
    //h2->Fill(N,Z,asy_moller);
    h2->Fill(N,Z,102-asy_smooth*100);
    //h3->Fill(N,Z,110-GetAsymmetryMoller(yield_filename, Z)*100);
    h3->Fill(N,Z,101-GetAsymmetrySmooth(yield_filename, Z)*100);
  }
 
  TExec* ex1 = new TExec("ex1","gStyle->SetPalette(kBird)");
  TExec* ex2 = new TExec("ex2","gStyle->SetPalette(kDarkBodyRadiator)");


  h2->GetXaxis()->SetRangeUser(90,150);
  h2->GetYaxis()->SetRangeUser(70,95);
  TCanvas* c1 = new TCanvas("c1","c1",1000,800);
  c1->cd(1);
  //hStableNuclei->Draw("a");
  h2->Draw("colz");
  ex2->Draw();
  hStableNuclei->Draw("colzsame");
  ex1->Draw();
  //h2->Draw("surf2");
  h2->Draw("colzsame");
  //hStableNuclei->Draw("colsame");
 
  TGraph* gX[200];
  TGraph* gY[200];
  for(int i=0; i<200; i++){
    gX[i] = new TGraph();
    gX[i]->SetPoint(0,i,0);
    gX[i]->SetPoint(1,i,150);

    gY[i] = new TGraph();
    gY[i]->SetPoint(0,0,i);
    gY[i]->SetPoint(1,150,i);
   
    gX[i]->SetLineStyle(1);
    gX[i]->SetLineWidth(1);
    gY[i]->SetLineStyle(1);
    gY[i]->SetLineWidth(1);
    
    gX[i]->Draw("lsame");
    gY[i]->Draw("lsame");
  }

  int magic[6] = {4,20,28,50,82,126};
  TGraph* gXmagic_down[6];
  TGraph* gXmagic_up[6];
  TGraph* gYmagic_down[6];
  TGraph* gYmagic_up[6];
  for(int i=0; i<6; i++){
    gXmagic_down[i] = new TGraph();
    gXmagic_down[i]->SetPoint(0,magic[i],0);
    gXmagic_down[i]->SetPoint(1,magic[i],150);

    gXmagic_up[i] = new TGraph();
    gXmagic_up[i]->SetPoint(0,magic[i]+1,0);
    gXmagic_up[i]->SetPoint(1,magic[i]+1,150);

    gYmagic_down[i] = new TGraph();
    gYmagic_down[i]->SetPoint(0,0,magic[i]);
    gYmagic_down[i]->SetPoint(1,150,magic[i]);
  
    gYmagic_up[i] = new TGraph();
    gYmagic_up[i]->SetPoint(0,0,magic[i]+1);
    gYmagic_up[i]->SetPoint(1,150,magic[i]+1);
  
    gXmagic_down[i]->SetLineStyle(1);
    gXmagic_down[i]->SetLineWidth(3);
    gXmagic_down[i]->SetLineColor(2);
    gYmagic_down[i]->SetLineStyle(1);
    gYmagic_down[i]->SetLineWidth(3);
    gYmagic_down[i]->SetLineColor(2);
 
    gXmagic_up[i]->SetLineStyle(1);
    gXmagic_up[i]->SetLineWidth(3);
    gXmagic_up[i]->SetLineColor(2);
    gYmagic_up[i]->SetLineStyle(1);
    gYmagic_up[i]->SetLineWidth(3);
    gYmagic_up[i]->SetLineColor(2);
  
    gXmagic_down[i]->Draw("lsame");
    gYmagic_down[i]->Draw("lsame");
    
    gXmagic_up[i]->Draw("lsame");
    gYmagic_up[i]->Draw("lsame");
  }
 

  h3->GetXaxis()->SetTitle("Neutrons, N");
  h3->GetYaxis()->SetTitle("Protons, Z");
  h3->GetZaxis()->SetTitle("Asymmetry");
  h3->GetXaxis()->CenterTitle();
  h3->GetYaxis()->CenterTitle();
  h3->GetZaxis()->CenterTitle();
  h3->GetXaxis()->SetTitleOffset(1.9);
  h3->GetYaxis()->SetTitleOffset(1.9);

  hStableNuclei->GetXaxis()->SetRangeUser(90,150);
  hStableNuclei->GetYaxis()->SetRangeUser(70,100);

  h3->SetFillColor(1);
  h3->SetLineColor(1);
  hStableNuclei->SetLineColor(0);
  //h3->Add(hStableNuclei);
 
  h3->GetXaxis()->SetRangeUser(90,150);
  h3->GetYaxis()->SetRangeUser(70,100);

  // X axis //
  h3->GetXaxis()->SetAxisColor(0);
  h3->GetXaxis()->SetTickSize(0);
  h3->GetXaxis()->SetLabelOffset(0.03);
  h3->GetXaxis()->SetLabelSize(0);
  h3->GetXaxis()->SetTitle("");
  
  // Y axis //
  h3->GetYaxis()->SetAxisColor(0);
  h3->GetYaxis()->SetTickSize(0);
  h3->GetYaxis()->SetLabelOffset(0.03);
  h3->GetYaxis()->SetLabelSize(0);
  h3->GetYaxis()->SetTitle("");
  
  // Z axis //
  h3->GetZaxis()->SetAxisColor(0);
  h3->GetZaxis()->SetTickSize(0);
  h3->GetZaxis()->SetLabelColor(0);
  h3->GetZaxis()->SetTitle("");

 TCanvas* c2 = new TCanvas("c2","c2",1600,1000);
  c2->cd();
  
  h3->Draw("0lego2fbbb");
  double x[2] = {85,130};
  double y[2] = {80.5,80.5};
  double z[2] = {0,0};

  TPolyLine3D* line = new TPolyLine3D(2,x,y,z);
  line->SetLineStyle(2);
  line->Draw();

  x[0] = 110;
  x[1] = 145;
  y[0] = 90.5;
  y[1] = 90.5;
  line = new TPolyLine3D(2,x,y,z);
  line->SetLineStyle(2);
  line->Draw();

  x[0] = 100.5;
  x[1] = 100.5;
  y[0] = 69;
  y[1] = 87;
  line = new TPolyLine3D(2,x,y,z);
  line->SetLineStyle(2);
  line->Draw();

  x[0] = 120.5;
  x[1] = 120.5;
  y[0] = 69;
  y[1] = 92;
  line = new TPolyLine3D(2,x,y,z);
  line->SetLineStyle(2);
  line->Draw();

  x[0] = 140.5;
  x[1] = 140.5;
  y[0] = 69;
  y[1] = 95;
  line = new TPolyLine3D(2,x,y,z);
  line->SetLineStyle(2);
  line->Draw();

  h3->Draw("0lego2fbbbsame");

  TExec* ex3 = new TExec("ex3","gStyle->SetPalette(kAvocado)");
  ex3->Draw();
 
  hStableNuclei->Draw("0lego2fbbbsame");

  ex1->Draw();
}

////////////////////////////////////////////
double GetAsymmetry(string filename, int Zb)
{
  double asymmetry;
  ifstream ifile;

  cout << "*** opening " << filename << endl;
  ifile.open(filename.c_str());

  vector<int> v_Z;
  vector<double> v_yield;
  vector<double> v_yield_err;

  int Z;
  double yield;
  double yield_err;
  while(!ifile.eof()){
    ifile >> Z >> yield >> yield_err;
    //cout << Z << " " << yield << " " << yield_err << endl;
    if(yield>0 && yield_err>0){
      v_Z.push_back(Z);
      v_yield.push_back(yield);
      v_yield_err.push_back(yield_err);
    }
  }
  ifile.close();
 
  unsigned int size = v_Z.size();
  double norm=0;
  double denum=0;
  for(unsigned int i=0; i<size; i++){
    norm += pow(v_yield[i],2)*pow(v_Z[i]-(double)Zb/2,2);
    denum += pow(v_yield[i],2);
  }

  asymmetry = sqrt(norm/denum);

  return asymmetry;
}

////////////////////////////////////////////
double GetAsymmetryMoller(string filename, int Zb)
{
  double asymmetry;
  ifstream ifile;

  ifile.open(filename.c_str());

  vector<int> v_Z;
  vector<double> v_yield;
  vector<double> v_yield_err;
  double yield_max=0;
  double yield_sym=0;
  int Zmax;

  int Z;
  double yield;
  double yield_err;
  while(!ifile.eof()){
    ifile >> Z >> yield >> yield_err;
    //cout << Z << " " << yield << " " << yield_err << endl;
    if(yield>0 && yield_err>0){
      v_Z.push_back(Z);
      v_yield.push_back(yield);
      v_yield_err.push_back(yield_err);
      if(yield>yield_max){
        yield_max = yield;
        Zmax = Z;
      }
      if(Z==int(Zb/2.)){
        yield_sym = yield;
      }
    }
  }
  ifile.close();

  asymmetry = yield_sym/yield_max;

  return asymmetry;
}

////////////////////////////////////////////
double GetAsymmetrySmooth(string filename, int Zb)
{
  double asymmetry;
  ifstream ifile;

  ifile.open(filename.c_str());

  vector<int> v_Z;
  vector<double> v_yield;
  vector<double> v_yield_err;
  double yield_max=0;
  double yield_sym=0;
  int Zmax;

  int Z;
  double yield;
  double yield_err;
  int i=0;
  int i_max=0;
  int i_sym=0;
  while(!ifile.eof()){
    ifile >> Z >> yield >> yield_err;
    //cout << Z << " " << yield << " " << yield_err << endl;
    if(yield>0 && yield_err>0){
      v_Z.push_back(Z);
      v_yield.push_back(yield);
      v_yield_err.push_back(yield_err);
      if(yield>yield_max){
        yield_max = yield;
        Zmax = Z;
        i_max = i;
      }
      if(Z==round(Zb/2)){
        yield_sym = yield;
        i_sym = i;
      }
      i++;
    }
  }
  ifile.close();

  yield_sym = v_yield[i_sym] + 0.5*(v_yield[i_sym-1] + v_yield[i_sym+1]);
  yield_max = v_yield[i_max] + 0.5*(v_yield[i_max-1] + v_yield[i_max+1]);


  asymmetry = yield_sym/yield_max;
  if(asymmetry>1)
    asymmetry=1.;

  cout << Zb << " " << asymmetry << endl;
  return asymmetry;
}



////////////////////////////////////////////
double GetRMS(string filename, int Zb)
{
  double rms;
  ifstream ifile;

  ifile.open(filename.c_str());

  TGraphErrors* gerr = new TGraphErrors();

  int Z;
  double yield;
  double yield_err;
  int i=0;
  while(!ifile.eof()){
    ifile >> Z >> yield >> yield_err;
    if(yield>0 && yield_err>0){
      gerr->SetPoint(i,Z,yield);
      gerr->SetPointError(i,0,yield_err);
      i++;
    }
  }
  ifile.close();
 
  rms = gerr->GetRMS();


  return rms;
}


