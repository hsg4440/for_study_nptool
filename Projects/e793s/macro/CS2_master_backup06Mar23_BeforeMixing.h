/* Predefine functions */
void CS();
void CS_Diagnosis();
vector<vector<double>> GetExpDiffCross(double Energy);
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize);
TH1F* PullThetaCMHist(int i, double minTheta, double gatesize);
TH1F* PullPhaseSpaceHist(int i, double minTheta, double gatesize);
void Scale(TGraph* g , TGraphErrors* ex);
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j);
double ToMininize(const double* parameter);
TGraph* FindNormalisation(TGraph* theory, TGraphErrors* experiment);
TList* peakFitList = new TList();

/* Global variables */
vector<Double_t> anglecentres, anglewidth;
TGraph* currentThry;
TGraphErrors* staticExp;
int indexE;
double globalS, globalSerr;
int numAngleBins;
int numAngleBinsBase;
int numPeaks_CS;
double widthAngleBins;
double firstAngle;
double CSangleL, CSangleH;
vector<double> means_CS;
bool GammaGate = false;
double globGammaGate;

/* Output volume toggle */
bool loud = 0;

/* Scale method toggle */
bool scaleTogether = 1;

/* Strings for image */
string orbitalname;
string orbital;

/* Strings for SolidAngle input file */
string statename;
string inputdate;


////////////////////////////////////////////////////////////////////////////////
void canclone(TCanvas* major, int padNum, string name){
  string minorName = "./CS2_Figures/" + name + ".root";
  TFile* minorFile = new TFile(minorName.c_str(),"READ");
  TCanvas* minor = (TCanvas*) minorFile->FindObjectAny("c_peakFits");
 
  major->cd(padNum);
  TPad *pad = (TPad*)gPad;
  minor->cd();
  TObject *obj, *clone;
  pad->Range(minor->GetX1(),minor->GetY1(),minor->GetX2(),minor->GetY2());
  pad->SetTickx(minor->GetTickx());
  pad->SetTicky(minor->GetTicky());
  pad->SetGridx(minor->GetGridx());
  pad->SetGridy(minor->GetGridy());
  pad->SetLogx(minor->GetLogx());
  pad->SetLogy(minor->GetLogy());
  pad->SetLogz(minor->GetLogz());
  pad->SetBorderSize(minor->GetBorderSize());
  pad->SetBorderMode(minor->GetBorderMode());
  minor->TAttLine::Copy((TAttLine&)*pad);
  minor->TAttFill::Copy((TAttFill&)*pad);
  minor->TAttPad::Copy((TAttPad&)*pad);

  TIter next(minor->GetListOfPrimitives());
  gROOT->SetSelectedPad(pad);
  while ((obj=next())) {
     clone = obj->Clone();
     obj->Draw("SAME");
     pad->GetListOfPrimitives()->Add(clone,obj->GetDrawOption());
  }
  pad->Modified();
  pad->Update();
  major->cd(padNum);
  pad->Draw();
}

////////////////////////////////////////////////////////////////////////////////
void CS(double Energy, double Spin, double spdf, double angmom){
  /* Clean global variables, in case of multiple runs */
  indexE = 100;
  anglecentres.clear();
  anglewidth.clear();
  globalS=0.;
  globalSerr=0.;
  peakFitList->Clear();
  numAngleBins=numAngleBinsBase;

  /* Assign local variables */
  // p3/2 -> spdf = 1, angmom = 1.5
  // J0 is incident spin, which is 47K g.s. therefore J0 = 1/2
  double J0 = 0.5;
  double ElasticNorm = 0.000234, ElasticNormErr = 0.000; //18Oct22
  inputdate = "18Oct22";

  string backupFileName = "SolidAngle_HistFiles/";
  string saFileName = "SolidAngle_HistFiles/SAHF_";
  saFileName.append(inputdate);

  if(reactionName=="47K(d,p)"){
    numAngleBins=20;
    numAngleBinsBase=20;
    numPeaks_CS = numPeaks;
    widthAngleBins=2.5;
    firstAngle=105.0;
    backupFileName.append("SolidAngle_HistFile_10Aug22_TrueStripRemoval.root");
    saFileName.append("_47Kdp_");
    CSangleL = 100.;
    CSangleH = 160.;
    means_CS=means;
  } else if(reactionName=="47K(d,t)"){
    numAngleBins=18;
    numAngleBinsBase=18;
    numPeaks_CS = numPeaks_dt;
    widthAngleBins=1.0;
    firstAngle=2.0;
    backupFileName.append("SolidAngle_HistFile_18Oct22_47Kdt.root");
    saFileName.append("_47Kdt_");
    CSangleL = 00.;
    CSangleH = 30.;
    means_CS=means_dt;
  }

  // Extract spin orbit name
  orbitalname.clear();
  orbital.clear();
  
  switch((int)spdf){
    case 0:
      orbitalname.append("s_{");
      orbital.append("s");
      break;
    case 1:
      orbitalname.append("p_{");
      orbital.append("p");
      break;
    case 2:
      orbitalname.append("d_{");
      orbital.append("d");
      break;
    case 3:
      orbitalname.append("f_{");
      orbital.append("f");
      break;
    default:
      cout << RED
	   <<"ERROR: SPDF SELECTION OUT OF RANGE" 
	   << endl;
  }
  orbitalname.append(to_string((int)(2*angmom)));
  orbitalname.append("/2}");
  orbital.append(to_string((int)(2*angmom)));
  orbital.append("2");

  // Number of nodes
  //CHECK THE NODES HERE ARE RIGHT!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double nodes;
  if(spdf==1 || spdf==0){
    nodes=1;
  }
  else if(spdf==3 || spdf==2){
    nodes=0;
  }
  else{
    cout << " INPUT NODES::" << endl;
    cin >> nodes;
  }

  /* Retrieve array index of the entered peak energy */
  /* numpeaks and Energy[] defined globally in KnownPeakFitter.h */
  bool found = 0;
  for(int i=0;i<numPeaks_CS;i++){
    if(abs(Energy-means_CS.at(i))<0.01){
      cout << "========================================================" << endl;
      cout << "Identified as state #" << i << ", E = " << means_CS.at(i) << endl;
      indexE = i;
      found = 1;
      stringstream ss;
      ss << setfill('0') << setw(4) << (int)(means_CS.at(i)*1000.);
      statename = ss.str();
    }
  }
  if(!found){
    cout << "========================================================" << endl;
    cout << "NO STATE AT THAT ENERGY INDENTIFIED!! CHECK KNOWN PEAKS!!" << endl;
    return;
  }
  saFileName.append(statename);
  saFileName.append(".root");

  /* Solid Angle (from simulation) */
  TFile* file;
  ifstream infile(saFileName);
  ifstream backup(backupFileName);
  if(infile.good()){
    cout << GREEN << "Opening file " << saFileName << RESET << endl;
      file = TFile::Open(saFileName.c_str());
  } else if (backup.good()){
    cout << BOLDRED << "FAILED TO OPEN " << saFileName << endl;
    cout << "Open standard backup file " << backupFileName << endl;
    cout << RESET;
      file = TFile::Open(backupFileName.c_str());
  } else {
    cout << BOLDRED << "FAILED TO OPEN MAIN OR BACKUP SOLID ANGLE FILE" << endl;
    cout << RED << "Check SolidAngle file exists..." << RESET << endl;
  }

  string objName;
  if(reactionName=="47K(d,p)"){
    objName.append("SolidAngle_Lab_MG");
  }
  else if(reactionName=="47K(d,t)"){
    objName.append("SolidAngle_CM_MM");
  }

  TH1F* SolidAngle = (TH1F*) file->FindObjectAny(objName.c_str());
  TCanvas* c_SolidAngle = new TCanvas("c_SolidAngle","c_SolidAngle",1000,1000);
  SolidAngle->Draw("HIST");
  SolidAngle->GetXaxis()->SetRangeUser(CSangleL,CSangleH);
  /* (canvas deleted after Area/SA calculation) */
 
  /* Area of experimental peaks */
  TCanvas* c_PeakArea = new TCanvas("c_PeakArea","c_PeakArea",1000,1000);
  vector<vector<double>> areaArray = GetExpDiffCross(means_CS.at(indexE));
  delete c_PeakArea;

  // Array: peakenergy, peakarea, areaerror, anglemin, anglemax
  if(loud){
    for(int i=0; i<areaArray.size();i++){
      cout << i << " " 
    	   << areaArray[i][0] << " " 
	   << areaArray[i][1] << " "
	   << areaArray[i][2] << " "
	   << areaArray[i][3] << " "
	   << areaArray[i][4] << endl;
    }
  }

  /* AoSA = Area/Solid Angle [#/msr] */
  /* dSdO = Experimental Diff. Cross Sect. (Area/Solid Angle *Norm) [mb/msr] */
  vector<Double_t> AoSA, AoSAerr;
  vector<Double_t> expDCS, expDCSerr;
  for(int i=0; i<areaArray.size();i++){
    int binmin = SolidAngle->FindBin(areaArray[i][3]+0.0001);
    int binmax = SolidAngle->FindBin(areaArray[i][4]-0.0001);

    anglecentres.push_back(((areaArray[i][4]-areaArray[i][3])/2.)+areaArray[i][3]);
    anglewidth.push_back(areaArray[i][4]-areaArray[i][3]);

    double tempsum=0, tempsumerr=0;
    for(int x=binmin;x<binmax+1;x++){
      tempsum += SolidAngle->GetBinContent(x);
      tempsumerr += SolidAngle->GetBinError(x);
      if(loud){cout << x << endl;}
    }
    if(loud){cout << "TEST CHECK " << tempsum << " +- " << tempsumerr << endl;}

    double SAerr;
    double SA = SolidAngle->IntegralAndError(binmin,binmax,SAerr);
    SA = SA*1000.;       //sr->msr
    SAerr = SAerr*1000.; //sr->msr

    /* Area over Solid angle ONLY */
    AoSA.push_back(areaArray[i][1]/SA);
    AoSAerr.push_back((areaArray[i][1]/SA) 
		    * sqrt( pow(areaArray[i][2]/areaArray[i][1],2) 
			    + pow(SAerr/SA,2)));

    /* Differential Cross Section */
    /* NOTE: DON'T INCLUDE NORM ERROR IN ERR PROPAGATION AS IT IS SYSTEMATIC! */
    double tempvalue = (areaArray[i][1]/SA)*ElasticNorm; 
    expDCS.push_back(tempvalue);
    double temperror = tempvalue
		     * sqrt( pow(areaArray[i][2]/areaArray[i][1],2)
		           + pow(SAerr/SA,2)
			   );
    expDCSerr.push_back(temperror);

    if(loud){
      cout << "Angle " << areaArray[i][3] << " to " << areaArray[i][4] 
	   << ", centre " << anglecentres[i]
	   << ": Area = " << areaArray[i][1] << " +- " << areaArray[i][2] << " cnts" 
	   << "  SA = " << SA << " +- " << SAerr << " msr" 
           << "  Area/SA = " << AoSA[i] << " +- " << AoSAerr[i] << " counts/msr"
	   << setprecision(5)
           << "  Norm = " << ElasticNorm << " +- " << ElasticNormErr
	   << " (don't include norm err in propagation)"
           << "  dSdO = " << tempvalue << " +- " << temperror  
	   << setprecision(3)
	   << endl;
    }

  }
  //delete c_SolidAngle;
  
  /* Graph of Area/Solid Angle*/
  TCanvas* c_AoSA = new TCanvas("c_AoSA","c_AoSA",1000,1000);
  c_AoSA->SetLogy();
  TGraphErrors* gAoSA = new TGraphErrors(
		  anglecentres.size(), //n
		  &(anglecentres[0]), &(AoSA[0]), //x, y
		  0, &(AoSAerr[0]) );  //errX, errY 
  gAoSA->SetTitle("Area/SolidAngle");
  gAoSA->GetXaxis()->SetTitle("ThetaLab [deg]");
  gAoSA->GetYaxis()->SetTitle("Counts/#Omega [counts/msr]");
  gAoSA->Draw();

  /* Graph of experimental diff. cross sect (dSigma/dOmega) */
  TCanvas* c_dSdO = new TCanvas("c_dSdO","c_dSdO",1000,1000);
  c_dSdO->SetLogy();
  TGraphErrors* gdSdO = new TGraphErrors(
		  anglecentres.size(),
		  &(anglecentres[0]), &(expDCS[0]),
		  0, &(expDCSerr[0]) );
  gdSdO->SetTitle("Differential Cross Section");
  if(reactionName=="47K(d,p)"){
    gdSdO->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  } else if (reactionName=="47K(d,t)"){
    gdSdO->GetXaxis()->SetTitle("#theta_{CM} [deg]");
  }
  gdSdO->GetXaxis()->SetTitleOffset(1.2);
  gdSdO->GetXaxis()->SetTitleSize(0.04);
  gdSdO->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/msr]");
  gdSdO->GetYaxis()->SetTitleOffset(1.2);
  gdSdO->GetYaxis()->SetTitleSize(0.04);
  gdSdO->Draw();
  c_dSdO->Update();

  /* TWOFNR diff. cross section, in mb/msr */ 
  TCanvas* c_TWOFNR = new TCanvas("c_TWOFNR","c_TWOFNR",1000,1000);
  c_TWOFNR->SetLogy();
  TGraph* TheoryDiffCross = TWOFNR(means_CS.at(indexE), J0, Spin, nodes, spdf, angmom); 
  TheoryDiffCross->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/msr]"); //msr set in func above
  TheoryDiffCross->GetXaxis()->SetTitle("ThetaLab [deg]");
  TheoryDiffCross->Draw();

  /** TEMP **/
  if(loud){
    cout << "UNSCALED THEORY DIFF CROSS SECTION EVALUATED AT DATA POINTS:::" << endl;
    cout << setprecision(6);
    for(int i=0; i<numAngleBins; i++){
      cout << anglecentres.at(i) << "\t" << TheoryDiffCross->Eval(anglecentres.at(i)) << endl;
    }
    cout << setprecision(3);
    cout << "......................" << endl;
  }
  /** **** **/

  /* Using Chi2 minimizaiton */
  if(loud){cout << "USING CHI2 MINIMIZAITON..." << endl;}
  TCanvas* c_Chi2Min = new TCanvas("c_Chi2Min","c_Chi2Min",1000,1000);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.03);

  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad1->SetLogy();
  pad1->SetTickx();
  pad1->SetTicky();
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);
  pad2->SetTickx();
  pad2->SetTicky();
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  TGraph* Final = FindNormalisation(TheoryDiffCross,gdSdO);
  gdSdO->SetLineColor(kRed);
  gdSdO->SetMarkerColor(kRed);
  gdSdO->SetMarkerStyle(21);
  /* Construct file name string */
  /**/  ostringstream tempstream;
  /**/  if(means_CS.at(indexE)<1.0){tempstream << 0;}
  /**/  tempstream << (int) (means_CS.at(indexE)*1000);
  /**/  tempstream << "_" << orbital; 
  /**/  tempstream << "_spin" << Spin;
  /**/  string tempstr = tempstream.str();
  /* Construct hist title string */
  /**/  ostringstream textstream;
  /**/  textstream << std::fixed << setprecision(3);
  /**/  textstream << "   " << means_CS.at(indexE);
  /**/  textstream << " MeV, ";
  /**/  textstream <<  orbitalname;
  /**/  textstream << ", spin " << (int)Spin;
  /**/	textstream << " --> S = " << globalS 
  /**/	           << " +- " << globalSerr;
  /**/  string textstring = textstream.str(); 

  gdSdO->SetTitle(textstring.c_str());
  gdSdO->GetYaxis()->SetTitleOffset(1.3);
  gdSdO->GetYaxis()->SetTitleSize(0.042);
  gdSdO->GetXaxis()->SetRangeUser(CSangleL,CSangleH); 
  gdSdO->GetYaxis()->SetRangeUser(5e-6,5e-3); //Same vertical range for all? WNC suggestion
  gdSdO->GetXaxis()->SetNdivisions(512,kTRUE);
  gdSdO->Draw("AP");
  Final->Draw("SAME");

  pad2->cd();
  TGraphErrors* gResid = new TGraphErrors(*gdSdO);
  for(int n=0; n < gResid->GetN(); n++){
    double x = gdSdO->GetPointX(n);
    double residual = gdSdO->GetPointY(n) - Final->Eval(x);
    gResid->SetPoint(n,x,residual);
    gResid->SetPointError(n,0,gdSdO->GetErrorY(n));
  }
  TLine* markzero = new TLine(CSangleL,0.,CSangleH,0.);
  gResid->SetTitle("");
  gResid->GetXaxis()->SetRangeUser(CSangleL,CSangleH);
  gResid->GetXaxis()->SetNdivisions(512,kTRUE);
  gResid->GetYaxis()->SetTitle("Residuals");
  gResid->GetYaxis()->SetTitleSize(0.15);
  gResid->GetYaxis()->SetTitleOffset(0.36);
  gResid->GetYaxis()->SetLabelSize(0.08);
  gResid->GetYaxis()->SetNdivisions(305);
  gResid->GetXaxis()->SetTitleSize(0.15);
  gResid->GetXaxis()->SetTitleOffset(0.8);
  gResid->GetXaxis()->SetLabelSize(0.1);
  gResid->GetXaxis()->SetTickLength(0.1);
  gResid->Draw();
  markzero->SetLineStyle(2);
  markzero->Draw("same");

  string savestring1 = "./CS2_Figures/"+tempstr+".root";
  string savestring2 = "./CS2_Figures/"+tempstr+".pdf";
  c_Chi2Min->SaveAs(savestring1.c_str());
  c_Chi2Min->SaveAs(savestring2.c_str());

  cout << YELLOW
       << " -----  C2S = (2J+1)S = (2*" << (int)Spin << "+1)S = (" 
       << (int)((2.*Spin)+1.) << ")S = " 
       << ((2.*Spin)+1.) * globalS << "  ----- " 
       << RESET << endl;

  //delete c_AoSA;
  //delete c_dSdO;
}

////////////////////////////////////////////////////////////////////////////////
vector<vector<double>> GetExpDiffCross(double Energy){
  cout << "========================================================" << endl;
  vector<vector<double>> AllPeaks_OneGate;
  vector<vector<double>> OnePeak_AllGates;
  /****CHANGE ANGLE GATING****/
  //double widthAngleBins = 2.5;
  //double firstAngle = 105.;
  /***************************/
  double x[numAngleBins], y[numAngleBins];
  //TList* list = new TList();

  double trackScale = 0.0;
  if(reactionName=="47K(d,p)"){
    /* Determine scaling factor for PhaseSpace */
    TCanvas* c_ExSubPSpace = new TCanvas("c_ExSubPSpace","c_ExSubPSpace",1000,1000);
    if(scaleTogether){
      TH1F* baseEx = PullThetaLabHist(0,firstAngle,widthAngleBins);
      TH1F* basePS = PullPhaseSpaceHist(0,firstAngle,widthAngleBins);
      for(int i=1; i<numAngleBins;i++){
        TH1F* addEx = PullThetaLabHist(i,firstAngle,widthAngleBins); baseEx->Add(addEx,1.);
        TH1F* addPS = PullPhaseSpaceHist(i,firstAngle,widthAngleBins); basePS->Add(addPS,1.);
      }
  
      /* Subtract flat background equal to smallest bin in range */
      baseEx->GetXaxis()->SetRange(baseEx->FindBin(-1.),baseEx->FindBin(1.));
      double minValueInRange = baseEx->GetBinContent(baseEx->GetMinimumBin());
      baseEx->GetXaxis()->UnZoom();
      cout << "Subtracting background of " << minValueInRange << endl;
      for(int b=1; b<baseEx->GetNbinsX() ; b++){
        baseEx->SetBinContent(b,baseEx->GetBinContent(b)-minValueInRange);
      }
  
      // Fix error where basePS = 600, baseEx=150
      if(loud){
        cout << "baseEx bins: " << baseEx->GetNbinsX() << endl;
        cout << "basePS bins: " << basePS->GetNbinsX() << endl;
      }
      if(basePS->GetNbinsX()!=baseEx->GetNbinsX()){
        cout << basePS->GetNbinsX() / baseEx->GetNbinsX() << endl;
        int ratio = basePS->GetNbinsX() / baseEx->GetNbinsX();
        basePS->Rebin(ratio);
      }
      if(loud){
        cout << "baseEx bins: " << baseEx->GetNbinsX() << endl;
        cout << "basePS bins: " << basePS->GetNbinsX() << endl;
      }

      /* Begin scaling within range, track changes */
      basePS->Scale(0.1);
      trackScale = 0.1;
      int numAngleBinsScale = baseEx->GetNbinsX();
      int nbinlow = basePS->FindBin(4.); int nbinhigh = basePS->FindBin(8.0);
      for(int b=nbinlow; b<nbinhigh; b++){
        if(baseEx->GetBinContent(b) > 0.0 && basePS->GetBinContent(b) > baseEx->GetBinContent(b)){
  	while(basePS->GetBinContent(b) > baseEx->GetBinContent(b)){
            basePS->Scale(0.99999);
            trackScale *= 0.99999;
          }
        }
      }

      baseEx->Add(basePS,-1.);
      baseEx->SetName("ExSubPSpace");
      baseEx->SetTitle("ExSubPSpace");
      baseEx->Draw();
      if(loud){cout << "PhaseSpace -> ExpData scaling = " << trackScale << endl;}
    }
  }

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  if(reactionName=="47K(d,p)"){
    if(means_CS.at(indexE)>3.5){
      //Find max angle, then round down to nearest integer
      double thetaMax = -14.36*means_CS.at(indexE) + 205.50;
      thetaMax = floor(thetaMax-5.0);//TEMP REMOVE MORE!!!!
      if(loud){cout << RED << thetaMax << " rounds by " << widthAngleBins << " to ";}
    
      //Round down to nearest angular bin (i.e. 2.5deg sections)
      thetaMax = floor((double)thetaMax/widthAngleBins)*widthAngleBins;
      if(loud){cout << thetaMax << RESET << endl;}

      numAngleBins = (int)(thetaMax - firstAngle)/2.5;
    }
  } 
  else if(reactionName=="47K(d,t)"){
    double thetaMax = 0.23*pow(means_dt[indexE],2) 
	            + 0.76*means_dt[indexE]
		    + 9.01;
    thetaMax = floor(thetaMax);

    double thetaMin = 0.10*pow(means_dt[indexE],2) 
	            + 0.06*means_dt[indexE]
		    + 1.58;
    thetaMin = ceil(thetaMin);
    
    numAngleBins = (int)thetaMax - (int)thetaMin;
    firstAngle = (int)thetaMin;

  } 

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


  for(int i=0; i<numAngleBins;i++){
    double min = firstAngle + (i*widthAngleBins);
    double max = min + widthAngleBins;
    cout << "===================================" << endl;
    cout << "min: " << min << " max: " << max << endl;
  
    stringstream tmp; tmp << fixed << setprecision(0); 
    tmp << "c_peakFits_" << min << "_" << max; 
    string tmp2 = tmp.str();
    TCanvas* c_peakFits = new TCanvas("c_peakFits",tmp2.c_str(),1000,1000);

    /* Retrieve theta-gated Ex TH1F from file GateThetaLabHistograms.root */
    /* To change angle gates, run GateThetaLab_MultiWrite() */
    TH1F *gate, *pspace;
    if(reactionName=="47K(d,p)"){
      gate = PullThetaLabHist(i,firstAngle,widthAngleBins);
      pspace = PullPhaseSpaceHist(i,firstAngle,widthAngleBins);

      /* Scale the Phase Space at this angle... */
      /* ... for all angles together */
      if(scaleTogether){
        if(loud){ cout << "gate bins: " << gate->GetNbinsX() << endl;
                  cout << "pspc bins: " << pspace->GetNbinsX() << endl; }
        // Fix error where basePS = 600, baseEx=150
        if(pspace->GetNbinsX()!=gate->GetNbinsX()){
          if(loud){cout << pspace->GetNbinsX() / gate->GetNbinsX() << endl;}
          int ratio = pspace->GetNbinsX() / gate->GetNbinsX();
          pspace->Rebin(ratio);
        }
        if(loud){ cout << "gate bins: " << gate->GetNbinsX() << endl;
                  cout << "pspc bins: " << pspace->GetNbinsX() << endl; }
        gate->Add(pspace,-trackScale);
      } 
      /* ... or seperately for each angular bin */
      /* NOTE THAT THIS DOES NOT ACCOUNT FOR FLAT BACKGROUND */
      else {
        if(pspace->Integral() > 50.){ // Non-garbage histogram
          pspace->Scale(0.01);
  	  trackScale=0.01;
          int numAngleBins = gate->GetNbinsX();
          for(int b=0; b<numAngleBins; b++){
	    if(loud){cout << " FROM " << pspace->GetBinContent(b) << 
	  	              " > " << gate->GetBinContent(b); 
	    }
            while(pspace->GetBinContent(b) > gate->GetBinContent(b)){
              pspace->Scale(0.9999);
	      trackScale*=0.9999;
	    }
	    if(loud){cout << " TO " << pspace->GetBinContent(b) << 
	    	             " > " << gate->GetBinContent(b) << endl;
	    }
          }
          if(loud){cout << " !!! SCALE FOR THIS ANGLE = " << trackScale << endl;}
          gate->Add(pspace,-1);
        }
      }
     /* Retrieve array containing all fits, for one angle gate. *
     * Specific peak of interest selected from the vector by   *
     * global variable indexE                                  */
     AllPeaks_OneGate = FitKnownPeaks_RtrnArry(gate, 0.0); if(loud){cout << "!!!!!!! NO SLIDING SHIFT!!!!!" << endl;}
    }
    else if (reactionName=="47K(d,t)"){
      gate = PullThetaCMHist(i,firstAngle,widthAngleBins);
      AllPeaks_OneGate = FitKnownPeaks_dt_RtrnArry(gate, 0.0); if(loud){cout << "!!!!!!! NO SLIDING SHIFT!!!!!" << endl;}
    }
    
    /* Write PS-subtracted spectrum to list */
    //list->Add(gate);
    //list->Add(c_peakFits);
    gate->GetXaxis()->SetRangeUser(-1.0,+8.0);
    string filename = "./CS2_Figures/";
    filename += tmp2 + ".root";
    c_peakFits->SaveAs(filename.c_str());
    //auto tempfile = new TFile(filename.c_str(),"UPDATE"); //Reopen newly made file
    //auto templist = new TList();
    //templist->Add();


    /* Check correct OneGate vector is selected */
    if(loud){
      cout << "area of " << means_CS.at(indexE) << " = "
           << AllPeaks_OneGate[indexE][1] 
	   << " +- " << AllPeaks_OneGate[indexE][2] 
	   << endl;
    }

    /* Add min and max angle to end of relevant OneGate vector */
    AllPeaks_OneGate[indexE].push_back(min);
    AllPeaks_OneGate[indexE].push_back(max);

    /* Push relevant OneGate vector to end of AllGates vector */
    OnePeak_AllGates.push_back(AllPeaks_OneGate[indexE]);
    delete c_peakFits; 
  }

  /* Write PS-subtracted spectrum to file */
  //TFile* outfile = new TFile("GateThetaLab_ExMinusGatePhaseSpace.root","RECREATE");
  //list->Write("GateThetaLab_ExMinusPhaseSpace",TObject::kSingleKey);

  return OnePeak_AllGates;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullThetaLabHist(int i, double minTheta, double gatesize){
  //TFile* file = new TFile("GateThetaLabHistograms_47Kdp_18Oct22_bin0p2.root","READ");
  string name = "GateThetaLabHistograms_47Kdp_18Oct22_";
  if(!GammaGate){
    name = name + "bin0p2.root";
  } else {
    name = name + to_string((int)(globGammaGate*1000))
	        + "GammaGate.root";
  }
  TFile* file = new TFile(name.c_str(),"READ");
  
  string histname = "cThetaLabGate_" 
	          + to_string((int)(minTheta+(i*gatesize))) + "-" 
		  + to_string((int)(minTheta+((i+1)*gatesize)));
  cout << "Loading " << histname << endl;
  TList *list = (TList*)file->Get("GateThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
  return hist;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullThetaCMHist(int i, double minTheta, double gatesize){
  //TFile* file = new TFile("GateThetaCMHistograms_47Kdt_18Oct22_bin0p1.root","READ");
  string name = "GateThetaCMHistograms_47Kdt_18Oct22_";
  if(!GammaGate){
    name = name + "bin0p1.root";
  } else {
    name = name + to_string((int)(globGammaGate*1000))
	        + "GammaGate.root";
  }
  TFile* file = new TFile(name.c_str(),"READ");

  string histname = "cThetaCMGate_" 
	          + to_string((int)(minTheta+(i*gatesize))) + "-" 
		  + to_string((int)(minTheta+((i+1)*gatesize)));
  cout << "Loading " << histname << endl;
  TList *list = (TList*)file->Get("GateThetaCMHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
  return hist;
}

////////////////////////////////////////////////////////////////////////////////
TH1F* PullPhaseSpaceHist(int i, double minTheta, double gatesize){
  //TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_ReadMe.root","READ");
  //TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_2p5degAngles_26May22v2.root","READ");
  //TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_11Jul22.root","READ");
  TFile* file = new TFile("GatePhaseSpaceThetaLabHistograms_29Aug22_TrueStripRemoval_0p05.root","READ");
  string histname = "cPSpaceThetaLabGate_" 
	          + to_string((int) (minTheta+(i*gatesize))) + "-" 
		  + to_string((int) (minTheta+((i+1)*gatesize)));
  cout << "Loading " << histname << endl;
  TList *list = (TList*)file->Get("GatePhaseSpaceThetaLabHistograms");
  TH1F* hist = (TH1F*)list->FindObject(histname.c_str());
  file->Close();
  return hist;
}

////////////////////////////////////////////////////////////////////////////////
void Scale(TGraph* g , TGraphErrors* ex){
  double scale;
  double mean = 0 ;
  double* eX = ex->GetX();
  double* eY = ex->GetY();
  double totalW = 0;
  double W = 0;
  for(int i = 0 ; i < ex->GetN() ; i++){
    if(eY[i]>1 && eY[i] <1e4){
      // Incremental Error weighted average
      W = 1./ex->GetErrorY(i);
      scale = eY[i]/g->Eval(eX[i]);
      totalW +=W;
      mean = mean + W*(scale - mean)/(totalW);
    }
  }

  //scaleTWOFNR = mean;
  if(loud){
    cout << "SCALED THEORY BY " << mean << endl;
    cout << " therefore S = " << 1/mean << " ??" << endl;
  }
  
  double* x = g->GetX();
  double* y = g->GetY();
  for(unsigned int i = 0 ; i < g->GetN() ; i++)
    g->SetPoint(i,x[i],y[i]*mean);
}

////////////////////////////////////////////////////////////////////////////////
TGraph* TWOFNR(double E, double J0, double J, double n, double l, double j){
  /* This function moves between directories in order to run TWOFNR in proper *
   * location. This is, weirdly, the least tempremental way of doing this.    */

  cout << "========================================================" << endl;
  
  double QValue, BeamEnergy =  7.7;
  int johnson, tandyval, modelA, modelB;
  string njjj;
  if(reactionName=="47K(d,p)"){
    cout << "Using Johnson-Soper ..."; johnson=5; tandyval=0;
    cout << " ... and Chapel-Hill." << endl; modelA=2; modelB=2;
    //QValue = 2.274 - E;
    QValue = 2.419 - E;
    njjj.append("24.jjj");
  } else if (reactionName=="47K(d,t)"){
    QValue = -2.112 - E;
    johnson=1; //just makess it work below with the dp/dt matching
    njjj.append("21.jjj");
  }

  char origDirchar[200];
  getcwd(origDirchar,200);
  string origDir{origDirchar};
  string twofnrDir = "/home/charlottepaxman/Programs/TWOFNR";
  cout << "Current directory    " << origDir << endl;
  cout << "Moving to directory  " << twofnrDir << endl;
  chdir(twofnrDir.c_str());
  //Check
  system("pwd"); 
  cout << "===================================" << endl;

  /* Delete existing tran.jjj & 24.jjj files */
  remove("tran.jjj");
  remove(njjj.c_str());

  std::ofstream Front_Input("in.front");
  Front_Input << "jjj" << std::endl;
  Front_Input << "pipo" << std::endl;
  if(reactionName=="47K(d,p)"){
    Front_Input << 2 << std::endl;
  } else if (reactionName=="47K(d,t)"){
    Front_Input << 5 << std::endl;
  }
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << BeamEnergy << std::endl;
  Front_Input << 47 << " " << 19 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << "0 0 0" << std::endl;
  Front_Input << l << " " << j << std::endl;
  Front_Input << n << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << QValue << std::endl; 
  Front_Input << 1 << std::endl;
  Front_Input << J0 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << johnson << std::endl;
  if(johnson==6){//JTandy selected, give version
    Front_Input << tandyval << std::endl;
  }
  Front_Input << 1 << std::endl;
  Front_Input << J << std::endl;
  if(reactionName=="47K(d,p)"){
    Front_Input << 1 << std::endl;
    Front_Input << modelA << std::endl;
    Front_Input << modelB << std::endl;
    Front_Input << 1 << std::endl;
    Front_Input << 1 << std::endl;
    Front_Input << 1 << std::endl;
  } else if (reactionName=="47K(d,t)"){
    Front_Input << 1 << std::endl;
    Front_Input << 2 << std::endl;
    Front_Input << 1 << std::endl;
    Front_Input << 1 << std::endl;
  }
  Front_Input << 1.25 << " " << 0.65 << std::endl;
  Front_Input << 6.0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;

  Front_Input.close();

  cout << "Filled front20 input file." << endl;
  cout << "Executing front20..." << endl;
  system("(exec ~/Programs/TWOFNR/front20 < in.front > /dev/null)"); 
    ifstream checkfront("tran.jjj");
    if(!checkfront){
      //front20 fail!
      cout << " !! ERROR !! \n front20 failed to complete" << endl;
      return 0;
    } else {
      cout << "-> tran.jjj generated (success not guaranteed)" << endl;
      checkfront.close();
    }

  /* in.twofnr instructs twofnr20 to evaluate tran.jjj */
  cout << "Executing twofnr20..." << endl;
  system("(exec ~/Programs/TWOFNR/twofnr20 < in.twofnr > /dev/null)");
    ifstream checktwofnr(njjj.c_str());
    if(!checktwofnr){
      //twofnr20 fail!
      cout << " !! ERROR !! \n twofnr20 failed to complete" << endl;
      terminate();
//      return 0;
    } else {
      cout << "-> twofnr20 complete!" << endl;
      checktwofnr.close();
    }

  TGraph* CS = new TGraph(njjj.c_str());

  //mb/sr->mb/msr is x1/1000
  for(int i=0; i<CS->GetN(); i++){
    double x, newy;
    CS->GetPoint(i,x,newy);
    newy=newy/1000;
    CS->SetPoint(i,x,newy);
  }


  cout << "===================================" << endl;
  cout << "Current directory    " << twofnrDir << endl;
  cout << "Moving to directory  " << origDir << endl;
  chdir(origDir.c_str());
  system("pwd"); 
  cout << "========================================================" << endl;

  return CS;
}

////////////////////////////////////////////////////////////////////////////////
double Chi2(TGraph* theory, TGraphErrors* exper){
  double Chi2 = 0;
  double chi = 0;

  //cout << setprecision(8);
  //for(int i = 1 ; i < exper->GetN() ; i++){
  for(int i = 0 ; i < exper->GetN() ; i++){
    if(exper->GetPointY(i)>1.0e-8){ //0){
      //chi=(exper->Eval(anglecentres[i])-theory->Eval(anglecentres[i]) ) / (exper->GetErrorY(i));
      chi=(exper->GetPointY(i) - theory->Eval(anglecentres[i]) ) / (exper->GetErrorY(i));
      //cout << "COMPARE::::: " << exper->Eval(anglecentres[i]) << " TO " << exper->GetPointY(i) << endl;
      Chi2 +=chi*chi;
    }
  }
  if(loud){cout << "Chi2 = " << Chi2 << endl;}
  return Chi2;
  //cout << setprecision(3);
}

////////////////////////////////////////////////////////////////////////////////
double ToMininize(const double* parameter){
  static int f = 0 ;
  TGraph* g = new TGraph();
  double* X = currentThry->GetX(); // gets valies from global g1 = tgraph passed to find norm
  double* Y = currentThry->GetY();
  for(int i = 0 ; i < currentThry->GetN() ; i++){
    g->SetPoint(g->GetN(),X[i],parameter[0]*Y[i]); // scales this tgraph by parameter
  }
  double chi2  = Chi2(g,staticExp);  //compares theory tgraph to experimental tgrapherrors by chi2
  return chi2;
}

////////////////////////////////////////////////////////////////////////////////
TGraph* FindNormalisation(TGraph* theory, TGraphErrors* experiment){
  /* (dSdO)meas = S * (dSdO)calc */

  // Set global variable
  currentThry = theory;
  staticExp = experiment;

  // Construct minimiser
  const char* minName ="Minuit";const char* algoName="Migrad";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // Number of parameters (should only be 1 for me)
  int mysize = 1;

  // Create funciton wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(&ToMininize,mysize);
  min->SetFunction(f);
  
  // Set range of parameter(??)
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = 0.8;
    char name[4];
    sprintf(name,"S%d",i+1);
    min->SetLimitedVariable(i,name,parameter[i],0.01,0,10);
  }
 

  ///// TO IMPROVE: FIND WAY OF OBTAINING NDF AND PRINT CHI2/NDF /////

  // Minimise
  min->Minimize();
  const double *xs = min->X();
  const double *err = min->Errors(); 

  // Write out
  for(int i = 0  ; i < mysize ; i++){
    cout << Form("S%d=",i+1) << xs[i] << "(" << err[i] << ")" << endl;
  }

  /* Store S value in global variable, to access for drawing on plots */
  globalS = xs[0];
  globalSerr = err[0];

  // Return the Fitted CS
  TGraph* g = new TGraph(); 
  double* X = theory->GetX();
  double* Y = theory->GetY();
  if(loud){
    cout << setprecision(8);
    cout << "Start: X[0] = " << theory->GetPointX(4) << " Y[0] = " << theory->GetPointY(4) << endl;
    cout << "multip by " << xs[0] << endl;
  }
  
  //for(int i=0; i<theory->GetN(); i++){ g->SetPoint(g->GetN(),X[i],xs[0]*Y[i]); }
  for(int i=0; i<theory->GetN(); i++){ g->SetPoint(i,X[i],xs[0]*Y[i]); }

  if(loud){
    cout << "End:   X[0] = " << g->GetPointX(4) << " Y[0] = " << g->GetPointY(4) << endl;
    cout << setprecision(3);
  }

  return g;
}

////////////////////////////////////////////////////////////////////////////////
void CS(double Energy, double Spin, double spdf, double angmom, double GammaGateEnergy){
  GammaGate=true;
  globGammaGate = GammaGateEnergy;

  CS(Energy, Spin, spdf, angmom);
}
