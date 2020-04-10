/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  F. Flavigny contact address: flavigny@lpccaen.in2p3.fr *
 *                                                                           *
 * Creation Date   : April 2019                                              *
 * Last update     : Nov 2019                                                *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class deal with Quasi Free Scattering Reaction in which a cluster   *
 *  or a nucleon is removed from a projectile  by interaction with a target  *
 *  nucleon (proton target in general)                                       *
 *                                                                           *
 *  First step (dissociation):  A -> B + c                                   *
 *  Second step (scattering) :  c + T -> 1 + 2                               *
 *  Labeling is:                                                             *
 *                                                                           *
 *              A --> T  ==> B + (c -> T) =>  B + 1 + 2                      *
 *                                                                           *
 *  where:                                                                   *
 *    +  A is the beam nucleus                                               *
 *    +  T is the target nucleon (proton)                                    *
 *                                                                           *
 *    +  B is the residual fragment (beam-like)                              *
 *    +  1 is the scattered target nucleon  (former T)                       *
 *    +  2 is the knocked-out cluster/nucleon (noted c) in the intermediate  *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *    +  Adapted from original event generator from V. Panin (R3B collab)    *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>

#include "NPQFS.h"
#include "NPCore.h"
#include "NPOptionManager.h"
#include "NPFunction.h"

// Use CLHEP System of unit and Physical Constant
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// ROOT
#include"TF1.h"
#include"TRandom3.h"

ClassImp(QFS)

////////////////////////////////////////////////////////////////////////////////

QFS::QFS(){

    //------------- Default Constructor -------------
    fVerboseLevel         = NPOptionManager::getInstance()->GetVerboseLevel();
    fBeamEnergy           = 0;
    fThetaCM              = 0;
    fPhiCM                = 0;
 
    fExcitationA          = 0;
    fExcitationB          = 0;
    fMomentumSigma        = 0;
    fshootB=false;
    fshoot1=true;
    fshoot2=true;
    fisotropic = true;

    fTheta2VsTheta1 = 0;
    fPhi2VsPhi1 = 0;

}

////////////////////////////////////////////////////////////////////////////////

QFS::~QFS(){

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void QFS::ReadConfigurationFile(string Path){
  ifstream ReactionFile;
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/" + Path;
  ReactionFile.open(Path.c_str());
  if (!ReactionFile.is_open()) {
    ReactionFile.open(StandardPath.c_str());
    if(ReactionFile.is_open()) {
      Path = StandardPath;
    }
    else {cout << "QFS File " << Path << " not found" << endl;exit(1);}
  }
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}

////////////////////////////////////////////////////////////////////////////////

void QFS::ReadConfigurationFile(NPL::InputParser parser){

  cout << " In QFS ReadConfiguration " << endl;
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("QFSReaction");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
      cout << endl << "\033[1;35m//// QFS reaction found " << endl;

  vector<string> token1 = {"Beam","Target","Scattered","KnockedOut","Heavy"};
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
      if(blocks[i]->HasTokenList(token1)){
          int v = NPOptionManager::getInstance()->GetVerboseLevel();
          NPOptionManager::getInstance()->SetVerboseLevel(0);
          fNucleiA.ReadConfigurationFile(parser);
          NPOptionManager::getInstance()->SetVerboseLevel(v);

          fBeamEnergy= fNucleiA.GetEnergy();
          GetNucleus(blocks[i]->GetString("Beam"),parser);
          fNucleiT = GetNucleus(blocks[i]->GetString("Target"),parser);
          fNucleiB = GetNucleus(blocks[i]->GetString("Heavy"),parser);
          fNuclei1 = GetNucleus(blocks[i]->GetString("Scattered"),parser);
          fNuclei2 = GetNucleus(blocks[i]->GetString("KnockedOut"),parser);
      }
      else{
          cout << "ERROR: check your input file formatting \033[0m" << endl;
          exit(1);
      }
      if(blocks[i]->HasToken("ExcitationEnergyBeam")){
          fExcitationA = blocks[i]->GetDouble("ExcitationEnergyBeam","MeV");
      }
      if(blocks[i]->HasToken("ExcitationEnergyHeavy")){
          fExcitationB = blocks[i]->GetDouble("ExcitationEnergyHeavy","MeV");
      }
      if(blocks[i]->HasToken("MomentumSigma")){
          fMomentumSigma = blocks[i]->GetDouble("MomentumSigma","MeV");
      }
      if(blocks[i]->HasToken("ShootHeavy")){
          fshootB = blocks[i]->GetInt("ShootHeavy");
      }
      if(blocks[i]->HasToken("ShootLight")){
          fshoot1 = blocks[i]->GetInt("ShootLight");
          fshoot2 = blocks[i]->GetInt("ShootLight");
      }
  }

  cout << "\033[0m" ;
}

////////////////////////////////////////////////////////////////////////////////
Nucleus QFS::GetNucleus(string name, NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("DefineNucleus",name);
  unsigned int size = blocks.size();
  if(size==0)
    return NPL::Nucleus(name);
  else if(size==1){
    cout << " -- User defined nucleus " << name << " -- " << endl;
    vector<string> token = {"SubPart","BindingEnergy"};
    if(blocks[0]->HasTokenList(token)){
      NPL::Nucleus N(name,blocks[0]->GetVectorString("SubPart"),blocks[0]->GetDouble("BindingEnergy","MeV"));
      if(blocks[0]->HasToken("ExcitationEnergy"))
        N.SetExcitationEnergy(blocks[0]->GetDouble("ExcitationEnergy","MeV"));
      if(blocks[0]->HasToken("SpinParity"))
        N.SetSpinParity(blocks[0]->GetString("SpinParity").c_str());
      if(blocks[0]->HasToken("Spin"))
        N.SetSpin(blocks[0]->GetDouble("Spin",""));
      if(blocks[0]->HasToken("Parity"))
        N.SetParity(blocks[0]->GetString("Parity").c_str());
      if(blocks[0]->HasToken("LifeTime"))
        N.SetLifeTime(blocks[0]->GetDouble("LifeTime","s"));

    cout << " -- -- -- -- -- -- -- -- -- -- --" << endl;
      return N;
    }
  }
  else{
    NPL::SendErrorAndExit("NPL::QFS","Too many nuclei define with the same name");
  }
}


////////////////////////////////////////////////////////////////////////////////////////////
void QFS::CalculateVariables(){

  if(fBeamEnergy < 0)
    fBeamEnergy = 0 ; 

    //cout<<"---- COMPUTE ------"<<endl;
   // cout<<"--CM--"<<endl; 

    mA = fNucleiA.Mass();            // Beam mass in MeV
    mT =  fNucleiT.Mass();           // Target mass in MeV 
    mB =  fNucleiB.Mass();           // Heavy residual mass in MeV 
    m1 =  mT;                        // scattered target nucleon (same mass);
    m2 =  fNuclei2.Mass();           // knocked out cluster mass in MeV 
    ma =  m2;                        // intermediate cluster mass in MeV (same);
 
    double TA = fBeamEnergy;                 // Beam kinetic energy
    double PA = sqrt(TA*(TA+2*mA));          // Beam momentum (norm)
    double EA = sqrt(mA*mA + PA*PA);         // Beam total energy
    fEnergyImpulsionLab_A = TLorentzVector(0.,0.,PA,EA);
    
    //Internal momentum of removed cluster/nucleon
    static TRandom3 r;
    //r.SetSeed(0);
    //double mom_sigma = 0; // MeV/c
    Pa.SetX(r.Gaus(0.,fMomentumSigma));
    Pa.SetY(r.Gaus(0.,fMomentumSigma));
    Pa.SetZ(r.Gaus(0.,fMomentumSigma));

    //Internal momentum of heavy recoil after removal
    PB.SetXYZ( (-Pa.X()) , (-Pa.Y()) , (-Pa.Z()) );

    //cout<<"Pa_cm=\t("<<Pa.Px()<<","<<Pa.Py()<<","<<Pa.Pz()<<")"<<endl;
    //cout<<"||Pa||^2=\t("<<Pa.Mag2()<<endl;
    //cout<<"mA^2=\t"<<mA*mA<<endl;
    //cout<<"mB^2=\t"<<mB*mB<<endl;
 

    // Off-shell mass of the bound nucleon from E conservation
    // in virtual dissociation of A -> B + a
    double buffer = mA*mA + mB*mB - 2*mA*sqrt(mB*mB+Pa.Mag2()) ; 
    if(buffer<=0) { cout<<"ERROR off shell mass ma_off=\t"<<buffer<<endl; return;}
    ma_off = sqrt(buffer);

    //deduced total energies of "a" and "B" in restframe of A
    double Ea = sqrt(ma_off*ma_off + Pa.Mag2());
    double EB = sqrt(mB*mB + PB.Mag2());

    //cout<<"ma_off^2=\t"<<buffer<<endl;
    //cout<<"Ea=\t"<<Ea<<endl;
    //cout<<"EB=\t"<<EB<<endl;

    fEnergyImpulsionLab_a = TLorentzVector(Pa,Ea);
    fEnergyImpulsionLab_B = TLorentzVector(PB,EB);
    fEnergyImpulsionLab_a.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    fEnergyImpulsionLab_B.Boost(0,0,fEnergyImpulsionLab_A.Beta());
    Ea_lab = fEnergyImpulsionLab_a.E();
    EB_lab = fEnergyImpulsionLab_B.E();
    Pa = fEnergyImpulsionLab_a.Vect();
    PB = fEnergyImpulsionLab_B.Vect();
   
    // Scattering part (2-body kinematics)
    // virtual cluster of mass "ma_off" scattering on target T
    // to give scattered  cluster with real mass (ma=m2)
    // and scattered target (mT=m1)

    fQValue =ma_off+mT-m1-m2;

    s = ma_off*ma_off + mT*mT + 2*mT*Ea_lab ; 
    fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s));
    fEcm = sqrt(s) - m1 -m2;
    if(fEcm<=0) { cout<<"ERROR Ecm negative =\t"<<fEcm<<endl; return;}

    ECM_a = (s + ma_off*ma_off - mT*mT)/(2*sqrt(s));
    ECM_T = (s + mT*mT - ma_off*ma_off)/(2*sqrt(s));
    ECM_1 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
    ECM_2 = (s + m2*m2 - m1*m1)/(2*sqrt(s));

    pCM_a = sqrt(ECM_a*ECM_a - ma_off*ma_off);
    pCM_T = sqrt(ECM_T*ECM_T - mT*mT);
    pCM_1 = sqrt(ECM_1*ECM_1 - m1*m1);
    pCM_2 = sqrt(ECM_2*ECM_2 - m2*m2);

    BetaCM = Pa.Mag() / (Ea_lab + mT);

    //if(ECM_a<0 || ECM_T<0 || ECM_1<0 || ECM_2<0 ) { 
    /*
    if(pCM_a<0 || pCM_T<0 || pCM_1<0 || pCM_2<0 ) { 
    cout<<"--LAB after Boost--"<<endl; 
    cout<<"Pa_lab=\t("<<Pa.Px()<<","<<Pa.Py()<<","<<Pa.Pz()<<")"<<endl;
    cout<<"PB_lab=\t("<<PB.Px()<<","<<PB.Py()<<","<<PB.Pz()<<")"<<endl;
    cout<<"Ea_lab=\t"<<Ea_lab<<endl;
    cout<<"EB_lab=\t"<<EB_lab<<endl;
    cout<<"--Back to CM--"<<endl; 
    cout<<"fQValue=\t"<<fQValue<<endl;
    cout<<"s=\t"<<s<<endl;
    cout<<"Ecm=\t"<<fEcm<<endl;
    cout<<"ea*=\t"<<ECM_a<<endl;
    cout<<"eT*=\t"<<ECM_T<<endl;
    cout<<"e1*=\t"<<ECM_1<<endl;
    cout<<"p1*=\t"<<pCM_1<<endl;
    cout<<"e2*=\t"<<ECM_2<<endl;
    cout<<"p2*=\t"<<pCM_2<<endl;
    cout<<"beta_cm=\t"<<BetaCM<<endl;
    }
    */
}

////////////////////////////////////////////////////////////////////////////////////////////

void QFS::KineRelativistic(double& ThetaLab1, double& PhiLab1, double& KineticEnergyLab1, double& ThetaLab2, double& PhiLab2, double& KineticEnergyLab2){

    CalculateVariables();

    //double thetaCM_1 = fThetaCM*TMath::Pi()/180.;
    double thetaCM_1 = fThetaCM;
    double thetaCM_2 =  M_PI - thetaCM_1;
    //double phiCM_1 = fPhiCM*TMath::Pi()/180.;
    double phiCM_1 = fPhiCM;
    double phiCM_2 =  2*M_PI - phiCM_1;

    TVector3 z_axis(0.,0.,1.);

    fEnergyImpulsionCM_2	= TLorentzVector(
                                        pCM_2*sin(thetaCM_2)*cos(phiCM_2),
                                        pCM_2*cos(thetaCM_2)*sin(phiCM_2),
                                        pCM_2*cos(thetaCM_2),
                                        ECM_2);
    fEnergyImpulsionCM_1	= fTotalEnergyImpulsionCM - fEnergyImpulsionCM_2;

    fEnergyImpulsionLab_1 = fEnergyImpulsionCM_1;
    fEnergyImpulsionLab_1.Boost(0,0,BetaCM);
    fEnergyImpulsionLab_2 = fEnergyImpulsionCM_2;
    fEnergyImpulsionLab_2.Boost(0,0,BetaCM);

    // Angle in the lab frame
    //ThetaLab1 = fEnergyImpulsionLab_1.Angle(fEnergyImpulsionLab_A.Vect());
    ThetaLab1 = fEnergyImpulsionLab_1.Angle(z_axis);
    if (ThetaLab1 < 0) ThetaLab1 += M_PI;
    //ThetaLab2 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_A.Vect());
    ThetaLab2 = fEnergyImpulsionLab_2.Angle(z_axis);
    if (fabs(ThetaLab1) < 1e-6) ThetaLab1 = 0;
    ThetaLab2 = fabs(ThetaLab2);
    if (fabs(ThetaLab2) < 1e-6) ThetaLab2 = 0;

    PhiLab1 = M_PI + fEnergyImpulsionLab_1.Vect().Phi(); 
    if (fabs(PhiLab1) < 1e-6) PhiLab1 = 0;
    PhiLab2 = M_PI + fEnergyImpulsionLab_2.Vect().Phi(); 
    if (fabs(PhiLab2) < 1e-6) PhiLab2 = 0;

    // Kinetic Energy in the lab frame
    KineticEnergyLab1 = fEnergyImpulsionLab_1.E() - m1;
    KineticEnergyLab2 = fEnergyImpulsionLab_2.E() - m2;

    // test for total energy conversion
    //if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_1.E()+fEnergyImpulsionLab_2.E())) > 1e-6)
    //    cout << "Problem for energy conservation" << endl;

/*
    cout<<"--KINE RELATIVISTIC--"<<endl;
    cout<<"theta_cm:"<<fThetaCM*180./TMath::Pi()<<endl;
    cout<<"phi_cm:"<<fPhiCM*180./TMath::Pi()<<endl;
    cout<<"P1_CM=\t("<<fEnergyImpulsionCM_1.Px()<<","<<fEnergyImpulsionCM_1.Py()<<","<<fEnergyImpulsionCM_1.Pz()<<")"<<endl;
    cout<<"P2_CM=\t("<<fEnergyImpulsionCM_2.Px()<<","<<fEnergyImpulsionCM_2.Py()<<","<<fEnergyImpulsionCM_2.Pz()<<")"<<endl;
    cout<<"P1_lab=\t("<<fEnergyImpulsionLab_1.Px()<<","<<fEnergyImpulsionLab_1.Py()<<","<<fEnergyImpulsionLab_1.Pz()<<")"<<endl;
    cout<<"P2_lab=\t("<<fEnergyImpulsionLab_2.Px()<<","<<fEnergyImpulsionLab_2.Py()<<","<<fEnergyImpulsionLab_2.Pz()<<")"<<endl;
    cout<<"Theta1:\t"<<fEnergyImpulsionLab_1.Vect().Theta()*180./TMath::Pi()<<endl;
    cout<<"Theta2:\t"<<fEnergyImpulsionLab_2.Vect().Theta()*180./TMath::Pi()<<endl;
    cout<<"Phi1:\t"<<M_PI+fEnergyImpulsionLab_1.Vect().Phi()*180./TMath::Pi()<<endl;
    cout<<"Phi2:\t"<<M_PI+fEnergyImpulsionLab_2.Vect().Phi()*180./TMath::Pi()<<endl;
*/

}

////////////////////////////////////////////////////////////////////////////////////////////

double QFS::ShootRandomThetaCM(){

  double theta; // CM
/*
  if(fDoubleDifferentialCrossSectionHist){
    // Take a slice in energy
    TAxis* Y = fDoubleDifferentialCrossSectionHist->GetYaxis();
    int binY;

    // Those test are there for the tail event of the energy distribution
    // In case the energy is outside the range of the 2D histo we take the
    // closest available CS
    if(Y->FindBin(fBeamEnergy) > Y->GetLast())
      binY = Y->GetLast()-1;

    else if(Y->FindBin(fBeamEnergy) < Y->GetFirst())
      binY = Y->GetFirst()+1;

    else
      binY = Y->FindBin(fBeamEnergy);

    TH1D* Proj = fDoubleDifferentialCrossSectionHist->ProjectionX("proj",binY,binY);
    SetThetaCM( theta=Proj->GetRandom()*deg );
  }
  else if (fLabCrossSection){
    double thetalab=-1;
    double energylab=-1;
    while(energylab<0){
      thetalab=fCrossSectionHist->GetRandom()*deg; //shoot in lab
      energylab=EnergyLabFromThetaLab(thetalab);   //get corresponding energy
    }
    theta = EnergyLabToThetaCM(energylab, thetalab); //transform to theta CM
    SetThetaCM( theta );
  }
  else{
    // When root perform a Spline interpolation to shoot random number out of
    // the distribution, it can over shoot and output a number larger that 180
    // this lead to an additional signal at 0-4 deg Lab, especially when using a
    // flat distribution.
    // This fix it.
    theta=181;
    if(theta/deg>180)
      theta=fCrossSectionHist->GetRandom();
    //cout << " Shooting Random ThetaCM "  << theta << endl;
    SetThetaCM( theta*deg );
  }
*/
  return theta;
}

////////////////////////////////////////////////////////////////////////////////////////////

TGraph* QFS::GetTheta2VsTheta1(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta1,phi1,E1,theta2,phi2,E2;

  for (double angle=0 ; angle <= 180 ; angle+=AngleStep_CM){
    SetThetaCM(angle*TMath::Pi()/180.);
    KineRelativistic(theta1, phi1, E1, theta2, phi2, E2);

    vx.push_back(theta1*180./M_PI);
    vy.push_back(theta2*180./M_PI);
  }
  fTheta2VsTheta1 = new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fTheta2VsTheta1);
}

////////////////////////////////////////////////////////////////////////////////////////////

TGraph* QFS::GetPhi2VsPhi1(double AngleStep_CM){

  vector<double> vx;
  vector<double> vy;
  double theta1,phi1,E1,theta2,phi2,E2;

  for (double theta=0 ; theta <= 180 ; theta+=AngleStep_CM){
      for (double angle=0 ; angle <= 180 ; angle+=AngleStep_CM){
          SetThetaCM(theta*TMath::Pi()/180.);
          SetPhiCM(angle*TMath::Pi()/180.);
          KineRelativistic(theta1, phi1, E1, theta2, phi2, E2);
          vx.push_back(phi1*180./M_PI);
          vy.push_back(phi2*180./M_PI);
      }
  }
  fPhi2VsPhi1 = new TGraph(vx.size(),&vx[0],&vy[0]);

  return(fPhi2VsPhi1);
}

///////////////////////////////////////////////////////////////////////////////
// Check whenever the reaction is allowed at a given energy
bool QFS::IsAllowed(){//double Energy){
  //double AvailableEnergy = Energy + fNuclei1.Mass() + fNuclei2.Mass();
  //double RequiredEnergy  = fNuclei3.Mass() + fNuclei4.Mass();

  //if(AvailableEnergy>RequiredEnergy)
    return true;
  //else
  //  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
///////////////// Old R3B method not using TLorentz Vector  ////////////////////////////////
/////////////////// (used as a reference)///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/*
void QFS::CalculateVariablesOld(){

    initializePrecomputeVariable();

    vector<double> theta1;
    vector<double> theta2;
    vector<double> phi1;
    vector<double> phi2;

    //for(int i=0; i<=180; i++){
    int i = 29;
        KineR3B(s, ma_off, mT, ma, (double)i);
        if(!good) { cout<<"ERROR CM calculations!!!"<<endl; return;}

        TVector3 P1_cm(0.,0.,1.), P2_cm(0.,0.,1.);
        P2_cm.SetMag(p_clust);
        P2_cm.SetTheta(theta_clust);
        TRandom3 ra;
        ra.SetSeed(0);
        //P2_cm.SetPhi(ra.Uniform(-1*TMath::Pi(),+1*TMath::Pi()));
        P2_cm.SetPhi(0.);
        P1_cm.SetX(-P2_cm.X());
        P1_cm.SetY(-P2_cm.Y());
        P1_cm.SetZ(-P2_cm.Z());
        cout<<"FFFFFFFFFFFFFFFFF"<<endl;
        cout<<"P1_CM=\t("<<P1_cm.X()<<","<<P1_cm.Y()<<","<<P1_cm.Z()<<")"<<endl;
        cout<<"P2_CM=\t("<<P2_cm.X()<<","<<P2_cm.Y()<<","<<P2_cm.Z()<<")"<<endl;
 

        // Calculate relative to direction of quasi-particle (cluster)

        double beta_cm = -Pa.Mag() / (Ea_lab + mT);
        double gamma_cm = 1/sqrt(1-beta_cm*beta_cm);

        pair<double,double> lor_a1 = Lorentz(gamma_cm,beta_cm,e_scat,P1_cm.Z());
        pair<double,double> lor_a2 = Lorentz(gamma_cm,beta_cm,e_clust,P2_cm.Z());

        P1_cm.SetZ(lor_a1.second);
        P2_cm.SetZ(lor_a2.second);

        //Rotating back to beam direction
        //TVector3 P1_L = Rotations(P1_cm, Pa);
        //TVector3 P2_L = Rotations(P2_cm, Pa);
        
        TVector3 P1_L = P1_cm;
        TVector3 P2_L = P2_cm;
        TVector3 direction = Pa.Unit();
        //P1_L.RotateUz(direction);
        //P1_L.RotateUz(direction);

        cout<<"----Calculate variables output------"<<endl;
       cout<<"--CM--"<<endl;
        cout<<"theta1*=\t"<<theta_scat*180/TMath::Pi()<<endl;
        cout<<"theta2*=\t"<<theta_clust*180/TMath::Pi()<<endl;
        cout<<"e1*=\t"<<e_scat<<endl;
        cout<<"p1*=\t"<<p_scat<<endl;
        cout<<"e2*=\t"<<e_clust<<endl;
        cout<<"p2*=\t"<<p_clust<<endl;
        cout<<"T=\t"<<T<<endl;
        cout<<"beta_cm=\t"<<beta_cm<<endl;
        cout<<"gamma_cm=\t"<<gamma_cm<<endl;

        cout<<"--LAB (cluster dir)--"<<endl;
        cout<<"P1_lab=\t("<<P1_cm.X()<<","<<P1_cm.Y()<<","<<P1_cm.Z()<<")"<<endl;
        cout<<"P2_lab=\t("<<P2_cm.X()<<","<<P2_cm.Y()<<","<<P2_cm.Z()<<")"<<endl;
        cout<<"Theta1:\t"<<P1_cm.Theta()*180./TMath::Pi()<<endl;
        cout<<"Theta2:\t"<<P2_cm.Theta()*180./TMath::Pi()<<endl;

        cout<<"--LAB--"<<endl;
        cout<<"Theta1L:\t"<<P1_L.Theta()*180./TMath::Pi()<<endl;
        cout<<"Theta2L:\t"<<P2_L.Theta()*180./TMath::Pi()<<endl;
        cout<<"Phi1L:\t"<<P1_L.Phi()*180./TMath::Pi()<<endl;
        cout<<"Phi2L:\t"<<P2_L.Phi()*180./TMath::Pi()<<endl;


        //cout<<P1_cm.Theta()*180./TMath::Pi()<<"\t"<<P2_cm.Theta()*180./TMath::Pi()<<endl;
        //cout<<P1_L.Phi()*180./TMath::Pi()<<"\t"<<P2_L.Phi()*180./TMath::Pi()<<endl;
        
       theta1.push_back(P1_L.Theta()*180./TMath::Pi());
       theta2.push_back(P2_L.Theta()*180./TMath::Pi());
      
       double temp_phi1 = P1_L.Phi(); 
       double temp_phi2 = P2_L.Phi(); 
       phi1.push_back(180. + P1_L.Phi()*180./TMath::Pi());
       phi2.push_back(180. + P2_L.Phi()*180./TMath::Pi());

    //}
    TGraph* fTheta2VsTheta1 = new TGraph(theta1.size(),&theta1[0],&theta2[0]);
    TGraph* fPhi2VsPhi1 = new TGraph(phi1.size(),&phi1[0],&phi2[0]);
    fTheta2VsTheta1->SetName("Theta2VsTheta1");
    fPhi2VsPhi1->SetName("Phi2VsPhi1");
    TFile* f = new TFile("graphs.root","RECREATE");
    fTheta2VsTheta1->Write();
    fPhi2VsPhi1->Write();
    f->Close();


}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
void QFS::KineR3B(double s,double m2off,double m1,double m2,double thetacm)
{
     if(thetacm>180 || thetacm<0){
        cout << "\nERROR! ThetaCM (in deg) should be between 0 and 180"<<endl;
        return;
    }

	e_clust = 0;
	p_clust = 0;
	theta_clust = 0;
	e_scat = 0;
	p_scat = 0;
	theta_scat = 0;
	T = 0;
	good = false;


	double X = s;
	double Y = m2off*m2off;
	double Z = m1*m1;
	double sqrt_s = sqrt(s);

	// Kinematics before the scattering process
	// (with one off-shell mass)
	double p2_off = sqrt(function(X,Y,Z))/2/sqrt_s;
	double p1_off = p2_off;
	// CM energies
	double e1_off = (s+Z-Y)/2/sqrt_s;
	double e2_off = (s+Y-Z)/2/sqrt_s;

	// Now take the real masses (after scattering)
	Y = m2*m2;  Z = m1*m1;
	//And check whether the kinematical function is ok
	//for this specific kinematical case
	double ERROR_CI = function(X,Y,Z);
	if(ERROR_CI <= 0.){
		cout << "\nERROR!!! Kinematical function is negative!";
		return;
	}

	// Kinematics after the scattering process
	// (with all real masses)
	double p2 = sqrt(function(X,Y,Z))/2/sqrt_s;
	double p1 = p2;
	double e1 = (s+Z-Y)/2/sqrt_s;
	double e2 = (s+Y-Z)/2/sqrt_s;

	// Let's consider momentum transfer <t> from the
	// target particle 1 to the cluster 2
    double t = 2*(m1*m1 - e1_off*e1 + p1_off*p1*cos(thetacm*TMath::Pi()/180.));

    //CM scattering angles
    double theta1 = thetacm*TMath::Pi()/180.;
    double theta2 = TMath::Pi() - theta1;

    e_clust = e2;
    p_clust = p2;
    theta_clust = theta2;

    e_scat = e1;
    p_scat = p1;
    theta_scat = theta1;

    T = t;
    good = true;


}



double QFS::function(double x,double y,double z)
{	
	double lambda = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
	return lambda;
}

//---- Two consecutive rotations 
//first around Z on <phi>, then around new X' on <theta> (1=Pcm, 2=Pa in lab)
TVector3 QFS::Rotations(TVector3 v1,TVector3 v2) 
{
	double CT = v2.Z()/v2.Mag(); // cos(theta) of v2 wrt. Z-axis
	double ST = sqrt(1-CT*CT);   // sin(theta)
	double CF = v2.X()/v2.Mag()/ST;
	double SF = v2.Y()/v2.Mag()/ST;

	TVector3 v3;
	double _v3x =  v1.X()*CT*CF - v1.Y()*SF + v1.Z()*ST*CF;
	double _v3y =  v1.X()*CT*SF + v1.Y()*CF + v1.Z()*ST*SF;
	double _v3z = -v1.X()*ST   +  v1.Z()*CT;
	v3.SetXYZ(_v3x,_v3y,_v3z);
	return v3;
}



pair<double, double> QFS::Lorentz(double gamma,double beta,double e,double p)
{
	double eL = gamma*e - gamma*beta*p;
	double pL = gamma*p - gamma*beta*e;
	return make_pair(eL, pL);
}
*/
