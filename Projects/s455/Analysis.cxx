/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace contact address: pierre.morfouace@cea.fr    *
 *                                                                           *
 * Creation Date  : June 2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Sofia analysis project                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
#include<algorithm>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPPhysicalConstants.h"
#include"NPGlobalSystemOfUnits.h"


////////////////////////////////////////////////////////////////////////////////
struct TofPair
{
  double x = -1000;
  double y = -1000;
  double tof = -1000;
  double velocity = -1;
  double beta = -1;
  double gamma = -1;
  double theta_in = -10;
  double theta_out = -10;
  double psi = -10;
  int plastic = -1;
  //
  int isLorR  = -1;
  // *** isLorR = 1 -> Left
  // *** isLorR = 2 -> Right 
  //
  int isUorD  = -1;
  // *** isUorD = 1 -> Up
  // *** isUorD = 2 -> Down
  int section = -1;
  double Esec = -1;
  double Z = 0;
  int iZ = 0;
  double AoQ = 0;
  double A = 0;
  double DT = -100;
  double x2twim = -1000;
  double x1 = -1000;
  double x2 = -1000;
  double x3 = -1000;
  double y3 = -1000;
  double x3lab = 0;
  double z3lab = 0;
  double xb = 0;
  double xc = 0;
  double xd = 0;
  double yc = 0;
  double zb = 0;
  double zc = 0;
  double zd = 0;
  double flight_path = 0;
  double Leff = 0;
  double rho = 0;
  double Brho = 0;
  double BrhoX = 0;
  double BrhoZ = 0;
  double omega = 0;
  double deff1 = 0;
  double deff2 = 0;
};


////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  SofBeamID = new TSofBeamID();
  SofFF = new TSofFissionFragment();
  SofSci= (TSofSciPhysics*) m_DetectorManager->GetDetector("SofSci");
  SofTrim= (TSofTrimPhysics*) m_DetectorManager->GetDetector("SofTrim");
  SofTwim= (TSofTwimPhysics*) m_DetectorManager->GetDetector("SofTwim");
  SofTofW= (TSofTofWPhysics*) m_DetectorManager->GetDetector("SofTofW");
  SofAt= (TSofAtPhysics*) m_DetectorManager->GetDetector("SofAt");
  SofMwpc= (TSofMwpcPhysics*) m_DetectorManager->GetDetector("SofMwpc");

  ReadAnalysisConfig();
  for(int i=0; i<3; i++){
    // i=0 Cahtode 1 -> Lead 1
    // i=1 Cathode 2 -> Carbon
    // i=2 Cathode 3 -> Lead 2
    fDistancePlasticToCathode[i] = fDistanceStartToFirstATCathode + i*fDistanceBetweenCathode;
    cout << "**** Distance Plastic to cathode " << i+1 << " = " << fDistancePlasticToCathode[i] << endl; 
  }

  InitParameter();
  InitOutputBranch();
  LoadSpline();
  LoadActiveTargetCuts();

  m_GladField = new GladFieldMap();
  m_GladField->SetCurrent(fGladCurrent);
  //m_GladField->SetGladEntrance(0, 0.02*m, 2.774*m + 0.5405*m);
  m_GladField->SetGladEntrance(0, 0, -1.1135*m);
  //m_GladField->SetGladTurningPoint(0, 0.02*m, 2.774*m  + 0.5405*m + 1.1135*m);
  m_GladField->SetGladTurningPoint(0, 0, 0);
  m_GladField->SetGladTiltAngle(-14.*deg);
  m_GladField->LoadMap("GladFieldMap_50mm.dat");
  m_GladField->SetBin(50);
  m_GladField->SetTimeStep(0.8);
  m_GladField->SetCentralTheta(-20.*deg);
  
  double Z_MWPC3 = fDistanceGToMW3;// * cos(m_GladField->GetCentralTheta());
  double X_MWPC3 = (Z_MWPC3 - m_GladField->GetGladTurningPoint().Z())*tan(m_GladField->GetCentralTheta());
  m_GladField->Set_MWPC3_Position(X_MWPC3,0,Z_MWPC3);

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  ReInitValue();
  //cout << "************" << endl;
  RunID = fRunID;
  BeamAnalysis();

  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double start_time = SofSci->TimeNs[1];
    SofTofW->SetTofAlignedValue(36);
    SofTofW->SetStartTime(start_time);
    SofTofW->BuildPhysicalEvent();

    if(SofAt->Energy.size()==3 || SofAt->Energy.size()==4){ 
      double Anode1 = SofAt->Energy[0];
      double Anode2 = SofAt->Energy[1];
      double Anode3 = SofAt->Energy[2];
      int k= fRunID-1;
      int which_cathode = 0;
      
      if(cut_Pb1[k]->IsInside(Anode1,Anode2)){
        which_cathode = 1;
      }
      else if(cut_Pb2[k]->IsInside(Anode2,Anode3)){
        which_cathode = 3;
      }  
      else if(cut_C[k]->IsInside(Anode2,Anode3)){
        which_cathode = 2;
      }
      else
        return;

      FissionFragmentAnalysis(which_cathode);
      //BeamFragmentAnalysis();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamFragmentAnalysis(){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();

  double L_CC = 8.45;
  TofPair TofHit;
  if(softofw_size==1){
    TofHit.plastic  = SofTofW->PlasticNbr[0];
    TofHit.x        = SofTofW->CalPosX[0];
    TofHit.y        = SofTofW->CalPosY[0];
    TofHit.tof      = SofTofW->CalTof[0];
    TofHit.velocity = L_CC/TofHit.tof;
    TofHit.beta     = TofHit.velocity * m/ns / NPUNITS::c_light;

    double Brho = 9.62543 + 0.0076642*TofHit.x;
    double Lfactor = 9.17/L_CC;
    double Beta = TofHit.beta*Lfactor;
    double Gamma1 = 1. / sqrt(1 - Beta * Beta);

    double AoQ = Brho / (3.10761 * Beta * Gamma1);

    SofFF->SetTOF(TofHit.tof);
    SofFF->SetTofPosX(TofHit.x);
    SofFF->SetTofPosY(TofHit.y);
    SofFF->SetPlastic(TofHit.plastic);

    SofFF->SetBeta(Beta);
    SofFF->SetGamma(Gamma1);
    SofFF->SetAoQ(AoQ);
    SofFF->SetBrho(Brho);


  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::FissionFragmentAnalysis(int which_cathode){
  unsigned int softofw_size = SofTofW->PlasticNbr.size();
  unsigned int softwim_size = SofTwim->SectionNbr.size();

  double E1 = -1;
  double E2 = -1;
  double E3 = -1;
  double E4 = -1;
  double DT1 = -1000;
  double DT2 = -1000;
  double DT3 = -1000;
  double DT4 = -1000;
  double Theta1 = -1000;
  double Theta2 = -1000;
  double Theta3 = -1000;
  double Theta4 = -1000;
  double L_CC = 8.45;
  double Beta_norm = 0.745;
  double Gamma1 = -1;
  double Gamma2 = -1;
  double Beta_Z1 = -1;
  double Beta_Z2 = -1;
  double Brho1 = -1;
  double Brho2 = -1;
  double AoQ1 = -1;
  double AoQ2 = -1;
  double A1 = -1;
  double A2 = -1;
  double Z1 = -1;
  double Z2 = -1;
  double Zsum = -1;
  int iZ1 = -1;
  int iZ2 = -1;
  int iZsum = -1;

  TofPair TofHit[2];
  if(softofw_size==2){
    for(unsigned int i=0; i<softofw_size; i++){
      TofHit[i].plastic  = SofTofW->PlasticNbr[i];
      TofHit[i].x        = SofTofW->CalPosX[i];
      TofHit[i].y        = SofTofW->CalPosY[i];
      TofHit[i].tof      = SofTofW->CalTof[i];
      TofHit[i].velocity = L_CC/TofHit[i].tof;
      TofHit[i].beta     = TofHit[i].velocity * m/ns / NPUNITS::c_light;
    }

    if(TofHit[0].x>TofHit[1].x){
      TofHit[0].isLorR = 1;
      TofHit[1].isLorR = 2;
    }
    else if(TofHit[0].x<TofHit[1].x){
      TofHit[0].isLorR = 2;
      TofHit[1].isLorR = 1;
    }

    if(TofHit[0].y>TofHit[1].y){
      TofHit[0].isUorD = 1;
      TofHit[1].isUorD = 2;
    }
    else if(TofHit[0].y<TofHit[1].y){
      TofHit[0].isUorD = 2;
      TofHit[1].isUorD = 1;
    }
  }

  vector<double> X1;
  vector<double> X2;
  vector<double> X3;
  vector<double> Y1;
  vector<double> Y2;
  vector<double> Y3;
  for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
    // *** MWPC1 *** //
    if(SofMwpc->DetectorNbr[i]==2){
      if(SofMwpc->PositionX1[i]!=-1000){
        X1.push_back(SofMwpc->PositionX1[i]);
      }
      if(SofMwpc->PositionX2[i]!=-1000){
        X1.push_back(SofMwpc->PositionX2[i]);
      }
      if(SofMwpc->PositionY[i]!=-1000){
        Y1.push_back(SofMwpc->PositionY[i]);
      }
    }
    // *** MWPC2 *** //
    if(SofMwpc->DetectorNbr[i]==3){
      if(SofMwpc->PositionX1[i]!=-1000){
        X2.push_back(SofMwpc->PositionX1[i]);
      }
      if(SofMwpc->PositionX2[i]!=-1000){
        X2.push_back(SofMwpc->PositionX2[i]);
      }
      if(SofMwpc->PositionY[i]!=-1000){
        Y2.push_back(SofMwpc->PositionY[i]);
      }
    }

    // *** MWPC3 *** //
    if(SofMwpc->DetectorNbr[i]==4){
      if(SofMwpc->PositionX1[i]!=-1000)
        X3.push_back(SofMwpc->PositionX1[i]);

      if(SofMwpc->PositionY[i]!=-1000)
        Y3.push_back(SofMwpc->PositionY[i]);
    }
  }


  for(unsigned int i=0; i<2; i++){
    double tofx = TofHit[i].x;

    for(unsigned int k=0; k<X3.size(); k++){
      double posx = X3[k];
      double posy = -1000;
      if(Y3.size()==X3.size())
        posy = Y3[k];
      if(abs(tofx-posx) < 150){
        if(abs(tofx-posx)<abs(tofx-TofHit[i].x3)){
          TofHit[i].x3 = posx;
          TofHit[i].y3 = posy;
        }
      }
    }
  }
  int ileft=0;
  int iright=0;
  if(TofHit[0].x3>TofHit[1].x3){
    ileft  = 0;
    iright = 1;
  }
  else{  
    ileft  = 1;
    iright = 0;
  }

  if(X1.size()==2){
    if(X1[0]>X1[1]){
      TofHit[ileft].x1 = X1[0];
      TofHit[iright].x1 = X1[1];
    }
    else if(X1[0]<X1[1]){
      TofHit[ileft].x1 = X1[1];
      TofHit[iright].x1 = X1[0];
    }
  }

  if(X2.size()==2){
    if(X2[0]>X2[1]){
      TofHit[ileft].x2 = X2[0];
      TofHit[iright].x2 = X2[1];
    }
    else if(X2[0]<X2[1]){
      TofHit[ileft].x2 = X2[1];
      TofHit[iright].x2 = X2[0];
    }
  }



  int mult1 = SofTwim->mult1;
  int mult2 = SofTwim->mult2;
  int mult3 = SofTwim->mult3;
  int mult4 = SofTwim->mult4;

  int multL = mult1 + mult2;
  int multR = mult3 + mult4;

  if(softwim_size>1){
    if( (mult1>1 && mult1<17) || (mult2>1 && mult2<17) || (mult3>1 && mult3<17) || (mult4>1 && mult4<17)){
      for(unsigned int i=0; i< softwim_size; i++){
        int sec = SofTwim->SectionNbr[i];

        if(sec==1){
          E1 = SofTwim->EnergySection[i];
          DT1 = SofTwim->DriftTime[i];
          Theta1 = SofTwim->Theta[i];
        }
        if(sec==2){
          E2 = SofTwim->EnergySection[i];
          DT2 = SofTwim->DriftTime[i];
          Theta2 = SofTwim->Theta[i];
        }
        if(sec==3){
          E3 = SofTwim->EnergySection[i]; 
          DT3 = SofTwim->DriftTime[i];
          Theta3 = SofTwim->Theta[i];
        }
        if(sec==4){
          E4 = SofTwim->EnergySection[i];     
          DT4 = SofTwim->DriftTime[i];
          Theta4 = SofTwim->Theta[i];
        }
      }

      if(softwim_size>2){
        if(E1>0 && E2>0){
          E1 = E1+E2;
          E2 = -1;
        }
        if(E3>0 && E4>0){
          E3 = E3+E4;
          E4 = -1;
        }
      }

      if(E1>0)
        E1 = E1/16;
      if(E2>0)
        E2 = E2/16;
      if(E3>0)
        E3 = E3/16;
      if(E4>0)
        E4 = E4/16;


      // *** case 1 *** //
      if(E1!=-1 && E2!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta1;

          TofHit[0].DT = DT2;
          TofHit[1].DT = DT1;

          TofHit[0].section = 2;
          TofHit[1].section = 1;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta2;

          TofHit[0].DT = DT1;
          TofHit[1].DT = DT2;

          TofHit[0].section = 1;
          TofHit[1].section = 2;
        }
      }

      // *** case 2 *** //
      if(E1!=-1 && E3!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta1;

          TofHit[0].DT = DT3;
          TofHit[1].DT = DT1;

          TofHit[0].section = 3;
          TofHit[1].section = 1;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E3;

          TofHit[0].DT = DT1;
          TofHit[1].DT = DT3;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta3;

          TofHit[0].section = 1;
          TofHit[1].section = 3; 
        }
      }

      // *** case 3 *** //
      if(E1!=-1 && E4!=-1){
        if(TofHit[0].isLorR==1 && TofHit[1].isLorR==2){
          TofHit[0].Esec = E1;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta1;
          TofHit[1].theta_in = Theta4;

          TofHit[0].DT = DT1;
          TofHit[1].DT = DT4;

          TofHit[0].section = 1;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isLorR==2 && TofHit[1].isLorR==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E1;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta1;

          TofHit[0].DT = DT4;
          TofHit[1].DT = DT1;

          TofHit[0].section = 4;
          TofHit[1].section = 1;
        }
      }

      // *** case 4 *** //
      if(E2!=-1 && E3!=-1){
        if(TofHit[0].isLorR==1 && TofHit[1].isLorR==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E3;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta3;

          TofHit[0].DT = DT2;
          TofHit[1].DT = DT3;

          TofHit[0].section = 2;
          TofHit[1].section = 3;
        }
        if(TofHit[0].isLorR==2 && TofHit[1].isLorR==1){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta2;

          TofHit[0].DT = DT3;
          TofHit[1].DT = DT2;

          TofHit[0].section = 3;
          TofHit[1].section = 2;
        }
      }

      // *** case 5 *** //
      if(E2!=-1 && E4!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E2;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta2;
          TofHit[1].theta_in = Theta4;

          TofHit[0].DT = DT2;
          TofHit[1].DT = DT4;

          TofHit[0].section = 2;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E2;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta2;

          TofHit[0].DT = DT4;
          TofHit[1].DT = DT2;

          TofHit[0].section = 4;
          TofHit[1].section = 2;
        }
      }

      // *** case 6 *** //
      if(E3!=-1 && E4!=-1){
        if(TofHit[0].isUorD==1 && TofHit[1].isUorD==2){
          TofHit[0].Esec = E3;
          TofHit[1].Esec = E4;

          TofHit[0].theta_in = Theta3;
          TofHit[1].theta_in = Theta4;

          TofHit[0].DT = DT3;
          TofHit[1].DT = DT4;

          TofHit[0].section = 3;
          TofHit[1].section = 4;
        }
        if(TofHit[0].isUorD==2 && TofHit[1].isUorD==1){
          TofHit[0].Esec = E4;
          TofHit[1].Esec = E3;

          TofHit[0].theta_in = Theta4;
          TofHit[1].theta_in = Theta3;

          TofHit[0].DT = DT4;
          TofHit[1].DT = DT3;

          TofHit[0].section = 4;
          TofHit[1].section = 3;
        }
      }


      // *** spline correction *** //
      for(int i=0; i<2; i++){
        int section = TofHit[i].section;

        double DT_eval;
        if(section<3) DT_eval=55;
        if(section>2) DT_eval=-55;
        if(section>0){
          TofHit[i].Esec = TofHit[i].Esec / fcorr_z_beta[section-1]->Eval(TofHit[i].beta) * fcorr_z_beta[section-1]->Eval(Beta_norm);

          TofHit[i].Esec = TofHit[i].Esec / fcorr_z_dt[section-1]->Eval(TofHit[i].DT) * fcorr_z_dt[section-1]->Eval(DT_eval);
        }
      }

  
      double ATToG = fDistanceStartToG - fDistancePlasticToCathode[which_cathode-1];
      double Theta0 = 20.*deg;//m_GladField->GetCentralTheta();
      double XA = 0;
      double ZA = fDistanceStartToA - fDistanceStartToG;
      double ZG = m_GladField->GetGladTurningPoint().Z();
      double ZMW3 = m_GladField->Get_MWPC3_Position().Z();
      double XMW3 = m_GladField->Get_MWPC3_Position().X();
      double X3lab = 0;
      double Z3lab = 0;;
      double Tilt = 14.*deg;//;
      TVector3 vZ = TVector3(0,0,1);
      TVector3 vStart = TVector3(0,0,-ATToG);
      double XC = 0;
      double ZC = 0;
      double XB = 0;
      double ZB = 0;
      double XD = 0;
      double ZD = 0;
      double MagB = m_GladField->GetB()/1000;
      for(int i=0; i<2; i++){
       XA = TofHit[i].DT;
        if(XA != -1e6 && SofBeamID->GetBeta()>0){
      
          TVector3 vG = TVector3(0,0,ZG);
          TVector3 vA = TVector3(XA,0,ZA);
          // *** Extroplate to C position *** //
          XC = (XA+(ZG-ZA)*tan(TofHit[i].theta_in)) / (1-tan(Tilt)*tan(TofHit[i].theta_in));
          ZC = ZG + XC*tan(Tilt);
          TVector3 vC    = TVector3(XC,0,ZC);
          TofHit[i].xc = XC;
          TofHit[i].yc = (ZC/8592.)*TofHit[i].y;
          TofHit[i].zc = ZC;

          int ix, iy;
          ix = (int)(TofHit[0].xc - m_GladField->GetXmin())/m_GladField->GetBin();
          iy = (int)(TofHit[0].yc - m_GladField->GetYmin())/m_GladField->GetBin();
          TofHit[i].Leff = m_GladField->GetLeff(ix,iy);

          X3lab = TofHit[i].x3*cos(Theta0) + XMW3;
          Z3lab = TofHit[i].x3*sin(Theta0) + ZMW3;
          TVector3 vE = TVector3(X3lab, TofHit[i].y, Z3lab);
          TVector3 dir = TVector3(sin(TofHit[i].theta_in), 0, cos(TofHit[i].theta_in));
          TofHit[i].x3lab = X3lab;
          TofHit[i].z3lab = Z3lab;
          TVector3 vOut  = TVector3(X3lab-XC,0,Z3lab-ZC);
          double angle = -vZ.Angle(vOut);
          TofHit[i].theta_out = angle;
          TofHit[i].psi = TofHit[i].theta_in - TofHit[i].theta_out;
          TofHit[i].rho = TofHit[i].Leff/(2*sin(0.5*TofHit[i].psi)*cos(Tilt-0.5*TofHit[i].psi));
          //TofHit[i].Brho = MagB*TofHit[i].rho;
          TofHit[i].Brho = m_GladField->FindBrho(vA,dir,vE);

          // *** Extrapolate to B position *** //
          double ZI = ZG - TofHit[i].Leff/(2*cos(Tilt));
          XB = (XA+(ZI-ZA)*tan(TofHit[i].theta_in)) / (1-tan(Tilt)*tan(TofHit[i].theta_in));
          ZB = ZI + XB*tan(Tilt);
          TofHit[i].xb = XB;
          TofHit[i].zb = ZB;
          TVector3 vB = TVector3(XB,0,ZB);

          // *** Extrapolate to D position *** //
          //XD = XC + TofHit[i].Leff/2*tan(angle+Tilt)*cos(Tilt);
          //ZD = ZC + TofHit[i].Leff/2*cos(Tilt);
          double phi = abs(angle) - Tilt;
          double psi = TMath::Pi()/2-abs(angle);
          double l = TofHit[i].Leff/(2*cos(phi));
          XD = XC - l*cos(psi);
          ZD = ZC + l*sin(psi);
          TofHit[i].xd = XD;
          TofHit[i].zd = ZD;
          TVector3 vD = TVector3(XD,0,ZD);

          // *** Extrapolate position of the rho point *** //
          TofHit[i].BrhoX = XB - TofHit[i].rho*cos(TofHit[i].theta_in); 
          TofHit[i].BrhoZ = ZB + TofHit[i].rho*sin(TofHit[i].theta_in); 
          TVector3 vBrho = TVector3(TofHit[i].BrhoX, 0, TofHit[i].BrhoZ);

          TVector3 v3lab = TVector3(X3lab,0,Z3lab);
          TVector3 v1 = TVector3(XB,0,ZB)-vStart;
          TVector3 v3 = TVector3(X3lab-XD,0,Z3lab-ZD);
          TofHit[i].omega = abs(2.*asin(sqrt(pow(XD-XB,2) + pow(ZD-ZB,2))/(2*TofHit[i].rho)));
          double Path1 = v1.Mag();
          double Path2 = TofHit[i].rho*TofHit[i].omega;
          double Path3 = v3.Mag();
          
          double delta_theta = TofHit[i].theta_out - m_GladField->GetCentralTheta();
          //double PathLength = Path1 + Path2 + Path3 + fDistanceMW3ToToF/cos(delta_theta);
          double PathLength = m_GladField->GetFlightPath(vStart, TofHit[i].Brho, vA, dir) + fDistanceMW3ToToF/cos(delta_theta);
          PathLength = PathLength/1000.;

          double BeamTimeOffset1 = 0;
          double BeamTimeOffset2 = 0;
          double BeamTimeOffset3 = 0;
          double BeamTimeOffset  = 0;
          double new_beta=0;

          // BeamTimeOffset between Start and first cathode //
          BeamTimeOffset1 = fDistancePlasticToCathode[0]/(SofBeamID->GetBeta()*NPUNITS::c_light);

          // BeamTimeOffset between first and second cathode //
          double par0= 3.442;
          double par1= -0.842;
          double par2= 1.070;
          new_beta =  par0*SofBeamID->GetBeta() + par1*exp(par2*SofBeamID->GetBeta());
          BeamTimeOffset2 = fDistanceBetweenCathode/(new_beta*NPUNITS::c_light);
          
          // BeamTimeOffset between second and third cathode //
          par0= 1.028;
          par1= -0.013;
          par2= 0.846;
          new_beta =  par0*new_beta + par1*exp(par2*new_beta);
          BeamTimeOffset3 = fDistanceBetweenCathode/(new_beta*NPUNITS::c_light);

          if(which_cathode==1)// 1st Lead
            BeamTimeOffset = BeamTimeOffset1;
          else if(which_cathode==2)// Carbon
            BeamTimeOffset = BeamTimeOffset1 + BeamTimeOffset2;
          else if(which_cathode==3)// 2nd Lead
            BeamTimeOffset = BeamTimeOffset1 + BeamTimeOffset2 + BeamTimeOffset3;
      
          TofHit[i].tof = TofHit[i].tof - BeamTimeOffset;
          TofHit[i].flight_path = PathLength;
          TofHit[i].velocity = PathLength/TofHit[i].tof;
          TofHit[i].beta     = TofHit[i].velocity * m/ns / NPUNITS::c_light;
          TofHit[i].gamma = 1. / sqrt(1 - pow(TofHit[i].beta,2));
          TofHit[i].AoQ = TofHit[i].Brho / (3.107 * TofHit[i].beta * TofHit[i].gamma);
          TofHit[i].x2twim = XA;

          double Z = TofHit[i].Esec;
          Z = fZff_p0 + fZff_p1*Z + fZff_p2*Z*Z;
          Z = sqrt(Z);
          int iZ = (int) round(Z);
          TofHit[i].Z = Z;
          TofHit[i].iZ = iZ;
          TofHit[i].A = TofHit[i].AoQ * TofHit[i].iZ;
          
          TVector3 vCG = vG - vC;
          TVector3 vCA = vA - vC;
          //TofHit[i].deff1 = vCG.Angle(vCA)*180./TMath::Pi();
          //TofHit[i].deff2 = vCG.Angle(vOut)*180./TMath::Pi();

          TVector3 vBD = vD-vB;
          double vBD_angle = vZ.Angle(vBD);

          double XO = XB - TofHit[i].Leff/(2*cos(vBD_angle - Tilt))*cos(TMath::Pi()/2-vBD_angle);
          double ZO = ZB + TofHit[i].Leff/(2*cos(vBD_angle - Tilt))*sin(TMath::Pi()/2-vBD_angle);
          TVector3 vO = TVector3(XO,0,ZO);
          TofHit[i].deff1 = (vO-vB).Mag();
          TofHit[i].deff2 = (vD-vO).Mag();

        }
      }

      // *** Filling the Fission Fragment Tree *** //
      for(int i=0; i<2; i++){
        SofFF->SetTOF(TofHit[i].tof);
        SofFF->SetTofPosX(TofHit[i].x);
        SofFF->SetTofPosY(TofHit[i].y);
        SofFF->SetPlastic(TofHit[i].plastic);
        SofFF->SetPosXB(TofHit[i].xb);
        SofFF->SetPosXC(TofHit[i].xc);
        SofFF->SetPosXD(TofHit[i].xd);
        SofFF->SetPosZB(TofHit[i].zb);
        SofFF->SetPosZC(TofHit[i].zc);
        SofFF->SetPosZD(TofHit[i].zd);
        SofFF->SetPosX1(TofHit[i].x1);
        SofFF->SetPosX2(TofHit[i].x2);
        SofFF->SetPosX3lab(TofHit[i].x3lab);
        SofFF->SetPosZ3lab(TofHit[i].z3lab);
        SofFF->SetThetaIn(TofHit[i].theta_in/deg);
        SofFF->SetThetaOut(TofHit[i].theta_out/deg);
        SofFF->SetBeta(TofHit[i].beta);
        SofFF->SetGamma(TofHit[i].gamma);
        SofFF->SetiZ(TofHit[i].iZ);
        SofFF->SetZ(TofHit[i].Z);
        SofFF->SetAoQ(TofHit[i].AoQ);
        SofFF->SetA(TofHit[i].A);
        SofFF->SetBrho(TofHit[i].Brho);
        SofFF->SetBrhoX(TofHit[i].BrhoX);
        SofFF->SetBrhoZ(TofHit[i].BrhoZ);
        SofFF->SetRho(TofHit[i].rho);
        SofFF->SetOmega(TofHit[i].omega);
        SofFF->SetDT(TofHit[i].DT);
        SofFF->SetSection(TofHit[i].section);
        SofFF->SetLeff(TofHit[i].Leff);
        SofFF->Setdeff1(TofHit[i].deff1);
        SofFF->Setdeff2(TofHit[i].deff2);
        SofFF->SetFlightPath(TofHit[i].flight_path);

      }
      SofFF->SetCathode(which_cathode);
      SofFF->SetZsum(TofHit[0].Z+TofHit[1].Z);
      SofFF->SetiZsum(TofHit[0].iZ+TofHit[1].iZ);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::BeamAnalysis(){
  unsigned int sofsci_size = SofSci->DetectorNbr.size();
  if(sofsci_size==2){
    double beta = SofSci->Beta[0];
    //cout << "Set beta to " << beta << endl;
    SofTrim->SetBeta(beta);
    SofTrim->BuildSimplePhysicalEvent();
    double Zbeam,Theta;
    if(SofTrim->EnergySection.size()>0){
      double Anode1 = SofTrim->EnergySection[0];
      double Anode2 = SofTrim->EnergySection[1];
      double Anode3 = SofTrim->EnergySection[2];
      if(fRunID==13)
        Zbeam = max(Anode1, Anode2);
      else 
        Zbeam = SofTrim->GetMaxEnergySection();

      Theta = SofTrim->Theta[0];

      double TofFromS2    = SofSci->CalTof[0];
      double velocity_mns = SofSci->VelocityMNs[0];
      double Beta         = SofSci->Beta[0];
      double XS2          = SofSci->PosMm[0];
      //double XCC          = SofSci->PosMm[1];
      double XCC=0;
      double YCC=0;
      for(unsigned int i=0; i<SofMwpc->DetectorNbr.size(); i++){
        if(SofMwpc->DetectorNbr[i]==1){
          XCC = SofMwpc->PositionX1[i];
          YCC = SofMwpc->PositionY[i];
        }
      }

      double LS2;
      LS2 = fLS2_0;//*(1 + fK_LS2*Theta);
      velocity_mns = LS2/TofFromS2;
      Beta = velocity_mns * m/ns / NPUNITS::c_light;
      double Gamma        = 1./(TMath::Sqrt(1 - TMath::Power(Beta,2)));
      double Brho = fBrho0 * (1 - XS2/fDS2 - XCC/fDCC);
      double AoQ  = Brho / (3.10716*Gamma*Beta);

      // Y dependence correction //
      double Y_p0 = 23943.8;
      double Y_p1 = 12.362;
      Zbeam = Zbeam/(Y_p0 + Y_p1*YCC)*Y_p0;


      // Z calibration //
      Zbeam = fZbeam_p0 + fZbeam_p1*Zbeam + fZbeam_p2*Zbeam*Zbeam;
      Zbeam = sqrt(Zbeam);
      double A = AoQ * round(Zbeam);

      // Last beta correction //
      double Beta_norm = 0.8355;
      Zbeam = Zbeam/(fZBeta_p0 + fZBeta_p1*Beta)*(fZBeta_p0 + fZBeta_p1*Beta_norm);

      // Filling Beam tree
      SofBeamID->SetZbeam(Zbeam);
      SofBeamID->SetAoQ(AoQ);
      SofBeamID->SetAbeam(A);
      SofBeamID->SetBeta(Beta);
      SofBeamID->SetGamma(Gamma);
      SofBeamID->SetBrho(Brho);
      SofBeamID->SetXS2(XS2);
      SofBeamID->SetXCC(XCC);
      SofBeamID->SetYCC(YCC);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadSpline(){
  TString input_path = "./calibration/SofTwim/spline/";

  TString rootfile = input_path + "spline_beta.root";
  TFile* ifile = new TFile(rootfile,"read");

  TString splinename;
  if(ifile->IsOpen()){
    cout << "**** Loading Beta spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_beta_sec%i",i+1);
      fcorr_z_beta[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

  //*** ***//
  rootfile = input_path + "spline_dt.root";
  ifile = new TFile(rootfile,"read");

  if(ifile->IsOpen()){
    cout << "**** Loading DT spline for fission fragment analysis..." << endl;
    for(int i=0; i<4; i++){
      splinename = Form("spline_dt_sec%i",i+1);
      fcorr_z_dt[i] = (TSpline3*) ifile->FindObjectAny(splinename);
    }
    ifile->Close();
  }
  else
    cout << "File " << rootfile << " not found!" << endl;

}


////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitParameter(){
  
  fK_LS2 = -30e-8;
  
  //fLS2_0 = 135.614;
  //fDS2   = 8000;
  //fDCC   = -10000;
  //fBrho0 = 12.3255;
  //fRunID = 12;

  // Beam parameter //
  fZBeta_p0 = 1;
  fZBeta_p1 = 0;
  fZbeam_p0 = 1651.57;
  fZbeam_p1 = 0.0876127;
  fZbeam_p2 = 4.02563e-6;

  // FF parameter //
  fZff_p0 = 2.80063;
  fZff_p1 = 6.91985e-2;
  fZff_p2 = 1.01598e-7;

  if(fRunID==1 || fRunID==2){
    //fBrho0 = 10.6813; // 180Hg
    fBrho0 = 10.6955; // 180Hg
    fZbeam_p0 = -5303.06;
    fZbeam_p1 = 0.674945;
    fZbeam_p2 = -8.32085e-6;

    fZBeta_p0 = 72.946;
    fZBeta_p1 = 6.0644;
  }
  if(fRunID==3){
    fBrho0 = 10.8183; // 182Hg
    fZbeam_p0 = -2737.25;
    fZbeam_p1 = 0.452017;
    fZbeam_p2 = -3.48831e-6;

    fZBeta_p0 = 76.6738;
    fZBeta_p1 = 1.60128;
  }
  if(fRunID==4){
    fBrho0 = 10.9558; // 184Hg
    fZbeam_p0 = -5044.61;
    fZbeam_p1 = 0.639986;
    fZbeam_p2 = -7.3077e-6;
  }
  if(fRunID==5){
    fBrho0 = 10.8138; // 187Pb
    fZbeam_p0 = -2858.72;
    fZbeam_p1 = 0.454064;
    fZbeam_p2 = -3.36443e-6;

    fZBeta_p0 = 71.0975;
    fZBeta_p1 = 10.7007;
  }
  if(fRunID==6){
    fBrho0 = 10.9476; // 189Pb
    fZbeam_p0 = 1590.66;
    fZbeam_p1 = 0.0956455;
    fZbeam_p2 = 3.84585e-6;

    fZBeta_p0 = 74.6063;
    fZBeta_p1 = 6.4635;
  }
  if(fRunID==7){
    fBrho0 = 10.6814; // 175Pt
    fZbeam_p0 = 459.68;
    fZbeam_p1 = 0.162277;
    fZbeam_p2 = 3.10164e-6;

    fZBeta_p0 = 66.9433;
    fZBeta_p1 = 10.8664;
  }
  if(fRunID==8){
    fBrho0 = 11.0864; // 204Fr
    fZbeam_p0 = 4122.94;
    fZbeam_p1 = -0.119867;
    fZbeam_p2 = 8.29115e-6;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==9){
    fBrho0 = 11.2712; // 207Fr
    fZbeam_p0 = -1752.27;
    fZbeam_p1 = 0.346018;
    fZbeam_p2 = -8.64673e-7;

    fZBeta_p0 = 63.9575;
    fZBeta_p1 = 25.1988;
  }
  if(fRunID==10){
    fBrho0 = 11.0955; // 199At run 423 & 424
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==11){
    fBrho0 = 10.9970; // 199At run 425 & 426
    fZbeam_p0 = -116.425;
    fZbeam_p1 = 0.218256;
    fZbeam_p2 = 1.62399e-6;

    fZBeta_p0 = 61.3889;
    fZBeta_p1 = 25.8908;
  }
  if(fRunID==12){
    fBrho0 = 10.8697; //197At
    fZbeam_p0 = -2683.52;
    fZbeam_p1 = 0.422551;
    fZbeam_p2 = -2.44471e-6;

    fZBeta_p0 = 62.9188;
    fZBeta_p1 = 22.8827;
  }
  if(fRunID==13){
    fBrho0 = 11.3418; // 216Th
    fZbeam_p0 = 1651.57;
    fZbeam_p1 = 0.0876127;
    fZbeam_p2 = 4.02563e-6;

    fZBeta_p0 = 38.879;
    fZBeta_p1 = 58.7667;
  }
  if(fRunID==14){
    fBrho0 = 11.5067; // 221Pa
    fZbeam_p0 = 186.892;
    fZbeam_p1 = 0.20739;
    fZbeam_p2 = 1.61797e-6;
  }
  if(fRunID==15){
    fBrho0 = 12.3352;
  }

}

////////////////////////////////////////////////////////////////////////////////
void Analysis::LoadActiveTargetCuts()
{
  TString element[14] = {"180Hg_1", "180Hg_2", "182Hg", "184Hg", "187Pb", "189Pb", "175Pt", "204Fr", "207Fr", "199At_1", "199At_2", "197At", "216Th", "221Pa"};

  TFile* file;
  TString filename1[14];
  TString filename2[14];
  TString filename3[14];
  cout << "/// Loading Active Target cuts..." << endl;
  for(int i=0; i<14; i++){
    filename1[i]= "./macro/cuts/"+element[i]+"/cut_Pb1.root";
    filename2[i]= "./macro/cuts/"+element[i]+"/cut_Pb2.root";
    filename3[i]= "./macro/cuts/"+element[i]+"/cut_C.root";

    file = new TFile(filename1[i],"read");
    cout << "- " << filename1[i] << endl;
    cut_Pb1[i] = (TCutG*) file->FindObjectAny("cut_Pb1");
    cut_Pb1[i]->SetName(element[i]+"_Pb1");

    file = new TFile(filename2[i],"read");
    cout << "- " << filename2[i] << endl;
    cut_Pb2[i] = (TCutG*) file->FindObjectAny("cut_Pb2");
    cut_Pb2[i]->SetName(element[i]+"_Pb2");

    file = new TFile(filename3[i],"read");
    cout << "- " << filename3[i] << endl;
    cut_C[i] = (TCutG*) file->FindObjectAny("cut_C");
    cut_C[i]->SetName(element[i]+"_C");
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReadAnalysisConfig()
{
  bool ReadingStatus = false;

  string filename = "AnalysisConfig.dat";
  
  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(filename.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No AnalysisConfig.dat found: Default parameter loaded for Analayis " << filename << endl;
    return;
  }
  cout << "**** Loading user parameter for Analysis from AnalysisConfig.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% AnalysisConfig.dat %%%");
  asciiConfig->Append(filename.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "AnalysisConfig";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="RunID") {
        AnalysisConfigFile >> DataBuffer;
        fRunID = atoi(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fRunID << endl;
      }

      else if (whatToDo=="Brho0") {
        AnalysisConfigFile >> DataBuffer;
        fBrho0 = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fBrho0 << endl;
      }

      else if (whatToDo=="LS2_0") {
        AnalysisConfigFile >> DataBuffer;
        fLS2_0 = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fLS2_0 << endl;
      }

      else if (whatToDo=="DispersionS2") {
        AnalysisConfigFile >> DataBuffer;
        fDS2 = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fLS2_0 << endl;
      }

      else if (whatToDo=="DispersionCC") {
        AnalysisConfigFile >> DataBuffer;
        fDCC = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDCC << endl;
      }    

      else if (whatToDo=="GladCurrent") {
        AnalysisConfigFile >> DataBuffer;
        fGladCurrent = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fGladCurrent << endl;
      }    

      else if (whatToDo=="DistanceStartToG") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceStartToG = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceStartToG << endl;
      }    

      else if (whatToDo=="DistanceStartToA") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceStartToA = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceStartToA << endl;
      }    

      else if (whatToDo=="DistanceBetweenATCathode") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceBetweenCathode = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceBetweenCathode << endl;
      }    

      else if (whatToDo=="DistanceMW3ToToF") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceMW3ToToF = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceMW3ToToF << endl;
      }    

      else if (whatToDo=="DistanceStartToFirstATCathode") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceStartToFirstATCathode = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceStartToFirstATCathode << endl;
      }    

      else if (whatToDo=="DistanceGToMW3") {
        AnalysisConfigFile >> DataBuffer;
        fDistanceGToMW3 = atof(DataBuffer.c_str());
        cout << "**** " << whatToDo << " " << fDistanceGToMW3 << endl;
      }    



      else {
        ReadingStatus = false;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  SofBeamID->Clear();
  SofFF->Clear();
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  RootOutput::getInstance()->GetTree()->Branch("RunID",&RunID,"RunID/I");
  RootOutput::getInstance()->GetTree()->Branch("SofBeamID","TSofBeamID",&SofBeamID);
  RootOutput::getInstance()->GetTree()->Branch("SofFissionFragment","TSofFissionFragment",&SofFF);

}
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy{
    public:
      proxy(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy p;
}


