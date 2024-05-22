/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Vendeta analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"RootInput.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  InitOutputBranch();

  Vendeta= (TVendetaPhysics*) m_DetectorManager->GetDetector("Vendeta");
  FC= (TFissionChamberPhysics*) m_DetectorManager->GetDetector("FissionChamber");

  neutron = new NPL::Particle("1n");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

  ReInitValue();

  unsigned int FC_mult = FC->AnodeNumber.size();
  unsigned int HF_mult = FC->Time_HF.size();

  // Run 11
  //double GammaOffset[11] = {971.37, 970.67, 972.73, 972.59, 989.33, 982.03, 970.73, 985.13, 980.03, 976.83, 983.93};

  // Run 25++
  //double GammaOffset[11] = {969.58, 969.11, 975.03, 974.75, 978.90, 978.16, 967.34, 977.2, 981.14, 983.04, 982.32};

  // Run 59++
  //		double GammaOffset[11] = {969.28, 968.81, 987.43, 974.45, 994.30, 977.86, 980.25, 986.97, 996.30, 998.58, 997.12};

  //Run 87
  //		double GammaOffset[11] = {973.3, 972.81, 989.33, 978.45, 997.40, 982.01, 993.82, 995.10, 997.25, 999.25, 994.8};

  //Run 94++
  //  double GammaOffset[11] = {977.27, 976.72, 995.38, 984.38, 1007.49, 995.91, 999.38, 995.30, 1003.12, 999.25, 998.84};

  //Run 102++
  //double GammaOffset[11] = {975.14, 974.69, 985.36, 974.37, 1001.16, 987.81, 999.6, 997.08, 1001.04, 1003.2, 996.54};

  //Run 110++
 /* double GammaOffset[11] = {975.14, 974.69, 985.36, 974.37, 1001.16, 987.81, 999.6, 993.04, 1001.04, 1003.2, 996.54}; */
  
	//Run 120++
   /* double GammaOffset[11] = {975.298, 974.871, 985.779, 974.52, 996.487, 987.831, 1001.131, 993.013, 1001.013, 1001.964, 995.565}; */

  //Run 131-133
  double GammaOffset[11]={975.1, 974.97, 984.78, 974.7, 998.4, 987.8, 1002.6, 993.73, 1000.3, 1001.16, 995.17};
  double FakeFission_Offset = 987.94;

	//Run 135-137 (2.3 us pulse)
	/* double GammaOffset[11] = {1511.931, 1511.552, 1521.525, 1511.366, 1534.876 ,1524.462,1539.23, 1530.338, 1536.912, 1537.882, 1531.938}; */

	//Run 140++ 
  // double GammaOffset[11] = {975.1,975,988.6,978.7,998.3,983.7,1002.4,997.5,996.1,997.4,991.1};

  //Run 144 part 4 ++  (2.3 us pulse)
  /* double GammaOffset[11] = {1507.98, 1507.94, 1525.7, 1515.59, 1535.15,1520.7,1530.48, 1534.12, 1532.71, 1534.08, 1528.07}; */

	double incomingDT=0;
	double incomingE=0;
	double flight_path = 21500.;
	/*for(unsigned int j=0; j<HF_mult; j++){
		for(unsigned int i=0; i<FC_mult; i++){
		incomingDT = FC->Time[i] - FC->Time_HF[j] - GammaOffset[FC->AnodeNumber[i]-1];
		if(incomingDT<0){
		incomingDT += 1790;
		}
		double length = flight_path;// + 6*FC->AnodeNumber[i];	
		neutron->SetBeta((length/incomingDT) / c_light);
		incomingE = neutron->GetEnergy();

		inToF.push_back(incomingDT);
		inEnergy.push_back(incomingE);
		}
		}*/

	/*if(FC_mult==2){
		double HF1, HF2;
		for(int i=0; i<2; i++){
		HF1 = FC->Time[0] - FC->Time_HF[0];
		HF2 = FC->Time[1] - FC->Time_HF[1];
		}
		if(FC->AnodeNumber[0]>0 && FC->AnodeNumber[1]>0 && HF1<1790 && HF2<1790)
		}*/

	if(FC_mult>0){
		int EventMax=0;
		
		// double AnodeMaxQ=0;
		// for(unsigned int i=0;i<FC_mult;i++)
		//   {
		//     if(FC->AnodeNumber[i]==7) EventMax=i;
		// 	   //	if(AnodeMaxQ<AnodeQ) {EventMax=i; AnodeMaxQ=AnodeQ;}
    //   }

    int anode = FC->AnodeNumber[EventMax];
    double Time_FC = FC->Time[EventMax];
    bool isFake = FC->isFakeFission[EventMax];
    double Time_HF = FC->Time_HF[EventMax];
    incomingDT = FC->Time[EventMax] - FC->Time_HF[EventMax];
    TVector3 AnodePos = FC->GetVectorAnodePosition(6);
   
    if(anode>0){ 
      incomingDT -= GammaOffset[anode-1];
      AnodePos = FC->GetVectorAnodePosition(anode);
    }
    else if(anode ==0){
      incomingDT -= FakeFission_Offset;
    }
    /* if(FC_mult==3){ */
    /*   double inDT2 = FC->Time[1] -FC->Time_HF[1] - GammaOffset[FC->AnodeNumber[1]-1]; */
    /*   double inDT3 = FC->Time[2] -FC->Time_HF[2] - GammaOffset[FC->AnodeNumber[2]-1]; */
    /*   if( incomingDT - inDT2 > 0 ) */
    /* } */

    if(incomingDT<0){
      incomingDT += 1790;
      /* incomingDT += 2325; */
    }

    double length = flight_path;// + 6*FC->AnodeNumber[i];	
    neutron->SetBeta((length/incomingDT) / c_light);
    incomingE = neutron->GetEnergy();

    inToF.push_back(incomingDT);
    inEnergy.push_back(incomingE);

    FC_Q1.push_back(FC->Q1[EventMax]);
    FC_Q2.push_back(FC->Q2[EventMax]);
    FC_Qmax.push_back(FC->Qmax[EventMax]);
    FC_DT.push_back(FC->DT_FC[EventMax]);
    FC_FakeFission.push_back(FC->isFakeFission[EventMax]);
    FC_Anode_ID.push_back(anode);
    
    /* FC_Time.push_back(FC->Time[EventMax]); */
    /* HF_Time.push_back(FC->Time_HF[EventMax]); */

    Vendeta->SetAnodeNumber(anode);
    Vendeta->BuildPhysicalEvent();

    // VENDETA LG 
    unsigned int Vendeta_LG_mult = Vendeta->LG_DetectorNumber.size();
    for(unsigned int i=0; i<Vendeta_LG_mult; i++){

      int DetNbr          = Vendeta->LG_DetectorNumber[i];
      double Time_Vendeta = Vendeta->LG_Time[i];
      TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);
      double Rdet         = (DetPos-AnodePos).Mag() + 2 + 0.5*5.1 ; // Aluminum window + cell thickness
      double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendlg[DetNbr-1];
      
      if(DT>-500){
        double DeltaTheta = atan(63.5/Rdet);
        double Theta_Vendeta = DetPos.Theta();
        double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
        //neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
        //neutron->SetTimeOfFlight(DT*1e-9/(0.55));
        neutron->SetBeta(  (Rdet/DT) / c_light); 

        double En = neutron->GetEnergy();

        // Filling output tree
        LG_Tof.push_back(DT);
        LG_ID.push_back(DetNbr);
        LG_ELab.push_back(En);
        LG_ThetaLab.push_back(Theta_random);
        LG_Q1.push_back(Vendeta->LG_Q1[i]);
        LG_Q2.push_back(Vendeta->LG_Q2[i]);
        LG_Qmax.push_back(Vendeta->LG_Qmax[i]);
        /* LG_Time.push_back(Time_Vendeta); */

      }
    }

    // VENDETA HG 
    unsigned int Vendeta_HG_mult = Vendeta->HG_DetectorNumber.size();
    for(unsigned int i=0; i<Vendeta_HG_mult; i++){
      int DetNbr          = Vendeta->HG_DetectorNumber[i];
      double Time_Vendeta = Vendeta->HG_Time[i];
      TVector3 DetPos     = Vendeta->GetVectorDetectorPosition(DetNbr);
      double Rdet         = (DetPos-AnodePos).Mag() + 2 + 0.5*5.1 ; // Aluminum window + cell thickness
      double DT = Time_Vendeta - Time_FC;// + ToF_Shift_Vendhg[DetNbr-1];

      if(DT>-500){

        double DeltaTheta = atan(63.5/Rdet);

        double Theta_Vendeta = DetPos.Theta();
        double Theta_random = ra.Uniform(Theta_Vendeta-DeltaTheta,Theta_Vendeta+DeltaTheta);
        //neutron->SetTimeOfFlight(DT*1e-9/(Rdet*1e-3));
        //neutron->SetTimeOfFlight(DT*1e-9/(0.55));
        neutron->SetBeta( (Rdet/DT) / c_light); 
        double En = neutron->GetEnergy();
        // Filling output tree
        HG_ID.push_back(DetNbr);
        HG_Tof.push_back(DT);
        HG_ELab.push_back(En);
        HG_ThetaLab.push_back(Theta_random);
        HG_Q1.push_back(Vendeta->HG_Q1[i]);
        HG_Q2.push_back(Vendeta->HG_Q2[i]);
        HG_Qmax.push_back(Vendeta->HG_Qmax[i]);
        /* HG_Time.push_back(Time_Vendeta); */
      }
    }

    //Highlight saturated detectors
    static vector<int> LG_Saturated, HG_Saturated, LG_T_sat, HG_T_sat; 
    LG_Saturated.clear(), HG_Saturated.clear(), LG_T_sat.clear(), HG_T_sat.clear();
    for(int j = 0; j < LG_Tof.size();j++){
      int lgID = LG_ID[j];
      if( (lgID %4 < 2) && LG_Q1[j]>500e3){ // 50 Ohm
        LG_Saturated.push_back(lgID);
        LG_T_sat.push_back(LG_Tof[j]);
      }  
      else if( lgID %4 > 1 && LG_Q1[j]>590e3){ // 70 Ohm
        LG_Saturated.push_back(lgID);
        LG_T_sat.push_back(LG_Tof[j]);
      }
    }
    for(int j = 0; j < HG_Tof.size();j++){
      int hgID = HG_ID[j] ;
      if( (hgID %4 ==0 || hgID %4 ==1) && HG_Qmax[j] > 23800){
        HG_Saturated.push_back(hgID);
        HG_T_sat.push_back(HG_Tof[j]);
      }  
      else if( hgID %4 !=0 && hgID %4 !=1 && HG_Qmax[j] > 24300){
        HG_Saturated.push_back(hgID);
        HG_T_sat.push_back(HG_Tof[j]);
      }       
    }
    ///////////////////////////////////////////////////


    //Removes all signals from saturated detectors 
    for(int j = 0; j < LG_Saturated.size(); j++){
      for(int k = 0; k < HG_Tof.size(); k++){
        if(HG_ID[k] == LG_Saturated[j] && HG_Tof[k] > LG_T_sat[j]){
          HG_ID[k]       = -1;
          HG_Tof[k]      =-100000;
          HG_ELab[k]     =-100000;
          HG_ThetaLab[k] =-100000;
          HG_Q1[k]       =-100000;
          HG_Q2[k]       =-100000;
          HG_Qmax[k]     =-100000;
        }
      }
      for(int k = 0; k < LG_Tof.size(); k++){
        if(LG_ID[k] == LG_Saturated[j] && LG_Tof[k] > LG_T_sat[j]){
          LG_ID[k]       = -1;
          LG_Tof[k]      =-100000;
          LG_ELab[k]     =-100000;
          LG_ThetaLab[k] =-100000;
          LG_Q1[k]       =-100000;
          LG_Q2[k]       =-100000;
          LG_Qmax[k]     =-100000;
        }
      }
    }

    for(int j = 0; j < HG_Saturated.size(); j++){
      for(int k = 0; k < HG_Tof.size(); k++){
        if(HG_ID[k] == HG_Saturated[j] && HG_Tof[k] > HG_T_sat[j]){
          HG_ID[k]       = -1;
          HG_Tof[k]      =-100000;
          HG_ELab[k]     =-100000;
          HG_ThetaLab[k] =-100000;
          HG_Q1[k]       =-100000;
          HG_Q2[k]       =-100000;
          HG_Qmax[k]     =-100000;
        }
      }
    }

    ///////////////////////////////////////////////////

    //Process coincidences signals in VENDETA LG / HG
    /* if(HG_Tof.size() > 0 && LG_Tof.size() > 0 ){ */
    /*   for(int j = 0; j < LG_Tof.size();j++){ */
    /*     for(int k = 0; k < HG_Tof.size(); k++){ */
    /*       if(abs(HG_Tof[k]-LG_Tof[j]) < 3 && HG_ID[k] == LG_ID[j] && LG_ID[j]>0){ */
    /*         if( HG_Q1[k]>=120e3){ */
    /*           HG_ID[k]        = -1; */
    /*           HG_Tof[k]       = -100000; */
    /*           HG_ELab[k]      = -100000; */
    /*           HG_ThetaLab[k]  = -100000; */
    /*           HG_Q1[k]        = -100000; */
    /*           HG_Q2[k]        = -100000; */
    /*           HG_Qmax[k]      = -100000; */
    /*           /1* HG_Time[k] = - 100000; *1/ */
    /*         } */  
    /*         else if( HG_Q1[k]<120e3){ */
    /*           LG_ID[j]        = -1; */
    /*           LG_Tof[j]       = -100000; */
    /*           LG_ELab[j]      = -100000; */
    /*           LG_ThetaLab[j]  = -100000; */
    /*           LG_Q1[j]        = -100000; */
    /*           LG_Q2[j]        = -100000; */
    /*           LG_Qmax[j]      = -100000; */
    /*           /1* LG_Time[j] = -100000; *1/ */

    /*         } */
    /*       } */
    /*     } */
    /*   } */
    /* } // if LG && HG */

  }// if FC = 1
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  // Incoming neutron
  RootOutput::getInstance()->GetTree()->Branch("inToF",&inToF);
  RootOutput::getInstance()->GetTree()->Branch("inEnergy",&inEnergy);

  // FissionChamber
  RootOutput::getInstance()->GetTree()->Branch("FC_Q1",&FC_Q1);
  RootOutput::getInstance()->GetTree()->Branch("FC_Q2",&FC_Q2);
  RootOutput::getInstance()->GetTree()->Branch("FC_Qmax",&FC_Qmax);
  RootOutput::getInstance()->GetTree()->Branch("FC_DT",&FC_DT);
  RootOutput::getInstance()->GetTree()->Branch("FC_FakeFission",&FC_FakeFission);
  RootOutput::getInstance()->GetTree()->Branch("FC_Anode_ID",&FC_Anode_ID);
  /* RootOutput::getInstance()->GetTree()->Branch("FC_Time",&FC_Time); */
  /* RootOutput::getInstance()->GetTree()->Branch("HF_Time",&HF_Time); */

  // LG 
  RootOutput::getInstance()->GetTree()->Branch("LG_ID",&LG_ID);
  RootOutput::getInstance()->GetTree()->Branch("LG_ThetaLab",&LG_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("LG_ELab",&LG_ELab);
  RootOutput::getInstance()->GetTree()->Branch("LG_Tof",&LG_Tof);
  RootOutput::getInstance()->GetTree()->Branch("LG_Q1",&LG_Q1);
  RootOutput::getInstance()->GetTree()->Branch("LG_Q2",&LG_Q2);
  RootOutput::getInstance()->GetTree()->Branch("LG_Qmax",&LG_Qmax);
  /* RootOutput::getInstance()->GetTree()->Branch("LG_Time",&LG_Time); */

  // HG
  RootOutput::getInstance()->GetTree()->Branch("HG_ID",&HG_ID);
  RootOutput::getInstance()->GetTree()->Branch("HG_ThetaLab",&HG_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("HG_ELab",&HG_ELab);
  RootOutput::getInstance()->GetTree()->Branch("HG_Tof",&HG_Tof);
  RootOutput::getInstance()->GetTree()->Branch("HG_Q1",&HG_Q1);
  RootOutput::getInstance()->GetTree()->Branch("HG_Q2",&HG_Q2);
  RootOutput::getInstance()->GetTree()->Branch("HG_Qmax",&HG_Qmax);
  /* RootOutput::getInstance()->GetTree()->Branch("HG_Time",&HG_Time); */
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  inToF.clear();
  inEnergy.clear();

  LG_ThetaLab.clear();
  LG_ELab.clear();
  LG_Tof.clear();
  LG_ID.clear();
  LG_Q1.clear();
  LG_Q2.clear();
  LG_Qmax.clear();
  /* LG_Time.clear(); */

  HG_ThetaLab.clear();
  HG_ELab.clear();
  HG_Tof.clear();
  HG_ID.clear();
  HG_Q1.clear();
  HG_Q2.clear();
  HG_Qmax.clear();
  /* HG_Time.clear(); */

  FC_Q1.clear();
  FC_Q2.clear();
  FC_Qmax.clear();
  FC_DT.clear();
  FC_FakeFission.clear();
  FC_Anode_ID.clear();
  /* FC_Time.clear(); */
  /* HF_Time.clear(); */
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
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

