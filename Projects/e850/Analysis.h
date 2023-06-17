#ifndef Analysis_h 
#define Analysis_h
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
 *  This class describe  PISTA analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "NPVAnalysis.h"
#include "TPISTAPhysics.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "TInteractionCoordinates.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "TRandom3.h"
#include "TChain.h"

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void TreatEvent();
    void End();
    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();

    static NPL::VAnalysis* Construct();

  private:
    double BeamEnergy;
    double Xcalc;
    double Ycalc;
    double Zcalc;
    double XTarget;
    double YTarget;
    double ZTarget;
    double Elab;
    double DeltaE;
    double DeltaEcorr;
    double Eres;
    double ThetaLab;
    double PhiLab;
    double ThetaCM;
    double Ex240Pu;
    double Ex236U;
    double PID;
    double Chio_DE;
    double Chio_E;

    NPL::Reaction* Transfer10Be;
    NPL::Reaction* Transfer14C;

    TRandom3 Rand;
    double ThetaNormalTarget;
    double ThetaDetectorSurface;
    double TargetThickness;

    NPL::EnergyLoss C12C;
    NPL::EnergyLoss Be10C;
    NPL::EnergyLoss U238C;

  private:
    TPISTAPhysics* PISTA;
    TChain* chain;
   
    Float_t fIC[11];
    Float_t fTP_X;
    Float_t fTP_Y;
    ULong64_t fTS_TMW;
    Float_t fTAC_MW1_PISTA;
    Float_t fTAC_TMW1_FPMW0;
    Float_t fTAC_TMW1_FPMW1;
    Float_t fTAC_TMW2_FPMW0;
    Float_t fTAC_TMW2_FPMW1;

    // ECOGAM //
   Int_t           Inner6MVM;
   Float_t         Inner6MV[12];   //[Inner6MVM]
   UShort_t        Inner6MVN[12];   //[Inner6MVM]
   ULong64_t       Inner6MVTS[12];   //[Inner6MVM]
   Int_t           Inner20MVM;
   Float_t         Inner20MV[12];   //[Inner20MVM]
   UShort_t        Inner20MVN[12];   //[Inner20MVM]
   ULong64_t       Inner20MVTS[12];   //[Inner20MVM]
   Int_t           DeltaTVM;
   Float_t         DeltaTV[12];   //[DeltaTVM]
   UShort_t        DeltaTVN[12];   //[DeltaTVM]
   ULong64_t       DeltaTVTS[12];   //[DeltaTVM]
   Int_t           OutersVM;
   Float_t         OutersV[24];   //[OutersVM]
   UShort_t        OutersVN[24];   //[OutersVM]

};
#endif
