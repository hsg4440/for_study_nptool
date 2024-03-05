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
#include "TFPMWPhysics.h"
#include "TICPhysics.h"
#include "TVamosReconstruction.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "TInteractionCoordinates.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "TRandom3.h"
#include "TSpline.h"
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
    void LoadCalibParameter();

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
    int Telescope;
    double DeltaE;
    double DeltaEcorr;
    double Eres;
    double ThetaLab;
    double PhiLab;
    double ThetaCM;
    double Ex240Pu;
    double Ex236U;
    double Ex238U;
    double PID;

    double FF_Brho;
    double FF_Path;
    double FF_D;
    double FF_T;
    double FF_V;
    double FF_AoQ;
    double FF_Beta;
    double FF_Gamma;
    double FF_Q;
    double FF_M1;
    double FF_Mass;
    double FF_Etot;
    int FPMW_Section;
    double FF_D14;
    double FF_T14;
    double FF_V14;
    double FF_Q14;
    double FF_M114;
    double FF_AoQ14;
    double FF_Mass14;
    double FF_Qav;
    double FF_Massav;

    NPL::Reaction* Transfer10Be;
    NPL::Reaction* Transfer14C;
    NPL::Reaction* Elastic;

    TRandom3 Rand;
    double ThetaNormalTarget;
    double ThetaDetectorSurface;
    double TargetThickness;

    NPL::EnergyLoss C12C;
    NPL::EnergyLoss Be10C;
    NPL::EnergyLoss U238C;

    float T02;
    float T03;
    float T12;
    float T13;
    UShort_t FPMWPat_0RawNr[20];
    Int_t FPMWPat_0RawM;
    int Exo_Mult;
    vector<float>* Exo_Energy;
    vector<int>* Exo_Crystal;
    double Exogam_Energy;
    int Exogam_Crystal;


    double m_Q_p0[20];
    double m_Q_p1[20];
  private:
    TPISTAPhysics* PISTA;
    TFPMWPhysics* FPMW;
    TICPhysics* IC;
    TVamosReconstruction* Tracking;
    TChain* chain;
};
#endif
