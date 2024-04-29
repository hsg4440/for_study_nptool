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
    void TwoAlphaAnalysis();
    void ReadAnalysisConfig();

    static NPL::VAnalysis* Construct();

  private:
    double m_XTarget_offset;
    double m_YTarget_offset;
    double m_ZTarget_offset;
    double m_Beam_ThetaX;
    double m_Beam_ThetaY;
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
    int strip_DE;
    int strip_E;
    double Time_E;
    double DeltaE;
    double DeltaEcorr;
    double Eres;
    double ThetaLab;
    double PhiLab;
    double ThetaCM;
    double Ex242Pu;
    double Ex240Pu;
    double Ex236U;
    double Ex238U;
    double PID;

    ULong64_t fVAMOS_TS_sec;
    ULong64_t fPISTA_TS_sec;
    double VAMOS_TS_hour;
    double PISTA_TS_hour;

    int FPMW_Section;
    double FF_Brho;
    double FF_Path;
    
    double FF_D13;
    double FF_T13;
    double FF_V13;
    double FF_AoQ13;
    double FF_Beta13;
    double FF_Gamma13;
    double FF_Q13;
    double FF_M113;
    double FF_Mass13;
    double FF_Etot13;
    
    double FF_D14;
    double FF_T14;
    double FF_V14;
    double FF_Q14;
    double FF_M114;
    double FF_AoQ14;
    double FF_Mass14;
    
    double FF_D23;
    double FF_T23;
    double FF_V23;
    double FF_Q23;
    double FF_M123;
    double FF_AoQ23;
    double FF_Mass23;

    double FF_D24;
    double FF_T24;
    double FF_V24;
    double FF_Q24;
    double FF_M124;
    double FF_AoQ24;
    double FF_Mass24;
  
    double FF_Qav;
    double FF_Massav;

    int m_2alpha;
    vector<double> Elab1;
    vector<double> Elab2;

    NPL::Reaction* Transfer10Be;
    NPL::Reaction* Transfer8Be;
    NPL::Reaction* Transfer14C;
    NPL::Reaction* Elastic;
    
    TVector3 PositionOnTarget;
    TRandom3 Rand;
    double ThetaNormalTarget;
    double ThetaDetectorSurface;
    double TargetThickness;

    NPL::EnergyLoss C12C;
    NPL::EnergyLoss Be10C;
    NPL::EnergyLoss U238C;

    float T13;
    float T14;
    float T23;
    float T24;
    UShort_t FPMWPat_0RawNr[20];
    Int_t FPMWPat_0RawM;
    int Exo_Mult;
    vector<float>* Exo_Energy;
    vector<int>* Exo_Crystal;
    double Exogam_Energy;
    int Exogam_Crystal;

  private:
    double Xmean;
    double Ymean;
    double Xmean_iter;
    double Ymean_iter;
    int iteration;

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
