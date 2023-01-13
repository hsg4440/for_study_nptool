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
 *  This class describe  Sofia analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "NPVAnalysis.h"
#include "TInitialConditions.h"
#include "TInteractionCoordinates.h"
#include "TFissionConditions.h"
#include "GladFieldMap.h"
#include "TRandom3.h"

class Analysis: public NPL::VAnalysis{
  public:
    Analysis();
    ~Analysis();

  public: 
    void Init();
    void InitInputBranch();
    void InitOutputBranch();
    void TreatEvent();
    void ReInit();
    void End();

    static NPL::VAnalysis* Construct();

  private:
    vector<double> m_TOF;
    vector<double> m_Brho;
    vector<int> m_A;
    vector<double> m_Brho_calc;
    vector<int> m_A_calc;
    vector<double> m_FlightPath;

  private:
    TInitialConditions* InitialConditions;
    TInteractionCoordinates* InteractionCoordinates;
    TFissionConditions* FissionConditions;

    TRandom3 ran;
    GladFieldMap* m_GladField;
};
#endif
