#ifndef NPTIMESTAMP_h
#define NPTIMESTAMP_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Hugo Jacob hjacob@ijclab.in2p3.fr                       *
 *                                                                           *
 * Creation Date   : May 2024                                                *
 * Last update     : May 2024                                                *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Resources useful for TS treatment                                        *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>

// NPL
#include "NPInputParser.h"
using namespace NPL;

// ROOT header
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TVector3.h"

using namespace std;

namespace NPL {

  class TimeStamp {

   public: // Constructors and Destructors
    TimeStamp();
    ~TimeStamp(){};

   public: // Various Method
    void ReadConfigurationFile();
    // First argument is TS name, second one is TS value
    void AddTimeStamp(std::string TSname, unsigned long long TS){TSMap[TSname] = TS;};
    // Clears the map of TS
    void ClearTimeStamps(){ this->TSMap.clear();};
    //  
    void FindPrompts();
    //
    bool MatchTS(map<std::string,unsigned long long>::iterator TS1,
      map<std::string,unsigned long long>::iterator TS2);
    bool MatchTS(std::string TS1, std::string TS2);

    std::vector<std::string> GetGoodTS(){return GoodTS;};
   private:
    // Parameters to initialize by reading the config
    std::map<std::string, unsigned long long> TSConditions;
    unsigned long long TSWidth;
    
    // Map of all timestamps added
    std::map<std::string, unsigned long long> TSMap;
    // Initialization of 1 TS that will be the reference for all other TS
    std::string TSReference;
    // Vector of all the good TS after FindPrompts()
    std::vector<std::string> GoodTS;
    
    int fVerboseLevel;

    ClassDef(TimeStamp, 0)
  };
} // namespace NPL
#endif
