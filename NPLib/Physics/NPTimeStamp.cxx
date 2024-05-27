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


#include "NPTimeStamp.h"

#include <fstream>
#include <iostream>


// ROOT
#include "TF1.h"

ClassImp(TimeStamp)

TimeStamp::TimeStamp(){

};

// Checking if the TS diff is in the prompt within a certain width
bool TimeStamp::MatchTS(map<std::string,unsigned long long>::iterator TS1,
  map<std::string,unsigned long long>::iterator TS2){
    return abs((long long)((*TS1).second -(*TS2).second - (TSConditions[(*TS1).first] - TSConditions[(*TS1).first]))) < TSWidth;
};

// Checking if the TS diff is in the prompt within a certain width
bool TimeStamp::MatchTS(std::string TS1, std::string TS2){
    return abs((long long)(TSMap[TS1]-TSMap[TS2] - (TSConditions[TS1] - TSConditions[TS2]))) < TSWidth;
};

void TimeStamp::FindPrompts(){
  auto TSReference_it = find(TSMap.begin(), TSMap.end(), TSReference); 
  // Checking that the TS reference exists, we assume that at least 1 TS is the ref and is absolutely necessary
  // in the TS treatment
  if(TSReference_it != TSMap.end()){
    // CLear the good TS vector and adding the reference TS
    GoodTS.clear();
    GoodTS.push_back(TSReference_it->first);
    for(auto TS_it = TSMap.begin(); TS_it != TSMap.end(); TS_it++){
      // If prompt with the reference TS, add the TS to the good TS vector
      // else do nothing
      if(MatchTS(TSReference_it,TS_it))
        GoodTS.push_back(TS_it->first);
    }
  }
  else{
    std::cout << "TS matching not done : TS reference not correctly found" << std::endl;
  }
}
void TimeStamp::ReadConfigurationFile() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/TSConfig.dat";

  // open analysis config file
  ifstream TSConfigFile;
  TSConfigFile.open(FileName.c_str());

  if (!TSConfigFile.is_open()) {
    cout << " No TSConfig.dat found: No parameters loaded"  << endl;
    return;
  }
  cout << " Loading user parameter for TS from TSConfig.dat " << endl;

  // read config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!TSConfigFile.eof()) {
    // Pick-up next line
    getline(TSConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigZDD";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      TSConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        TSConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="ADD_TS") {
        TSConfigFile >> DataBuffer;
        std::string Temp = DataBuffer;
        TSConfigFile >> DataBuffer;
        TSMap[Temp] = atoll(DataBuffer.c_str());
        cout << whatToDo << " " << Temp << " with prompt value: " << DataBuffer << endl;
      }
      else if (whatToDo=="TS_WIDTH") {
        TSConfigFile >> DataBuffer;
        TSWidth = atoll(DataBuffer.c_str());
        cout << whatToDo << " " <<  DataBuffer << endl;
      }
      else if (whatToDo=="TS_REF") {
        TSConfigFile >> DataBuffer;
        TSReference = DataBuffer;
        cout << whatToDo << " " <<  DataBuffer << endl;
      }
      else {
        ReadingStatus = false;
      }
    }
  }
}