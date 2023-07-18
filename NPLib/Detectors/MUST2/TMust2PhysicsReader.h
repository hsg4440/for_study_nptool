#ifndef TMUST2PHYSICSREADER_H
#define TMUST2PHYSICSREADER_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    : July 2021
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold must2 treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <map>
#include <vector>
// NPL
#include "TMust2Data.h"
//#include "TMust2Physics.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TTree.h"
#include "NPVTreeReader.h"

// ROOT
//#include "TObject.h"

using namespace std;

// Forward Declaration

//class TMust2Physics;
//class TMust2PhysicsReader : public TMust2Physics  {
class TMust2PhysicsReader : public NPL::VTreeReader {
public:
  TMust2PhysicsReader();
  ~TMust2PhysicsReader(){};
  //TMust2PhysicsReader(TMust2Data* EventData){r_EventData = EventData;};
 // TMust2Data* GetData(){return r_EventData;}; 
  
public:
  void r_SetTreeReader(TTreeReader* TreeReader);
  //void r_BuildPhysicalEvent();
  //void r_ClearEventPhysics();
  //void r_ClearEventData();
  //void r_InitializeRootInputRaw();
  //void r_InitializeRootInputPhysics();  
  //void r_InitializeRootOutput();  

private:
  TTreeReader *dummy = new TTreeReader();
  //TTreeReaderValue<int>* r_Inner6MVM;
public:
  TTreeReaderValue<TMust2Data>* r_ReaderEventData = new TTreeReaderValue<TMust2Data>(*dummy,"");

public:
  ClassDef(TMust2PhysicsReader,0);

};


#endif
