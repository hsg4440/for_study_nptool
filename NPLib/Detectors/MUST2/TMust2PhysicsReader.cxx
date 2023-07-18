/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
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
 *  This class hold must2 TreeReader                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TMust2PhysicsReader.h"
#include "NPDetectorFactory.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TMust2PhysicsReader)

TMust2PhysicsReader::TMust2PhysicsReader()
{
};

void TMust2PhysicsReader::r_SetTreeReader(TTreeReader* TreeReader){
std::cout << "test TreeReader " << r_ReaderEventData  << std::endl;
r_ReaderEventData = new TTreeReaderValue<TMust2Data>(*TreeReader,"MUST2");
std::cout << "test TreeReader 2 " << r_ReaderEventData << std::endl;
//r_Inner6MVM = new TTreeReaderValue<int>(*TreeReader,"Inner6MVM");
}; 

//void TMust2PhysicsReader::r_BuildPhysicalEvent(){
//    std::cout << "Test 1 " << GetRawData() << std::endl;
//    SetRawDataPointer(&(**r_ReaderEventData));
//    std::cout << "Test 2 " << GetRawData() << std::endl;
//    (&(**r_ReaderEventData))->Dump();
//    BuildPhysicalEvent();
//};
//
//void TMust2PhysicsReader::r_ClearEventPhysics(){
//    TMust2Physics::ClearEventPhysics();
//};
//
//void TMust2PhysicsReader::r_ClearEventData(){
//    TMust2Physics::ClearEventData();
//};

//void TMust2PhysicsReader::r_InitializeRootInputPhysics(){
//  std::cout << "*************************** TEST HUGO 2 r ***************************************" << std::endl; 
//  TChain* inputChain = RootInput::getInstance()->GetChain();
//  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
//  inputTreeReader = new TTreeReader(inputChain);
//};
//
//void TMust2PhysicsReader::r_InitializeRootInputRaw(){
//  std::cout << "*************************** TEST HUGO r ***************************************" << std::endl; 
//  TChain* inputChain = RootInput::getInstance()->GetChain();
//  inputChain->Print();
//  std::cout << "adresse chain " << inputChain << std::endl;
//  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
//  inputTreeReader->SetTree(inputChain);
//  std::cout << "Test initialize root input 1" << std::endl;
//};

//void TMust2PhysicsReader::r_InitializeRootOutput() {
//  TMust2Physics::InitializeRootOutput();
//}