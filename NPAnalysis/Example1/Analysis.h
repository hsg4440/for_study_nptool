// You can use this file to declare your spectra, file, energy loss , ... and whatever you want.
// This way you can remove all unnecessary declaration in the main programm.
// In order to help debugging and organizing we use Name Space.

/////////////////////////////////////////////////////////////////////////////////////////////////
// -------------------------------------- VARIOUS INCLUDE ---------------------------------------

// NPL
#include "DetectorManager.h"
#include "NPOptionManager.h"
#include "NPReaction.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TMust2Physics.h"
#include "TSSSDPhysics.h"
#include "TInitialConditions.h"
#include "NPEnergyLoss.h"
using namespace NPL ;
// STL C++
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
using namespace std;
// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TObject.h>

// ----------------------------------------------------------------------------------------------
void InitOutputBranch() ;
void InitInputBranch() ;
void ReInitValue() ;
/////////////////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------- DOUBLE, INT, BOOL AND MORE -------------------------------
namespace VARIABLE{
	double Ex;
  double ELab;
  double ThetaLab;
  TInitialConditions* Init = new TInitialConditions();
}
	 
using namespace VARIABLE ;
// ----------------------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------ENERGY LOSS----------------------------------------------
namespace ENERGYLOSS
	{
	
		//	Declare your Energy loss here	:
			EnergyLoss He3CD2 = EnergyLoss 	("He3_CD2.G4table","G4Table",1 );
													         
		   EnergyLoss He3Al = EnergyLoss 	("He3_Aluminium.G4table","G4Table",1);
		   
		   EnergyLoss He3Si = EnergyLoss 	("He3_Si.G4table","G4Table",1);
	
      EnergyLoss Li11CD2 = EnergyLoss 	("Li11[0.0]_CD2.G4table","G4Table",100);
}
	
using namespace ENERGYLOSS ;
// ----------------------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////////////////


