#include "TReconstruction.h"

Int_t MRec[1100][550][3];

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;

ClassImp(TReconstruction);

//////////////////////////////////////////////////////////
TReconstruction::TReconstruction(){

}

//////////////////////////////////////////////////////////
TReconstruction::~TReconstruction(){

}

//////////////////////////////////////////////////////////
void TReconstruction::RandomInit()
{  
  Array = Ptr = new Float_t[255];
  for(UShort_t i=0;i<255;i++)
    *(Ptr++) = ((Float_t) i)/254.0 - 0.5;
  Ptr = Array;
}

//////////////////////////////////////////////////////////
float TReconstruction::GetRandom()
{
  float value = *Ptr;
 
  if(Ptr < Array+254)
    Ptr++;
  else
    Ptr = Array;
 
  return value;
}

//////////////////////////////////////////////////////////
void TReconstruction::ReadMatrix(Char_t * FName)
{
  int Len=255;
  int i,j;
  char Line[100];
  stringstream *InOut;
  ifstream File;

   ifstream IF;
   IF.open(FName);
   
   while(IF.getline(Line,Len))
     {
       InOut = new stringstream();
       *InOut << Line;
       *InOut >> i;
       *InOut >> j;
       *InOut >> MRec[i][j][0];
       *InOut >> MRec[i][j][1];
       *InOut >> MRec[i][j][2];
       delete InOut;
     } 
}

 
 
//////////////////////////////////////////////////////////
void Reconstruction::CalculateReconstruction(float Xf, float Tf, float BrhoRef, float&Brho, float&Theta, float&Path)
{
  Double_t Brhot,Thetat,Patht;
  Int_t VecI[2];
    
  Brho = Theta = Path = -500;                  

  Brhot  = 0.;
  Thetat = 0.;
  Patht  = 0.;

  VecI[0] = (int) (Xf+600.);
  VecI[1] = (int) (Tf+200.);
  if(
     (VecI[0] >=0 && VecI[0] < 1100)
     &&
     (VecI[1] >=0 && VecI[1] < 550)
     )
    {
      
      Brhot = (((Double_t) MRec[VecI[0]][VecI[1]][0]) + GetRandom())/1000.;
      Thetat = (((Double_t) MRec[VecI[0]][VecI[1]][1]) + GetRandom()) - 200.;
      Thetat *=-1;
      Patht = (((Double_t) MRec[VecI[0]][VecI[1]][2]) + GetRandom()) /10.;
    }

  if(Brhot >0.001 && Thetat > -300. && Thetat < 300.
                                                && Patht >0 && Patht < 2000.)
    {
      
      Brho = BrhoRef*((Float_t) Brhot);
      Theta = (Float_t) Thetat*-1;
      Path = (Float_t) Patht;
    }
}

