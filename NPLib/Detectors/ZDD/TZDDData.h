#ifndef __ZDDDATA__
#define __ZDDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <iostream>
using namespace std;

// ROOT
#include "TObject.h"


class TZDDData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    double counter = 0; 
    // Energy IC
    vector<UShort_t>   fZDD_E_IC_Nbr;
    vector<Double_t>   fZDD_IC_Energy;

    // Time IC
    vector<UShort_t>   fZDD_T_IC_Nbr;
    vector<Double_t>   fZDD_IC_Time;
    
    // Energy Plastic
    vector<UShort_t>   fZDD_E_Plastic_Nbr;
    vector<Double_t>   fZDD_Plastic_Energy;
    
    // Time Plastic
    vector<UShort_t>   fZDD_T_Plastic_Nbr;
    vector<Double_t>   fZDD_Plastic_Time;

    // Energy
    vector<UShort_t>   fZDD_E_DetectorNbr;
    vector<Double_t>   fZDD_Energy;

    // Time
    vector<UShort_t>   fZDD_T_DetectorNbr;
    vector<Double_t>   fZDD_Time;

    // DriftTime in DC
    vector<UShort_t> fZDD_Drift_DetectorNbr;
    vector<Double_t> fZDD_DriftTime;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TZDDData();
    ~TZDDData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // Energy
    inline void SetEnergy(const UShort_t& DetNbr,const Double_t& Energy){
      fZDD_E_DetectorNbr.push_back(DetNbr);
      fZDD_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetTime(const UShort_t& DetNbr,const Double_t& Time)	{
      fZDD_T_DetectorNbr.push_back(DetNbr);     
      fZDD_Time.push_back(Time);
    };//!

    // IC Energy
    inline void Set_IC_Energy(const UShort_t& IC_Nbr,const Double_t& IC_Energy){
      fZDD_E_IC_Nbr.push_back(IC_Nbr);
      fZDD_IC_Energy.push_back(IC_Energy);
    };//!
    
    // IC Time
    inline void Set_IC_Time(const UShort_t& IC_Nbr,const Double_t& IC_Time)	{
      fZDD_T_IC_Nbr.push_back(IC_Nbr);     
      fZDD_IC_Time.push_back(IC_Time);
    };//!

    // Plastic Energy
    inline void Set_Plastic_Energy(const UShort_t& Plastic_Nbr,const Double_t& Plastic_Energy){
      fZDD_E_Plastic_Nbr.push_back(Plastic_Nbr);
      fZDD_Plastic_Energy.push_back(Plastic_Energy);
    };//!
    
    // Plastic Time
    inline void Set_Plastic_Time(const UShort_t& Plastic_Nbr,const Double_t& Plastic_Time)	{
      fZDD_T_Plastic_Nbr.push_back(Plastic_Nbr);     
      fZDD_Plastic_Time.push_back(Plastic_Time);
    };//!
    
    // Position DriftTime and X in DC
    inline void Set_DC_Time(const UShort_t& DetNbr, const Double_t& DriftTime) {
        fZDD_Drift_DetectorNbr.push_back(DetNbr);
        fZDD_DriftTime.push_back(DriftTime);
    }; //!
    
    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy(std::string Detector) const
      {if(Detector == "IC")
        return fZDD_E_IC_Nbr.size();
      else if(Detector == "Plastic")
        return fZDD_E_Plastic_Nbr.size();
      else{
        std::cout << "Detector should be either IC or Plastic" << std::endl;
        return -1;
      }
      }
    
    inline UShort_t GetE_DetectorNbr(std::string Detector, const unsigned int &i) const
      {if(Detector == "IC")
        return fZDD_E_IC_Nbr[i];
      else if(Detector == "Plastic")
        return fZDD_E_Plastic_Nbr[i];
      else{
        std::cout << "Detector should be either IC or Plastic" << std::endl;
        return -1;
      }
      }
    inline Double_t Get_Energy(std::string Detector, const unsigned int &i) const
      {if(Detector == "IC")
        return fZDD_IC_Energy[i];
      else if(Detector == "Plastic")
        return fZDD_Plastic_Energy[i];
      else{
        std::cout << "Detector should be either IC or Plastic" << std::endl;
        return -1;
      }
      }
    
    // Time
    inline UShort_t GetMultTime(std::string Detector) const
      {if(Detector == "IC")
        return fZDD_T_IC_Nbr.size();
      else if(Detector == "Plastic")
        return fZDD_T_Plastic_Nbr.size();
      else if(Detector == "DC")
        return fZDD_Drift_DetectorNbr.size();
      else{
        std::cout << "Detector should be either IC, DC or Plastic" << std::endl;
        return -1;
      }
      }
    
    inline UShort_t GetT_DetectorNbr(std::string Detector, const unsigned int &i) const
      {if(Detector == "IC")
        return fZDD_T_IC_Nbr[i];
      else if(Detector == "Plastic")
        return fZDD_T_Plastic_Nbr[i];
      else if(Detector == "DC")
        return fZDD_Drift_DetectorNbr[i];
      else{
        std::cout << "Detector should be either IC, DC or Plastic" << std::endl;
        return -1;
      }
      }
    inline Double_t Get_Time(std::string Detector, const unsigned int &i) const
      {if(Detector == "IC")
        return fZDD_IC_Time[i];
      else if(Detector == "Plastic")
        return fZDD_Plastic_Time[i];
      else if(Detector == "DC")
        return fZDD_DriftTime[i];
      else{
        std::cout << "Detector should be either IC, DC or Plastic" << std::endl;
        return -1;
      }
      }
    // Position
    inline UShort_t GetMultDrift() const { return fZDD_Drift_DetectorNbr.size(); }
    inline UShort_t GetDrift_DetectorNbr(const unsigned int& i) const {
        return fZDD_Drift_DetectorNbr[i];
    } //!
    inline Double_t Get_DriftTime(const unsigned int& i) const {
        return fZDD_DriftTime[i];
    } //!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TZDDData,1)  // ZDDData structure
};

#endif