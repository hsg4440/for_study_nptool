#ifndef TSAMURAIFDC0PHYSICS_H
#define TSAMURAIFDC0PHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2020                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SamuraiFDC0 treated data                                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
#include <map>

// NPL
#include "TSamuraiFDC0Data.h"
#include "SamuraiDCIndex.h"
//#include "TSamuraiFDC0Spectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"
#include "NPDCReconstruction.h"
// ROOT 
#include "TVector3.h" 

// Forward declaration
//class TSamuraiFDC0Spectra;


using namespace std ;

class TSamuraiFDC0Physics : public TObject, public NPL::VDetector{
  public:
    TSamuraiFDC0Physics();
    ~TSamuraiFDC0Physics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    //   Provide Physical Multiplicity
    vector<double> DriftLength;
    vector<int> Detector;
    vector<int> Layer;
    vector<int> Wire;
    vector<double> Time;
    vector<double> ToT;
    vector<bool>   Matched;
    // Computed variable
    vector<TVector3> ParticleDirection;
    vector<TVector3> MiddlePosition;

    double PosX;
    double PosY;
    double ThetaX;
    double PhiY;
    double devX,devY;
    TVector3 Dir;
    int Mult;
    int MultMean;
    int PileUp;

  public:
    // Projected position at given Z plan
    TVector3 ProjectedPosition(double Z);

  private: // Charateristic of the DC 
    void AddDC(string name, NPL::XmlParser&);//! take the XML file and fill in Wire_X and Layer_Angle
    map<SamuraiDCIndex,double> Wire_X;//! X position of the wires
    map<SamuraiDCIndex,double> Wire_Z;//! Z position of the wires
    map<SamuraiDCIndex,double> Wire_Angle;//! Wire Angle (0 for X, 90 for Y, U and V are typically at +/-30)
  
  private: // Analysis
    double ToTThreshold_H;//! a ToT Low threshold to remove noise
    double ToTThreshold_L;//! a ToT High threshold to remove noise
    // since the calibration is a sigmoid there quite a few event at the edge 
    double DriftLowThreshold;//! Minimum Drift length to keep the hit 
    double DriftUpThreshold;//! Maximum Drift length to keep the hit
    double PowerThreshold;//! Maximum P2 minimisation value to keep the track   
    void RemoveNoise();
    // Construct the 2D track and ref position at Z=0 and Z=100 based on X,Z and Radius provided

    // Object use to perform the DC reconstruction
    NPL::DCReconstruction m_reconstruction;//!

  public: //   Innherited from VDetector Class

    // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;


    // Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() ;      

    // Activated associated Branches and link it to the private member DetectorData address
    // In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    // Activated associated Branches and link it to the private member DetectorPhysics address
    // In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    // Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent() ;

    // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    // This method aimed to be used for analysis performed during experiment, when speed is requiered.
    // NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() ;

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();
    // Write Spectra to file
    void WriteSpectra();

  public:      //   Specific to SamuraiFDC0 Array

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    // Retrieve raw and pre-treated data
    TSamuraiFDC0Data* GetRawData()        const {return m_EventData;}
    TSamuraiFDC0Data* GetPreTreatedData() const {return m_PreTreatedData;}
  
    double GetPosX(){return PosX;}
    double GetPosY(){return PosY;}
    double GetThetaX(){return ThetaX;}
    double GetDevX(){return devX;}
    double GetDevY(){return devY;}
    int GetPileUp(){return PileUp;}

  private:   //   Root Input and Output tree classes
    TSamuraiFDC0Data*         m_EventData;//!
    TSamuraiFDC0Data*         m_PreTreatedData;//!
    TSamuraiFDC0Physics*      m_EventPhysics;//!


  private: // Spectra Class
   // TSamuraiFDC0Spectra* m_Spectra; // !

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TSamuraiFDC0Physics,1)  // SamuraiFDC0Physics structure
};

#endif
