#ifndef TFPMWPHYSICS_H
#define TFPMWPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace  contact address: pierre.morfouace@cea.fr   *
 *                                                                           *
 * Creation Date  : October 2023                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FPMW Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <set>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
#include "TSpline.h"
// NPTool headers
#include "TFPMWData.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"



class TFPMWPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TFPMWPhysics();
    ~TFPMWPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      DetectorNbr;
    vector<double>   PositionX;
    vector<double>   PositionY;
    vector<double>   ChargeX;
    vector<double>   ChargeY;
  
  private:
    std::set<int> DetectorHitX;//!
    std::set<int> DetectorHit;//!
    vector<double>   DetPosX;//!
    vector<double>   DetPosY;//!
    vector<double>   DetPosZ;//!
    vector<vector<double>> Buffer_X_Q;//!
    vector<vector<double>> Buffer_Y_Q;//!
    
    map<int, vector<pair<int,double>>> MapX;//!
    map<int, vector<pair<int,double>>> MapY;//!
    map<int, pair<int, double>> MaxQX;//!
    map<int, pair<int, double>> MaxQY;//!
    map<int, double> QSumX;//!
    map<int, double> QSumY;//!


    /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 Pos); 
  void AddDetector(double R, double Theta, double Phi); 
  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    // add parameters to the CalibrationManger
    void AddParameterToCalibrationManager();

    // method called event by event, aiming at extracting the 
    // physical information from detector
    void BuildPhysicalEvent();

    // same as BuildPhysicalEvent() method but with a simpler
    // treatment
    void BuildSimplePhysicalEvent();

    // same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // activate raw data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputRaw();

    // activate physics data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputPhysics();

    // create branches of output ROOT file
    void InitializeRootOutput();

    // clear the raw and physical data objects event by event
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   


  //////////////////////////////////////////////////////////////
  // specific methods to FPMW array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TFPMWData object to TFPMWPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TFPMWData* rawDataPointer) {m_EventData = rawDataPointer;}

    double fThresholdX;
    double fThresholdY;

    double AnalyticHyperbolicSecant(std::pair<int,double>& MaxQ, std::vector<std::pair<int, double>>& Map);


  // objects are not written in the TTree
  private:
    TFPMWData*         m_EventData;        //!
    TFPMWData*         m_PreTreatedData;   //!
    TFPMWPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TFPMWData* GetRawData()        const {return m_EventData;}
    TFPMWData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    double m_E_Threshold;     //!

  // number of detectors
  private:
    int m_NumberOfDetectors;  //!

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TFPMWPhysics,1)  // FPMWPhysics structure
};
#endif
