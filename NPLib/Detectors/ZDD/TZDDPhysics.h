#ifndef TZDDPHYSICS_H
#define TZDDPHYSICS_H
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
 *  This class hold ZDD Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
// NPTool headers
#include "TZDDData.h"
#include "TZDDSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TZDDSpectra;



class TZDDPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TZDDPhysics();
    ~TZDDPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      IC_DetectorNumber;
    vector<double>   IC_Energy;
    vector<double>   IC_Time;

    vector<int>      Plastic_DetectorNumber;
    vector<double>   Plastic_Energy;
    vector<double>   Plastic_Time;

    vector<int>      DC_DetectorNumber;
    vector<double>   DC_DriftTime;
  
  /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 POS, string shape); 
  void AddDetector(double R, double Theta, double Phi, string shape); 
  
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

    // methods related to the TZDDSpectra class
    // instantiate the TZDDSpectra class and 
    // declare list of histograms
    void InitSpectra();

    // fill the spectra
    void FillSpectra();

    // used for Online mainly, sanity check for histograms and 
    // change their color if issues are found, for example
    void CheckSpectra();

    // used for Online only, clear all the spectra
    void ClearSpectra();

    // write spectra to ROOT output file
    void WriteSpectra();
    
    void Add_ZDD(double R, double Theta) {
  /*    cout << "Balise 2" << endl;
        m_R     = R;
        m_Theta = Theta;
  */
  m_NumberOfDetectors++;
    }

    void Add_Ionisation_Chamber(double Z, double Thickness, string Gas, double Pressure,
            double Temperature) {
        m_Ionisation_Chamber_Z.push_back(Z);
        m_Ionisation_Chamber_Thickness.push_back(Thickness);
        m_Ionisation_Chamber_Gas.push_back(Gas);
        m_Ionisation_Chamber_Pressure.push_back(Pressure);
        m_Ionisation_Chamber_Temperature.push_back(Temperature);
  
  m_NumberOfDetectors++;
    }

    void Add_Plastic(string Material, double Width, double Length,
            double Thickness, std::vector<int> Pos) {
/*        m_Plastic_Material.push_back(Material);
        m_Plastic_Width.push_back(Width);
        m_Plastic_Length.push_back(Length);
        m_Plastic_Thickness.push_back(Thickness);
        m_Plastic_Position.push_back(Pos);
  */
  m_NumberOfDetectors++;
    }
    void Add_Gas_Gap(double Z,double Thickness,string Gas,double Pressure,double Temperature){
        m_Gas_Gap_Z.push_back(Z);
        m_Gas_Gap_Thickness.push_back(Thickness);
        m_Gas_Gap_Gas.push_back(Gas);
        m_Gas_Gap_Pressure.push_back(Pressure);
        m_Gas_Gap_Temperature.push_back(Temperature);
      m_NumberOfDetectors++;
    };

    void Add_AC(double Z, double Thickness, string Material){
      AC_Material.push_back(Material);
      AC_Thickness.push_back(Thickness);
      AC_Z.push_back(Z);
      m_NumberOfDetectors++;
    };
    
    void Add_Entry_Exit(double Z, double Thickness, string Material){
      Entry_Exit_Material.push_back(Material);
      Entry_Exit_Thickness.push_back(Thickness);
      Entry_Exit_Z.push_back(Z);
      m_NumberOfDetectors++;
    };
    
    void Add_Drift_Chamber(double Z, double Thickness, string Gas, double Pressure,
        double Temperature) {
    m_Drift_Chamber_Z.push_back(Z);
    m_Drift_Chamber_Thickness.push_back(Thickness);
    m_Drift_Chamber_Gas.push_back(Gas);
    m_Drift_Chamber_Pressure.push_back(Pressure);
    m_Drift_Chamber_Temperature.push_back(Temperature);
  
  m_NumberOfDetectors++;
    }
  
  
  //////////////////////////////////////////////////////////////
  // specific methods to ZDD array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();


    void Treat_DC();
    // Matching Time and E for IC and Plastic
    void Match_E_T(std::string Detector);

    // PreTreating Energy for IC and Plastic 
    void PreTreatEnergy(std::string Detector, CalibrationManager* Cal);
    
    // Same for time
    void PreTreatTime(std::string Detector, CalibrationManager* Cal);

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TZDDData object to TZDDPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TZDDData* rawDataPointer) {m_EventData = rawDataPointer;}

    void quicksort(std::vector<double>* Z, std::vector<double>* Thickness, std::vector<string>* Material);  
    void quicksort(std::vector<double>* Z, std::vector<double>* Thickness, std::vector<string>* Gas, std::vector<double>* Pressure, std::vector<double>* Temperature);  

  // objects are not written in the TTree
  private:
    TZDDData*         m_EventData;        //!
    TZDDData*         m_PreTreatedData;   //!
    TZDDPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TZDDData* GetRawData()        const {return m_EventData;}
    TZDDData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!
    
    // Nb of detectors counter
    int ICcounter;
    int ACcounter;
    int GGcounter;
    int Entry_Exit_counter;
    int Plasticcounter;

  public:
    int GetICcounter(){return ICcounter;};

  // number of detectors
  private:
    int m_NumberOfDetectors;  //!
    std::vector<string> Entry_Exit_Material; //!
    std::vector<double> Entry_Exit_Thickness; //!
    std::vector<double> Entry_Exit_Z; //!
    
    std::vector<string> AC_Material; //!
    std::vector<double> AC_Thickness; //!
    std::vector<double> AC_Z; //!

    std::vector<double> m_Ionisation_Chamber_Z; //!
    std::vector<double> m_Ionisation_Chamber_Thickness; //!
    std::vector<string> m_Ionisation_Chamber_Gas; //!
    std::vector<double> m_Ionisation_Chamber_Pressure; //!
    std::vector<double> m_Ionisation_Chamber_Temperature; //!
    
    std::vector<double> m_Drift_Chamber_Z; //!
    std::vector<double> m_Drift_Chamber_Thickness; //!
    std::vector<string> m_Drift_Chamber_Gas; //!
    std::vector<double> m_Drift_Chamber_Pressure; //!
    std::vector<double> m_Drift_Chamber_Temperature; //!

    std::vector<double> m_Gas_Gap_Z; //!
    std::vector<double> m_Gas_Gap_Thickness; //!
    std::vector<string> m_Gas_Gap_Gas; //!
    std::vector<double> m_Gas_Gap_Pressure; //!
    std::vector<double> m_Gas_Gap_Temperature; //!

/*    std::vector<double> CATS_X_1;
    std::vector<double> CATS_X_2;
    std::vector<double> CATS_Y_1;
    std::vector<double> CATS_Y_2;

    std::vector<double> CATS_TS_1;
    std::vector<double> CATS_TS_2;*/

  public:
    std::vector<string> Get_Entry_Exit_Material(){return Entry_Exit_Material;};
    std::vector<double> Get_Entry_Exit_Thickness(){return Entry_Exit_Thickness;};
    std::vector<double> Get_Entry_Exit_Z(){return Entry_Exit_Z;};
    
    std::vector<string> Get_AC_Material(){return AC_Material;};
    std::vector<double> Get_AC_Thickness(){return AC_Thickness;};
    std::vector<double> Get_AC_Z(){return AC_Z;};

    std::vector<double>Get_m_Ionisation_Chamber_Z(){return m_Ionisation_Chamber_Z;};
    std::vector<double>Get_m_Ionisation_Chamber_Thickness(){return m_Ionisation_Chamber_Thickness;};
    std::vector<string>Get_m_Ionisation_Chamber_Gas(){return m_Ionisation_Chamber_Gas;};
    std::vector<double>Get_m_Ionisation_Chamber_Pressure(){return m_Ionisation_Chamber_Pressure;};
    std::vector<double>Get_m_Ionisation_Chamber_Temperature(){return m_Ionisation_Chamber_Temperature;};

    std::vector<double>Get_m_Gas_Gap_Z(){return m_Gas_Gap_Z;};
    std::vector<double>Get_m_Gas_Gap_Thickness(){return m_Gas_Gap_Thickness;};
    std::vector<string>Get_m_Gas_Gap_Gas(){return m_Gas_Gap_Gas;};
    std::vector<double>Get_m_Gas_Gap_Pressure(){return m_Gas_Gap_Pressure;};
    std::vector<double>Get_m_Gas_Gap_Temperature(){return m_Gas_Gap_Temperature;};
    
    std::vector<double>Get_m_Drift_Chamber_Z(){return m_Drift_Chamber_Z;};
    std::vector<double>Get_m_Drift_Chamber_Thickness(){return m_Drift_Chamber_Thickness;};
    std::vector<string>Get_m_Drift_Chamber_Gas(){return m_Drift_Chamber_Gas;};
    std::vector<double>Get_m_Drift_Chamber_Pressure(){return m_Drift_Chamber_Pressure;};
    std::vector<double>Get_m_Drift_Chamber_Temperature(){return m_Drift_Chamber_Temperature;};


  // spectra class
  private:
    TZDDSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TZDDPhysics,1)  // ZDDPhysics structure
};
#endif
