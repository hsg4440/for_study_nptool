#ifndef TEXOGAMPHYSICS_H
#define TEXOGAMPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: S. Giron   contact address: giron@ipno.in2p3.fr          *
 *                  B. Le Crom                  lecrom@ipno.in2p3.fr         *
 * Creation Date  : march 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold exogam treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/


// STL
#include <vector>
#include <map>
#include <algorithm>

// NPL
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPVTreeReader.h"
#include "Geometry_Clover_Exogam.h"
#include "TExogamData.h"
#include "TExogamSpectra.h"
#include "TExogamPhysicsReader.h"
#include "NPInputParser.h"
#include "RootHistogramsCalib.h"
#include "TSpectrum.h"
// ROOT 
#include "TVector2.h" 
#include "TVector3.h" 
#include "TObject.h"
//#include "TH1.h"
//Cubix
#if CUBIX
#include "CXRecalEnergy.h"
#endif

using namespace std ;
// Forward Declaration
class TExogamSpectra;

class TExogamPhysics : public TObject, public NPL::VDetector, public TExogamPhysicsReader{
 public:
  TExogamPhysics()	;
  ~TExogamPhysics() {};

  
 public: 
  void Clear()	              ;	
  void Clear(const Option_t*) {};

  // Only threshold and cal applied to exogam
  std::vector<double> E_Cal;
  std::vector<double> EHG_Cal;
  std::vector<double> Outer1_Cal;
  std::vector<double> Outer2_Cal;
  std::vector<double> Outer3_Cal;
  std::vector<double> Outer4_Cal;
  // std::vector<double> E_Doppler;
  std::vector<unsigned int> FlangeN;
  std::vector<unsigned int> CrystalN;
  
  // Previous treatment + Add_Back (size of vectors are not the same because of AB !)
  std::vector<double> E_AB;
  std::vector<unsigned int> FlangeN_AB;
  std::vector<unsigned int> Size_AB;
  std::vector<unsigned int> CrystalN_ABD;
  std::vector<int> OuterN_ABD;
  std::vector<double> E_ABD;




  Int_t	 		EventMultiplicity	; //!
  Int_t                 ECC_Multiplicity        ; //!
  Int_t                 GOCCE_Multiplicity      ; //!
  Int_t                 NumberOfClover          ; //!
  
  // Clover
  Int_t                 NumberOfHitClover       ; //!
  Int_t                 NumberOfHitCristal      ; //!
  vector<int>		ECC_CloverNumber		;    //!
  vector<int>		ECC_CristalNumber		; //!
  vector<int>		GOCCE_CloverNumber		;    //!
  vector<int>		GOCCE_CristalNumber		; //!
  vector<int>		GOCCE_SegmentNumber		; //!
    	
  //	ECC
  vector<double>	ECC_E				; //!
  vector<double>	ECC_T				; //!
  
  //	GOCCE
  vector<double>	GOCCE_E				; //!
 
  //  Add-Back and Doppler correction
  
  vector<int>      CristalNumber                     ; //!
  vector<int>      SegmentNumber                     ; //!
  vector<int>      CloverNumber                      ; //!
  int              CloverMult                        ; //!
  
  vector<double>   TotalEnergy_lab                   ; //!
  vector<double>   Time                              ; //!
  vector<double>   DopplerCorrectedEnergy            ; //!
  vector<double>   Position                          ; //!
  vector<double>   Theta                             ; //!

  vector < vector < vector < vector <double> > > > Clover_Angles_Theta_Phi;   //!
 
  /* 
  TH1F*                 clover_mult                  ;  
  TH1F*                 cristal_mult                 ;  
  */

 public:		//	Innherited from VDetector Class
			
  //	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
  void ReadConfiguration(NPL::InputParser) 				; //!
		

  //	Add Parameter to the CalibrationManger
  void AddParameterToCalibrationManager()	;		 //!
			
		
  //	Activated associated Branches and link it to the private member DetectorData address
  //	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
  void InitializeRootInputRaw() 					; //!

    //   Activated associated Branches and link it to the private member DetectorPhysics address
    //   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ; //!

  //	Create associated branches and associated private member DetectorPhysics address
  void InitializeRootOutput() 		 		; //!
		
  //	This method is called at each event read from the Input Tree. Aim is to build a tree of calibrated data.
  void PreTreat()			       ; //!

  //	This method is called at each event read from the Input Tree. 
  void BuildPhysicalEvent()					; //!
		
		
  //	Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
  //	This method aimed to be used for analysis performed during experiment, when speed is requiered.
  //	NB: This method can eventually be the same as BuildPhysicalEvent.
  void BuildSimplePhysicalEvent()	       ; //!

  double DopplerCorrection(double Energy, double Theta); //!

  //	Those two method all to clear the Event Physics or Data
  void ClearEventPhysics()		{Clear();}		 //!
  void ClearEventData()			{m_EventData->Clear();}	 //!
  void ClearPreTreatedData()	        {m_PreTreatedData->Clear();} //!

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra(); //!
    // Fill the spectra hold by the spectra class
    void FillSpectra(); //!
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra(); //!
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra(); //!
    
    void SetTreeReader(TTreeReader* TreeReader); //!
  
  
    void InitializeRootHistogramsCalib(); //!
    
    void FillHistogramsCalib(); //!
    
    void DoCalibration();//!

    void InitializeRootHistogramsE_F(unsigned int Detector_Nbr); //!

    void InitializeRootHistogramsEHG_F(unsigned int Detector_Nbr); //!

    // void InitializeRootHistogramsT_F(unsigned int Detector_Nbr); //!
    
    void InitializeRootHistogramsOuter_F(unsigned int Detector_Nbr, unsigned int Outer_Nbr); //!
    
    void FillRootHistogramsCalib_F(); //!
  
    void DoCalibrationE_F(unsigned int  Detector_Nbr,std::string CalibType, ofstream* calib_file, ofstream* dispersion_file, unsigned int Threshold); //!
    
    void DoCalibrationEHG_F(unsigned int  Detector_Nbr, ofstream* calib_file, ofstream* dispersion_file){}; //!
    
    void DoCalibrationOuter_F(unsigned int  Detector_Nbr, unsigned int Outer_Nbr, ofstream* calib_file, ofstream* dispersion_file){}; //!
      
    void MakeInitialCalibFolder(std::string make_folder); //!
    
    void MakeECalibFolders(std::string make_folder); //!
    
    void MakeEHGCalibFolders(std::string make_folder); //!
    
    void MakeOuterCalibFolders(std::string make_folder); //!

    void DefineCalibrationSource(std::string source); //!

    void CreateCalibrationEFiles(ofstream* calib_file,ofstream* dispersion_file); //!
    
    void CreateCalibrationEHGFiles(ofstream* calib_file,ofstream* dispersion_file); //!
    
    void CreateCalibrationOuterFiles(ofstream* calib_file,ofstream* dispersion_file); //!


    void ReadDoCalibration(NPL::InputParser parser); //!

    void WriteHistogramsCalib(); //!
    
    void WriteHistogramsE(); //!


 private:	//	Root Input and Output tree classes

 				
  TExogamData*         m_EventData;        //!
  TExogamData*         m_PreTreatedData;   //!
  TExogamPhysics*      m_EventPhysics;     //!
  

 public:		//	Specific to EXOGAM Array
  //	Add a Clover
  // void AddClover(string AngleFile);
  void AddClover(int Board, int Flange, int Channel0, int Channel1); //!

  Int_t GetClover_Mult()    { return(CloverNumber.size()); } //!
  //  Int_t GetECC_Mult()   { return(ECC_CristalNumber.size()); }
  //  Int_t GetGOCCE_Mult() { return(GOCCE_SegmentNumber.size()); }

  Double_t GetSegmentAnglePhi(int Clover, int Cristal, int Segment)    {return(Clover_Angles_Theta_Phi[Clover][Cristal][Segment][1]);}; //!
  Double_t GetSegmentAngleTheta(int Clover, int Cristal, int Segment)  {return(Clover_Angles_Theta_Phi[Clover][Cristal][Segment][0]);}; //!
 
  // Give and external TMustData object to TExogamPhysics. Needed for online analysis for example.
  void SetRawDataPointer(TExogamData* rawDataPointer) {m_EventData = rawDataPointer;} //!
  // Retrieve raw and pre-treated data
  TExogamData* GetRawData()        const {return m_EventData;} //!
  TExogamData* GetPreTreatedData() const {return m_PreTreatedData;} //!
  
  void ResetPreTreatVariable(); //!

  void ReadAnalysisConfig(); //!

  double ComputeMeanFreePath(double GammaEnergy); //!

  unsigned int GetFlangeNbr(unsigned int crystal_nbr); //!

  int GetMaxOuter(unsigned int EventId); //!
  
  double GetDoppler(double Energy, unsigned int Flange, unsigned int Crystal, unsigned int Outer); //!

  private: // Variables for analysis

  unsigned int m_EXO_Mult;//!
  double m_EXO_OuterUp_RAW_Threshold;//!
  double m_EXO_E_RAW_Threshold;//!
  double m_EXO_E_Threshold;//!
  double m_EXO_EHG_RAW_Threshold;//!
  double m_EXO_TDC_RAW_Threshold;//!
  double EXO_E;//!
  double EXO_EHG;//!
  double EXO_TDC;//!
  double EXO_Outer1;//!
  double EXO_Outer2;//!
  double EXO_Outer3;//!
  double EXO_Outer4;//!
  unsigned int flange_nbr;//!
  unsigned int crystal_nbr;//!
  double Beta = 0.257;//!

  Clover_struc Exogam_struc;//!
  
  std::map<unsigned int,std::pair<unsigned int,unsigned int>> MapCrystalFlangeCLover;//! Map key is raw crystal nbr, pair associated is flange nbr and crystal nbr in the flange

  double GeDensity = 0.005323; //! g/mm3
  std::map<double, double> Map_PhotonCS;//!
 
  map<int, bool> DoCalibrationE;                                    //!
  map<int, bool> DoCalibrationEHG;                                  //!
  map<int, bool> DoCalibrationT;                                    //!
  map<int, map<int,bool>> DoCalibrationOuter;                       //!

  unsigned int Threshold_E_Cal;//!
  unsigned int Threshold_EHG_Cal;//!
  unsigned int Threshold_Outers_Cal;//!
  
  std::vector<std::string> Source_isotope;                          //!
  std::vector<double> Source_E;                                     //!
  std::vector<double> Source_Sig;                                   //!
  std::vector<double> Source_branching_ratio;                       //!
  std::string Source_name;//!
  unsigned int FitPolOrder;//!
  
#if CUBIX
  CXRecalEnergy* CubixEnergyCal = new CXRecalEnergy();
#endif

  private: // Spectra Class   
    TExogamSpectra*      m_Spectra;//! 

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 		 //!

  public: // Static constructor to be passed to the Detector Factory
     static NPL::VDetector* Construct(); //!
     static NPL::VTreeReader* ConstructReader(); //!
     
     ClassDef(TExogamPhysics,1)  // ExogamPhysics structure
    };

namespace EXOGAM_LOCAL
{
  double fEXO_E(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_EHG(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_T(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_Outer(const TExogamData* m_EventData, const unsigned int& i, const unsigned int OuterNumber);
   const double Threshold_ECC   = 50;
   const double Threshold_GOCCE = 0;
   const double RawThreshold_ECC   = 0;
   const double RawThreshold_GOCCE = 0;
   
   
   //	tranform an integer to a string
   string itoa(int value);

 
}


#endif
