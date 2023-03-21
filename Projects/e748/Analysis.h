#ifndef Analysis_h 
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2025                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include"NPVAnalysis.h"
#include "TMust2Physics.h"
#include "TCATSPhysics.h"
#include "TModularLeafPhysics.h"
#include "TInitialConditions.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "TRandom3.h"
#include <Math/Point3D.h>

#include <Math/Vector3D.h>
#include <RtypesCore.h>
#include <TAttMarker.h>
#include <unordered_map>
#include <vector>
class Analysis: public NPL::VAnalysis{
public:
    Analysis();
    ~Analysis();

    //for points
    using XYZPoint = ROOT::Math::XYZPoint;
    using XYZVector = ROOT::Math::XYZVector;
public: 
    void Init();
    void TreatEvent();
    void End();
    void InitOutputBranch();
    void InitInputBranch();
    void ReInitValue();
    static NPL::VAnalysis* Construct();
private:
    void ResetMVariables();
    void TreatBeam();
    double ComputeXYZVectorAngle(const XYZVector& v1,
                                 const XYZVector& v2);
    double ComputeGamma(double vInSIUnits);
    double ComputeELossInCATS(double initEBeam, double beamAngle);
    double ComputeTimeCorrectionInCATS(double EBeam, double MBeam,
                                       const XYZPoint& p0, const XYZPoint& p1);
private:
    //Naming convention:
    //-- m members are going to be written to .root file
    //-- f members just used in the class
    //////////////////////////////////////
    // Kinematic variables for each track
    std::vector<double> mEx;
    std::vector<double> mELab;
    std::vector<double> mThetaLab;
    std::vector<double> mThetaCM;
    //std::vector<double> mNormalThetaM2;
    // Kinematics and variables for beam
    double mEBeam;
    double mThetaBeam;
    //double mNormalThetaTarget;
    double mTimeCorr;
    // MUST2 hits;
    // WARNING: Not all info is contained in TMUST2Physics (//!)
    std::vector<short> mMust2Telescopes;
    std::vector<double> mMust2SiE;
    std::vector<double> mMust2CsIE;
    std::vector<double> mMust2SiT;
    //points are decomposed by coords. bc ROOT doesnt
    // know how to read and RVec of XYZPoints and I dont want to play
    //with dictionaries...
    std::vector<double> mMust2PointsX;
    std::vector<double> mMust2PointsY;
    std::vector<double> mMust2PointsZ;
    //Lengths and TOF
    // double mBeamLength;
    // double mTrackLength;
    int mRunMinor;
    int mRunMajor;
    short mCATS1Calibrated;
    short mGATCONF;
    
    ///////////////////////////////////////
    NPL::Reaction* fReaction;
    NPL::EnergyLoss fHe3CD2;
    NPL::EnergyLoss fHe3Al;
    NPL::EnergyLoss fHe3Si;
    NPL::EnergyLoss fBeamCD2;
    NPL::EnergyLoss fBeamMylar;
    NPL::EnergyLoss fBeamIsobutane;

    double fTargetThickness;
    double kBeamMass {12.026922 * 931.494 - 4 * 0.511};//12Be
    
    TMust2Physics* fM2;
    TCATSPhysics* fCATS;
    TModularLeafPhysics* fModularLeaf;
    TInitialConditions* fInitial;
    // //other variables 
    // Short_t         vADC_CHIO_V15;
    // Short_t         vADC_VOIE_29;
    // Short_t         vCHIO;
    // Short_t         vCONFDEC;
    // Short_t         vCONFDEC_AGAVA;
    // Short_t         vDATATRIG;
    // Short_t         vDATATRIG_CHIO;
    // Short_t         vE1D6;
    // Short_t         vE2D6;
    // Short_t         vED4;
    // Short_t         vEXL_HF;
    // Short_t         vGALD4X;
    // Short_t         vGALD4Y;
    // Short_t         vGATCONF;
    // Short_t         vQCaviar;
    // Short_t         vQPlast;
    // Short_t         vTCAVHF;
    // Short_t         vTE1D6CAV;
    // Short_t         vTE1D6GAL;
    // Short_t         vTE1D6HF;
    // Short_t         vTED4HF;
    // Short_t         vTGALD4HF;
    // Short_t         vT_CATS1_2;
    // Short_t         vT_CATS1_CAV;
    // Short_t         vT_CATS1_CAV_Cal;
    // Short_t         vT_MUVI_CATS1;
    // Short_t         vT_PL_CATS1;
    // Short_t         vT_PL_CATS2;
    // Short_t         vT_PL_CHIO;
    // Short_t         vT_PLchCATS1;


};
#endif
