#include <Math/Vector3D.h>
#include <Math/Point3D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

#include "TCATSPhysics.h"
#include "TMust2Data.h"
#include "TMust2Physics.h"

using XYZVector = ROOT::Math::XYZVector;
using XYZPoint  = ROOT::Math::XYZPoint;

double LandTtoGamma(double L, double T)
{
    //MUST be in SI units
    double v {L / T};
    double beta {v / TMath::C()};
    std::cout<<"Beta = "<<beta<<'\n';
    double gamma {1. / TMath::Sqrt(1. - beta * beta)};
    return gamma;
}

double ConvertTOFToMass(double Tbeam,
                        double xM2, double yM2, double zM2, double tM2,
                        double tCorr,
                        TCATSPhysics& cats)
{
    double timeWindow {600.};// ns
    double lengthFactor {180.};//mm
    //Length
    XYZPoint vertex {cats.PositionOnTargetX, cats.PositionOnTargetY, 0.};
    //get MUST2 point
    XYZPoint must2Point {xM2, yM2, zM2 + lengthFactor};
    //Time (including correction)
    double tof {tM2 + tCorr};
    std::cout<<"tCorr = "<<tCorr<<'\n';
    auto LAfterReaction {(must2Point - vertex).R()};//should be in mm!
    //let's ignore by now distance cats (or whatever provides STOP signal) to vertex
    auto L {(LAfterReaction) * 1.E-3};//mm to m
    auto T {TMath::Abs(tof - timeWindow) * 1.E-9}; //ns to s
    std::cout<<"L = "<<L<<" T = "<<T<<'\n';
    auto gamma {LandTtoGamma(L, T)};
    std::cout<<"Gamma = "<<gamma<<'\n';
    return Tbeam / (gamma - 1);
}

void ImprovePID()
{
    ROOT::RDataFrame d("GatedTree", "./RootFiles/GatedTree_12Be.root");
    auto description {d.Describe()};
    description.Print();std::cout<<'\n';

    auto df = d.Redefine("TOF", "TOF");
    auto hPID {df.Define("x", "MUST2.Si_E")
        .Histo2D({"hPID", "PID;E_{Si} [MeV];TOF [au]", 1000, 0., 30., 1000, 460., 580.}, "x", "TOF")};
    df = df.Define("ReconstructedMass", ConvertTOFToMass, {"BeamEnergy",
                                                           "X_M2", "Y_M2", "Z_M2", "T_M2",
                                                           "TimeCorr",
                                                           "CATS"});

    auto hMass {df.Histo1D( "ReconstructedMass")};

    auto hTelN {df.Define("TelescopeNumber", [](const TMust2Physics& must2)
    {
        return must2.TelescopeNumber;
    },
            {"MUST2"})
        .Histo1D("TelescopeNumber")};
    
    //plotting
    auto* c1 {new TCanvas("c1")};
    c1->DivideSquare(2);
    c1->cd(1);
    hPID->DrawClone("col");
    c1->cd(2);
    hMass->DrawClone();

    auto* c2 {new TCanvas("c2")};
    c2->DivideSquare(2);
    c2->cd(1);
    hTelN->DrawClone();
}
