#ifndef ISOLDEMERGING3_HH
#define ISOLDEMERGING3_HH

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TGraph.h>
#include <random>
#include <unistd.h>

#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>
#include <utility> // For std::pair
#include <cstdio>  // for std::remove()
#include <filesystem>
#include <iomanip> // For std::setw and std::setfill

#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include <TTreeIndex.h>
#include "TMath.h"
#include "TMUSETTData.h"
#include "TH1F.h"
#include <TRandom3.h>





// Declare the global variables as extern
double t_max = 1e+09;
double coinc_time = 1000;
double T_ADC = 1000;
std::string FileName = "/Users/lh270370/Software/nptool/Outputs/Simulation/Isolde/BranchingRatio_test/222Ra_ChangedBR_p3_1e6__energyBeam30";

// Function to read the configuration file
void readConfig(const std::string& filename);



class MusettSimuEvent {
private:
    // Private data members
    double energy;
    double time;
    double totalTime;
    double timeStamp;
    int strip;
    int det;
    bool x;
    bool y;
    bool treated;
    bool interstrip;

public:
    // Constructor
    MusettSimuEvent() : energy(0), time(0), totalTime(0),timeStamp(0), strip(0), det(0), x(false), y(false), treated(false), interstrip(false) {}

    // Public setter methods
    void SetEnergy(double e) { energy = e; }
    void SetTime(double t) { time = t; }
    void SetTimeStamp(double t) { timeStamp = t; }
    void SetTotalTime(double tt) { totalTime = tt; }
    void SetStrip(int s) { strip = s; }
    void SetDet(int d) { det = d; }
    void SetX(bool flag) { x = flag; }
    void SetY(bool flag) { y = flag; }
    void SetInterstrip(bool i){interstrip = i;}
    void SetTreated(bool t){treated = t;}

    void AddEnergy(double e) { energy += e; }

    // Public getter methods
    double GetEnergy() const { return energy; }
    double GetTime() const { return time; }
    double GetTimeStamp() const { return timeStamp; }
    double GetTotalTime() const { return totalTime; }
    int GetStrip() const { return strip; }
    int GetDet() const { return det; }
    bool IsX() const { return x; }
    bool IsY() const { return y; }
    bool IsInterstrip() const { return interstrip; }
    bool IsEmpty() const {return x == y;}
    bool IsTreated() const {return treated ;}

    // FillEvent methods for X and Y events
    void FillEvent_X(const TMUSETTData* Mumu, int index) {
        SetEnergy(Mumu->GetDSSDXEEnergy(index));
        SetTime(Mumu->GetDSSDXTTime(index));
        SetTimeStamp(t_max*Mumu->GetDSSDX_TStamp(index));
        SetTotalTime(GetTime() + t_max * Mumu->GetDSSDX_TStamp(index));
        SetStrip(Mumu->GetDSSDXEStripNbr(index));
        SetDet(Mumu->GetDSSDXEDetectorNbr(index));
        SetX(true);
        SetY(false);
        SetInterstrip(Mumu->GetDSSDX_Interstrip(index));
    }

    void FillEvent_Y(const TMUSETTData* Mumu, int index) {
        SetEnergy(Mumu->GetDSSDYEEnergy(index));
        SetTime(Mumu->GetDSSDYTTime(index));
        SetTimeStamp(t_max*Mumu->GetDSSDY_TStamp(index));
        SetTotalTime(GetTime() + t_max * Mumu->GetDSSDY_TStamp(index));
        SetStrip(Mumu->GetDSSDYEStripNbr(index));
        SetDet(Mumu->GetDSSDYEDetectorNbr(index));
        SetX(false);
        SetY(true);
        SetInterstrip(Mumu->GetDSSDY_Interstrip(index));
    }

    void Clear(){
      SetEnergy(0);
      SetTime(0);
      SetTimeStamp(0);
      SetTotalTime(0);
      SetStrip(0);
      SetDet(0);
      SetX(false);
      SetY(false);
      SetTreated(false);
      SetInterstrip(false);
    }
};

class TimeComparator {
public:
    bool operator() (const MusettSimuEvent* lhs, const MusettSimuEvent* rhs) const {
        return lhs->GetTotalTime() < rhs->GetTotalTime();
    }
};

double sumE2(double E1, double E2, double dt) {
  // E1 = trigger signal
  // E2 = piling-up signal
  // dt = T(E2) - T(E1)
  if(dt < T_ADC)
    //return E1 + E2;

    return  E1 + E2*(1-TMath::Exp(5*(dt/T_ADC-1)))/(1-TMath::Exp((double)-5));
  else
    return E1;
}

double now() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}


const int NDet = 4; // Assuming there are 4 detectors
const int NStrip = 128; // Assuming each detector has 128 strips
constexpr int MaxMultiplicity = 10; // Maximum number of events you expect in any vector


std::vector<std::vector<std::vector<MusettSimuEvent*>>> SameTimeEvents;
std::vector<std::set<std::pair<int, int>>> HitStrips;

#endif // ISOLDEMERGING3_HH
