#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>

#include "TMust2Physics.h"
#include "Utils.cpp"

void CalibrateBeam()
{
    ROOT::EnableImplicitMT();
    //we can use all the data, because its measured between Caviar and CATS
    auto df {ReadAll12BeData()};
    
    //plot original values without telescope hits
    auto hOri {df.Filter([](const TMust2Physics& must2){return must2.TelescopeNumber.size() == 0;}, {"MUST2"})
        .Histo1D({"hOri", "Raw TOF", 1000, 100., 350.}, "T_CATS1_CAV")};
    
    hOri->DrawClone();   
}
