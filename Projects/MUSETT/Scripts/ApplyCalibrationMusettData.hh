#include "TMUSETTData.h"
#include "NPCalibrationManager.h"
#include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"

const std::string calibrationFileName = "calibrations/Calibration222Ra.txt";
const std::string inputFileName = "RootR/run_0271.root";
const std::string outputFileName = "RootRCal/run_0271_cal.root";
