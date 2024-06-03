#!/opt/homebrew/opt/python@3.11/libexec/bin/python
import os
import subprocess
import numpy as np


def format_value(value):
    """ Format the value by replacing '-' with 'm' and '.' with '_' """
    return str(value).replace('-', 'm').replace('.', '_')

def modify_files(params, outputFileName, outputDir, beamFile,macroFile):
    # Default parameters
    defaults = {'thetaTarget': 0,'offsetTargetX': 0,'offsetTargetY' :0,
        'offsetBeamX' : 0, 'offsetBeamY' : 0, 'sigmaBeamX': 0, 'sigmaBeamY': 0,
    'holderShift' : 0, 'energyBeam' : 0 }

    # Update params with any defaults that aren't specified
    for key, value in defaults.items():
        params.setdefault(key, value)

    # Define file paths
    musett_file = "/Users/lh270370/Software/nptool/NPSimulation/Detectors/MUSETT/MUSETT.cc"
    beam_file = "./Beam/" + beamFile
    routine_file = "./routine.sh"

    # Function to preserve indentation
    def replace_line(line, new_content):
        leading_spaces = len(line) - len(line.strip())
        return ' ' * leading_spaces + new_content + '\n'

    # Read and modify MUSETT.cc
    with open(musett_file, 'r') as file:
        lines = file.readlines()
    with open(musett_file, 'w') as file:
        for line in lines:
            if "double rotationAngle =" in line:
                new_line = replace_line(line, f"double rotationAngle = {params['thetaTarget']:.8f} * deg;")
                file.write(new_line[1:])
            elif "double offsetTargetX =" in line:
                new_line = replace_line(line, f"double offsetTargetX = {params['offsetTargetX']:.8f} * mm;")
                file.write(new_line[1:])
            elif "double offsetTargetY =" in line:
                new_line = replace_line(line, f"double offsetTargetY = {params['offsetTargetY']:.8f} * mm;")
                file.write(new_line[1:])
            elif "double holder_shift =" in line:
                new_line = replace_line(line, f"double holder_shift = {params['holderShift']:.8f} * mm;")
                file.write(new_line[1:])
            else:
                file.write(line)

    # Read and modify Ra222.beam
    with open(beam_file, 'r') as file:
        lines = file.readlines()
    with open(beam_file, 'w') as file:
        for line in lines:
            if "x0=" in line:
                new_line = replace_line(line, f"x0= {params['offsetBeamX']:.1f}")
                file.write(new_line[1:])
            elif "y0=" in line:
                new_line = replace_line(line, f"y0= {params['offsetBeamY']:.1f}")
                file.write(new_line[1:])
            elif "sigmaX =" in line:
                new_line = replace_line(line, f"sigmaX = {params['sigmaBeamX']:.2f}")
                file.write(new_line[1:])
            elif "sigmaY =" in line:
                new_line = replace_line(line, f"sigmaY = {params['sigmaBeamY']:.2f}")
                file.write(new_line[1:])
            elif "EnergyLow" in line:
                new_line = replace_line(line, f"EnergyLow = {(params['energyBeam']/1000):.4f}")
                file.write(new_line[1:])
            elif "EnergyHigh" in line:
                new_line = replace_line(line, f"EnergyHigh = {(params['energyBeam']/1000):.4f}")
                file.write(new_line[1:])
            else :
                file.write(line)

    # Read and modify routine.sh
    with open(routine_file, 'r') as file:
        lines = file.readlines()
    with open(routine_file, 'w') as file:
        for line in lines:
            if 'SIM_FILE=' in line:
                parts = []
                for param, value in params.items():
                    if value != defaults[param] and param == "energyBeam":
                        format_val = format_value(value)
                        parts.append(f"{param}{format_val}")
                new_sim_file = outputFileName + '__'
                for x in parts :
                    new_sim_file+=x
                new_line = replace_line(line, f'SIM_FILE="{new_sim_file}"')
                file.write(new_line[1:])
            elif 'DIR=' in line :
                new_line = replace_line(line, 'DIR="'+outputDir + '"')
                file.write(new_line[1:])
            elif 'MACRO=' in line:
                new_line = replace_line(line, 'MACRO="'+macroFile + '"')
                file.write(new_line[1:])
            elif 'EVENT_SOURCE=' in line :
                new_line = replace_line(line, 'EVENT_SOURCE="'+beamFile + '"')
                file.write(new_line[1:])
            else:
                file.write(line)

outputFileName = "222Ra_ChangedBR_p3_1e6"
outputDir = "BranchingRatio_test"
beamFile = "Ra222.beam"
macroFile = "run2e5"

params = {'offsetTargetX': 0, 'offsetTargetY': 0,
    'offsetBeamX' : 0,'offsetBeamY' : 0, 'thetaTarget' : 0,
'sigmaBeamX' : 0, 'sigmaBeamY' : 0, 'energyBeam ': 0 }


xBeam = -6
yBeamImpact = 6
yBeam = -150
#energyBeam = 50 # in keV
radius_target = 6
sigmaBeam = 0.2
holder_shift = 0.5
theta = 0
offsetX = - 3
offsetY = yBeamImpact - (offsetX-xBeam)*np.tan(theta*np.pi/180)
if(np.sqrt((offsetY-yBeamImpact)**2+(offsetX-xBeam)**2))>radius_target:
    print("BEAM OUTSIDE TARGET ! ")

for energyBeam in [30]:


    params = {'offsetTargetX': offsetX, 'offsetTargetY': offsetY,
            'offsetBeamX' : xBeam,'offsetBeamY' : yBeam, 'thetaTarget' : theta,
            'sigmaBeamX' : sigmaBeam, 'sigmaBeamY' : sigmaBeam,
             'holderShift' : holder_shift, 'energyBeam' : energyBeam }
    #print(params)
    modify_files(params, outputFileName,outputDir, beamFile, macroFile)
    subprocess.run('./routine.sh', shell=True, executable="/bin/bash", check=True)
