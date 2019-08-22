

header_string = """# Macro file for example B1\n
# \n
# Can be run in batch, without graphic\n
# or interactively: Idle> /control/execute run1.mac\n
#\n
# Change the default number of workers (in multi-threading mode) \n
#/run/numberOfWorkers 4\n
#\n
# Initialize kernel\n
\n
\n
# Set low energy EM processes, including fluorescence\n
/process/em/fluo true\n
/process/em/auger false\n
/process/em/augerCascade false\n
/process/em/pixe true\n
/process/em/deexcitationIgnoreCut false\n
\n
\n
#/cuts/setMaxCutEnergy 50 eV\n
/run/initialize\n
\n
\n
/control/verbose 0\n
/run/verbose 0\n
/event/verbose 0\n
/tracking/verbose 0\n
\n
# Maximum run size:\n
# Summit INT_MAX=2147480000\n"""

numBGparts = "1775000000";
numLCparts = "825000000";
numSigparts = ["12570", "9051", "14834"];

index = ["1","2","3","4"]
energy = ["100keV", "200keV","300keV"]
energy_str = ["100","200","300"]
for e in range(0, len(energy)):
    for i in range(0, len(index)):
        with open("runBeamOn" + energy[e] + "_" + index[i] + ".mac", "w") as f:
            f.write(header_string)
            f.write("#Set hit and signal file names\n/dataCollection/setHitFileName hit" + energy[e] + 
                    "_" + index[i] + ".csv\n/dataCollection/setSignalFileName signal" + energy[e] + "_" + index[i] + ".csv\n")
            f.write("# E0 of background distribution out of {100,200,300} keV\n/particleSource/setFoldingEnergy " + energy_str[e] + " # [keV]")
            f.write("# 0 = background trapped electrons\n/particleSource/setBackgroundType 0\n/run/beamOn " + numBGparts + "\n")
            f.write("# 1 = background loss cone electrons\n/particleSource/setBackgroundType 1\n/run/beamOn "+ numLCparts + "\n")
            f.write("# 2 = bremsstrahlung signal photons\n/particleSource/setBackgroundType 2\n/run/beamOn "+ numSigparts[e] + " \n")

