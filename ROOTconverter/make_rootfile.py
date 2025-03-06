import ROOT
import json
import os
import numpy as np
import time
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')

start_time = time.time()
date = "20230722"

pathToDataFolder = f"/eos/project/f/flic-bi207/bi207/DATA/{date}/jsonData/"
pathToSaveFolder = "/eos/user/j/jcapotor/bi207Data/"

def create_root_file(date, path, fileName):
    jsonFile = open(path + fileName)
    dataFile = json.load(jsonFile)
    print(dataFile.keys())
    rootFile = ROOT.TFile(f"{pathToSaveFolder}{date}.root", "RECREATE")
    rootTree = ROOT.TTree("events", "Tree containing events info")
    channelInfo = {}

    # Iterate over the keys in the first data event to determine array sizes
    first_data_event = dataFile["all"][0]
    for key, value in first_data_event.items():
        if isinstance(value, list):
            array_size = len(value)
        else:
            array_size = 1
        if key == "runUnixTime":
            value = int(value)
        if key == "runDate":
            y, m, d = (value.split("-"))
            value = int(y+m+d)
        if key == "runTime":
            h, min, s = value.split(" ")[3].split(":")
            value = int(h+min+s)
        # Create a NumPy array with appropriate data type and size
        if isinstance(value, list):
            channelInfo[key] = np.zeros(array_size, dtype="int32")
        elif isinstance(value, int):
            if key == "runUnixTime":
                channelInfo[key] = np.zeros(array_size, dtype="int64")
            else:
                channelInfo[key] = np.zeros(array_size, dtype="int32")

        # Create TBranches based on the key and data type
        rootTree.Branch(key, channelInfo[key], f"{key}[{array_size}]/I")

    return rootFile, rootTree, channelInfo

rootFile, rootTree, channelInfo = create_root_file(date, pathToDataFolder, f"{date}_0002.json")
for (dirpath, dirnames, filenames) in os.walk(pathToDataFolder):
    dirnames.sort()
    filenames.sort()
    for fname in filenames:
        print(fname)
        with open(pathToDataFolder + fname) as f:
            data = json.load(f)
            for ievt, dataEvent in enumerate(data["all"]):
                for key in data["all"][ievt].keys():
                    if key == "runUnixTime":
                        data["all"][ievt][key] = int(data["all"][ievt][key])
                    if key == "runDate":
                        y, m, d = (data["all"][ievt][key].split("-"))
                        data["all"][ievt][key] = int(y+m+d)
                    if key == "runTime":
                        h, min, s = data["all"][ievt][key].split(" ")[3].split(":")
                        data["all"][ievt][key] = int(h+min+s)
                    if key in channelInfo:
                        if isinstance(data["all"][ievt][key], list):
                            if len(data["all"][ievt][key]) != len(channelInfo[key]):
                                continue
                            channelInfo[key][:] = data["all"][ievt][key]
                        elif isinstance(data["all"][ievt][key], int):
                            channelInfo[key][0] = data["all"][ievt][key]
                rootTree.Fill()

rootTree.Write()
rootFile.Close()

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
