import ROOT
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')
import matplotlib.pyplot as plt
from eventHandling.event import Event
import numpy as np
import math 

pathToData = "/eos/user/j/jcapotor/bi207Data/"
fName = "20220511.root"
rootFile = ROOT.TFile(f"{pathToData}{fName}", "READ")
energy = {22:[], 23:[], 24:[], 25:[], 26:[], 27:[]}
channelNumber, x, y = [], [], []
for event in rootFile.Get("events"):
    event0 = Event(event)
    print(f"EventId: {event0.eventID}")
    if event0.eventID > 20000:
        break
    event0.make_voxels()
    event0.make_clusters()
    # print(f"Number of clusters: {event0.nClusters} -> {event0.clusters}")
    if event0.clusters is None:
        continue
    for clusterId, cluster in event0.clusters.items():
        if cluster.nVoxels > 2:
            continue
        for voxel in cluster.voxels:
            if voxel.hitColl is None:
                continue
            elif voxel.hitColl is not None:
                if (voxel.hitColl.channelNumber >= 22) and (voxel.hitColl.channelNumber <= 27):
                    print(voxel.hitColl.channelNumber, voxel.xCord, voxel.yCord)
                    x.append(float(voxel.xCord))
                    y.append(float(voxel.yCord))
                    channelNumber.append(voxel.hitColl.channelNumber)
                    energy[voxel.hitColl.channelNumber].append(voxel.hitColl.charge)
plt.subplot(2,3,1)
plt.hist(energy[22], bins=50)
plt.title("Channel 22")
plt.subplot(2,3,2)
plt.hist(energy[23], bins=50)
plt.title("Channel 23")
plt.subplot(2,3,3)
plt.hist(energy[24], bins=50)
plt.title("Channel 24")
plt.subplot(2,3,4)
plt.hist(energy[25], bins=50)
plt.title("Channel 25")
plt.subplot(2,3,5)
plt.hist(energy[26], bins=50)
plt.title("Channel 26")
plt.subplot(2,3,6)
plt.hist(energy[27], bins=50)
plt.title("Channel 27")
plt.tight_layout()
plt.savefig("equal.png")
