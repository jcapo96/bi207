import ROOT
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from eventHandling.event import Event
from peakFinderPerformance import savePeakFinderPerformance

pathToData = "/eos/user/j/jcapotor/bi207Data/"
fName = "20230722.root"
rootFile = ROOT.TFile(f"/Users/jcapo/cernbox/DUNE-IFIC/Software/bi207/20220515.root", "READ")

for event in rootFile.Get("events"):
    event0 = Event(event)
    print(f"EventId: {event0.eventID}")
    # if event0.eventID <= 11:
    #     continue
    # if event0.eventID > 12:
    #     break
    event0.make_voxels()
    event0.make_clusters()
    print(f"Number of clusters: {event0.nClusters} -> {event0.clusters}")
    color = plt.cm.rainbow(np.linspace(0, 1, len(event0.clusters)))
    color_map = {clusterId: color[i] for i, clusterId in enumerate(event0.clusters.keys())}
    if event0.clusters is None:
        continue
    # fig, axes = makeFigure(event=event0)
    fig = plt.figure()
    axes = fig.add_subplot(projection='3d') # fig.add_subplot(111, projection='3d')
    axes.set_xlim(0, 5*48)
    axes.set_ylim(0, 320)
    axes.set_zlim(0, 344*1.55)

    axes.set_xlabel("X (cm)")
    axes.set_ylabel("Y (cm)")
    axes.set_zlabel("Z (cm)")

    axes.set_title(f"Event {event0.eventID}")
    for clusterId, cluster in event0.clusters.items():
        print(clusterId)
        for voxel in cluster.voxels:
            # Original 3D points
            axes.scatter(voxel.xCord, voxel.yCord, voxel.zCord,
                         color=color_map[clusterId], alpha=0.7)

            # Add 2D projections as "shadows"
            axes.scatter(voxel.xCord, voxel.yCord, 0,
                         color="gray", alpha=0.3, marker="o")  # Shadow on XY plane
            axes.scatter(voxel.xCord, 320, voxel.zCord,
                         color="gray", alpha=0.3, marker="o")  # Shadow on XZ plane
            axes.scatter(0, voxel.yCord, voxel.zCord,
                         color="gray", alpha=0.3, marker="o")  # Shadow on YZ plane

    plt.show()
    #fig.savefig(f"exampleCLUST_event{event0.eventID}.png", dpi=300)