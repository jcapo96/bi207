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

def makeFigure(event):
    VOXEL_FINDER_PERFORMANCE = plt.figure(figsize=(15, 8))

    # Create a GridSpec with 3 rows and 3 columns
    gs_master = GridSpec(4, 4, height_ratios=[1, 1, 1, 4], width_ratios=[1, 1, 1, 1])

    # Create subplots for the first three rows sharing the same X-axis
    gs_upper = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_master[:3, :2], hspace=0)
    COLLECTION_PLANE_WAVEFORM = plt.subplot(gs_upper[0, 0])
    INDUCTION1_PLANE_WAVEFORM = plt.subplot(gs_upper[1, 0], sharex=COLLECTION_PLANE_WAVEFORM)
    INDUCTION2_PLANE_WAVEFORM = plt.subplot(gs_upper[2, 0], sharex=COLLECTION_PLANE_WAVEFORM)
    TOP_AXES = [COLLECTION_PLANE_WAVEFORM, INDUCTION1_PLANE_WAVEFORM, INDUCTION2_PLANE_WAVEFORM]
    for ax in TOP_AXES[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    gs_upper_right = GridSpecFromSubplotSpec(1, 1, subplot_spec=gs_master[:3, 2:])
    MESH_INTERSECTION = plt.subplot(gs_upper_right[0, 0], projection='3d')
    MESH_INTERSECTION.set_proj_type('persp')

    # Create subplots for the three columns sharing the same Y-axis
    gs_lower = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_master[3, :], height_ratios=[1,1])
    XY_MAP = plt.subplot(gs_lower[1, 0])
    XZ_MAP = plt.subplot(gs_lower[0, 0], sharex=XY_MAP)
    YZ_MAP = plt.subplot(gs_lower[1, 1], sharey=XZ_MAP)
    DEGENERACY_MAP = plt.subplot(gs_lower[0, 1])
    MAP_AXES = [XY_MAP, XZ_MAP, YZ_MAP]
    for i, ax in enumerate(MAP_AXES):
        if i == 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        if i == 0:
            continue
        if i == 2:
            plt.setp(ax.get_yticklabels(), visible=False)

    # Create titles, labels and limits for the axes
    XY_MAP.set_xlim(0, 5*48)
    XY_MAP.set_ylim(0, 320)
    XZ_MAP.set_ylim(0, 344*1.55)
    YZ_MAP.set_xlim(0, 320)
    MESH_INTERSECTION.set_ylim(0, 320)
    MESH_INTERSECTION.set_xlim(0, 5*48)
    INDUCTION2_PLANE_WAVEFORM.set_xlabel(fr"Time Ticks (0.5 $\mu$s)")
    VOXEL_FINDER_PERFORMANCE.tight_layout(pad=0.5)
    axes = [COLLECTION_PLANE_WAVEFORM, XY_MAP, INDUCTION1_PLANE_WAVEFORM,
            XZ_MAP, INDUCTION2_PLANE_WAVEFORM, YZ_MAP,
            DEGENERACY_MAP, MESH_INTERSECTION]

    return VOXEL_FINDER_PERFORMANCE, axes


for event in rootFile.Get("events"):
    event0 = Event(event)
    print(f"EventId: {event0.eventID}")
    if event0.eventID <= 11:
        continue
    if event0.eventID > 12:
        break
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