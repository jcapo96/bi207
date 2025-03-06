import ROOT
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from eventHandling.event import Event
from peakFinderPerformance import savePeakFinderPerformance
from sympy import symbols, lambdify

pathToData = "/eos/user/j/jcapotor/bi207Data/"
fName = "20230722.root"
rootFile = ROOT.TFile(f"{pathToData}{fName}", "READ")

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
    MESH_INTERSECTION = plt.subplot(gs_upper_right[0, 0])

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
    if event0.eventID < 45:
        continue
    if event0.eventID > 45:
        break
    nCollHits, deltaC, nInd1Hits, deltaI1, nInd2Hits, deltaI2, deltaT = 0, [], 0, [], 0, [], []
    collCharge, ind1Charge, ind2Charge = 0, 0, 0
    event0.make_voxels()
    if len(event0.voxels) < 30:
        continue
    fig, axes = makeFigure(event=event0)
    color = plt.cm.rainbow(np.linspace(0, 1, len(event0.voxels)))
    savePeakFinderPerformance(event0=event0, pathToData=pathToData, fName=fName)
    for nVoxel, voxel in enumerate(event0.voxels):
        for channel, wire in voxel.mesh.wires.items():
            if (type(wire.function)) == int:
                axes[7].axvline(wire.function, color=color[nVoxel], linewidth=1.0)
            else:
                x = symbols('x')
                xValues = np.linspace(0, 47*5, 1000)
                func = lambdify(x, wire.function, 'numpy')
                axes[7].plot(xValues, func(xValues), color=color[nVoxel], linewidth=1.0)
        if voxel.hitColl is not None:
            nCollHits += 1
            axes[0].plot(voxel.hitColl.Peak.waveform.xPoints,
                voxel.hitColl.Peak.waveform.waveform - voxel.hitColl.Peak.waveform.baseline,
                color="black", linewidth=1.0)
            axes[0].fill_between(voxel.hitColl.Peak.waveform.xPoints,
                     voxel.hitColl.Peak.waveform.waveform -  voxel.hitColl.Peak.waveform.baseline,
                     0,
                     where = [voxel.hitColl.Peak.tini < x < voxel.hitColl.Peak.tend for x in voxel.hitColl.Peak.waveform.xPoints],
                     alpha=0.5)
            for hit in voxel.hits:
                if hit is not None:
                    if hit.channel != voxel.hitColl.channel:
                        axes[6].scatter(voxel.hitColl.channelNumber, hit.channelNumber, color="black")
        if voxel.hitInd1 is not None:
            nInd1Hits += 1
            axes[2].plot(voxel.hitInd1.Peak.waveform.xPoints,
                voxel.hitInd1.Peak.waveform.waveform - voxel.hitInd1.Peak.waveform.baseline,
                color="black", linewidth=1.0)
            axes[2].fill_between(voxel.hitInd1.Peak.waveform.xPoints,
                     voxel.hitInd1.Peak.waveform.waveform -  voxel.hitInd1.Peak.waveform.baseline,
                     0,
                     where = [voxel.hitInd1.Peak.tini < x < voxel.hitInd1.Peak.tend for x in voxel.hitInd1.Peak.waveform.xPoints],
                     alpha=0.5)
        if voxel.hitInd2 is not None:
            nInd2Hits += 1
            axes[4].plot(voxel.hitInd2.Peak.waveform.xPoints,
                voxel.hitInd2.Peak.waveform.waveform - voxel.hitInd2.Peak.waveform.baseline,
                color="black", linewidth=1.0)
            axes[4].fill_between(voxel.hitInd2.Peak.waveform.xPoints,
                     voxel.hitInd2.Peak.waveform.waveform -  voxel.hitInd2.Peak.waveform.baseline,
                     0,
                     where = [voxel.hitInd2.Peak.tini < x < voxel.hitInd2.Peak.tend for x in voxel.hitInd2.Peak.waveform.xPoints],
                     alpha=0.5)

        if None not in voxel.positionVector:
            axes[1].scatter(voxel.positionVector[0], voxel.positionVector[1],
                            s=20, c=voxel.charge, cmap='viridis', ec='k',
                            norm = mcolors.Normalize(vmin=0, vmax=1000))
            axes[3].scatter(voxel.positionVector[0], voxel.positionVector[2],
                            s=20, c=voxel.charge, cmap='viridis', ec='k',
                            norm = mcolors.Normalize(vmin=0, vmax=1000))
            axes[5].scatter(voxel.positionVector[1], voxel.positionVector[2],
                            s=20, c=voxel.charge, cmap='viridis', ec='k',
                            norm = mcolors.Normalize(vmin=0, vmax=1000))
            axes[1].errorbar(x=voxel.positionVector[0], y=voxel.positionVector[1],
                            xerr=voxel.positionVectorErr[0], yerr=voxel.positionVectorErr[1],
                            fmt="None", ecolor="black", capsize=5)
            axes[3].errorbar(x=voxel.positionVector[0], y=voxel.positionVector[2],
                            xerr=voxel.positionVectorErr[0], yerr=voxel.positionVectorErr[2],
                            fmt="None", ecolor="black", capsize=5)
            axes[5].errorbar(x=voxel.positionVector[1], y=voxel.positionVector[2],
                            xerr=voxel.positionVectorErr[1], yerr=voxel.positionVectorErr[2],
                            fmt="None", ecolor="black", capsize=5)
            axes[7].scatter(voxel.positionVector[0], voxel.positionVector[1], color=color[nVoxel], marker="*", ec='k')


    collectionWaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nCollHits}")
    axes[0].legend(handles=[collectionWaveformPatch])
    induction1WaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nInd1Hits}")
    axes[2].legend(handles=[induction1WaveformPatch])
    induction2WaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nInd2Hits}")
    axes[4].legend(handles=[induction2WaveformPatch])
    fig.savefig(f'{pathToData}/Plots/{fName.split(".")[0]}/voxelisationPerformance/display{event0.eventID}.pdf', dpi=300, format="pdf")
    for ax in axes:
        ax.cla()
    plt.close(fig)
