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

def makeFigure(event):
    PEAK_FINDER_PERFORMANCE = plt.figure(figsize=(15, 8))

    # Create a GridSpec with 3 rows and 3 columns
    gs_master = GridSpec(4, 3, height_ratios=[1, 1, 1, 4], width_ratios=[1, 1, 1])

    # Create subplots for the first three rows sharing the same X-axis
    gs_upper = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_master[:3, :], hspace=0)
    COLLECTION_PLANE_WAVEFORM = plt.subplot(gs_upper[0, 0])
    INDUCTION1_PLANE_WAVEFORM = plt.subplot(gs_upper[1, 0], sharex=COLLECTION_PLANE_WAVEFORM)
    INDUCTION2_PLANE_WAVEFORM = plt.subplot(gs_upper[2, 0], sharex=COLLECTION_PLANE_WAVEFORM)
    TOP_AXES = [COLLECTION_PLANE_WAVEFORM, INDUCTION1_PLANE_WAVEFORM, INDUCTION2_PLANE_WAVEFORM]
    for ax in TOP_AXES[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    # Create subplots for the three columns sharing the same Y-axis
    gs_lower = GridSpecFromSubplotSpec(2, 3, subplot_spec=gs_master[3, :], height_ratios=[1,2], wspace=0, hspace=0)
    COLLECTION_PLANE_MAP = plt.subplot(gs_lower[1, 0])
    INDUCTION1_PLANE_MAP = plt.subplot(gs_lower[1, 1], sharey=COLLECTION_PLANE_MAP)
    INDUCTION2_PLANE_MAP = plt.subplot(gs_lower[1, 2], sharey=COLLECTION_PLANE_MAP)
    COLLECTION_PLANE_CHARGE = plt.subplot(gs_lower[0, 0], sharex=COLLECTION_PLANE_MAP)
    INDUCTION1_PLANE_CHARGE = plt.subplot(gs_lower[0, 1], sharex=INDUCTION1_PLANE_MAP, sharey=COLLECTION_PLANE_CHARGE)
    INDUCTION2_PLANE_CHARGE = plt.subplot(gs_lower[0, 2], sharex=INDUCTION2_PLANE_MAP, sharey=COLLECTION_PLANE_CHARGE)
    BOTTOM_UP_AXES = [COLLECTION_PLANE_CHARGE, INDUCTION1_PLANE_CHARGE, INDUCTION2_PLANE_CHARGE]
    BOTTOM_DOWN_AXES = [COLLECTION_PLANE_MAP, INDUCTION1_PLANE_MAP, INDUCTION2_PLANE_MAP]
    for ax in BOTTOM_DOWN_AXES[1:]:
        plt.setp(ax.get_yticklabels(), visible=False)
    for i, ax in enumerate(BOTTOM_UP_AXES):
        if i != 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)

    # Create titles, labels and limits for the axes

    INDUCTION2_PLANE_WAVEFORM.set_xlabel(fr"Time Ticks (0.5 $\mu$s)")
    COLLECTION_PLANE_MAP.set_ylabel(fr"Time Ticks (0.5 $\mu$s)")
    COLLECTION_PLANE_CHARGE.set_ylabel(fr"Charge (ADC)")
    PEAK_FINDER_PERFORMANCE.suptitle(fr"Event #{event.eventID}")
    COLLECTION_PLANE_MAP.set_xlim(0, 47)
    INDUCTION1_PLANE_MAP.set_xlim(48, 87)
    INDUCTION2_PLANE_MAP.set_xlim(88, 127)
    PEAK_FINDER_PERFORMANCE.tight_layout(pad=0.5)
    axes = [COLLECTION_PLANE_WAVEFORM, COLLECTION_PLANE_MAP, INDUCTION1_PLANE_WAVEFORM,
            INDUCTION1_PLANE_MAP, INDUCTION2_PLANE_WAVEFORM, INDUCTION2_PLANE_MAP,
            COLLECTION_PLANE_CHARGE, INDUCTION1_PLANE_CHARGE, INDUCTION2_PLANE_CHARGE]

    return PEAK_FINDER_PERFORMANCE, axes

def savePeakFinderPerformance(event0, pathToData, fName):
    fig, axes = makeFigure(event=event0)
    eventHits = event0.hits
    nCollHits, deltaC, nInd1Hits, deltaI1, nInd2Hits, deltaI2, deltaT = 0, [], 0, [], 0, [], []
    collCharge, ind1Charge, ind2Charge = 0, 0, 0
    for channelNumber, hits in eventHits.items():
        for nHit, hit in hits.items():
            if channelNumber in event0.collectionChannelsID:
                nCollHits += 1
                collCharge += hit.charge
                deltaC.append(hit.channelNumber)
                deltaT.append(hit.hitTime)
                axes[0].plot(hit.Peak.waveform.xPoints,
                             hit.Peak.waveform.waveform - hit.Peak.waveform.baseline,
                             color="black", linewidth=1.0)
                axes[0].fill_between(hit.Peak.waveform.xPoints,
                     hit.Peak.waveform.waveform -  hit.Peak.waveform.baseline,
                     0,
                     where = [hit.Peak.tini < x < hit.Peak.tend for x in hit.Peak.waveform.xPoints],
                     alpha=0.5)

                axes[1].scatter(hit.channelNumber, hit.hitTime, s=hit.charge, c=hit.charge, cmap='viridis', ec='k',
                                norm = mcolors.Normalize(vmin=0, vmax=1000))
                axes[6].bar(x=hit.channelNumber, height=hit.charge, width=1,
                            color="green", alpha=0.5)
            elif channelNumber in event0.induction1ChannelsID:
                nInd1Hits += 1
                ind1Charge += hit.charge
                deltaI1.append(hit.channelNumber)
                axes[2].plot(hit.Peak.waveform.waveform, color="black", linewidth=1.0)
                axes[2].fill_between(hit.Peak.waveform.xPoints,
                     hit.Peak.waveform.waveform,
                     hit.Peak.waveform.baseline,
                     where = [hit.Peak.tini < x < hit.Peak.tend for x in hit.Peak.waveform.xPoints],
                     alpha=0.5)
                axes[3].scatter(hit.channelNumber, hit.hitTime, s=hit.charge, c=hit.charge, cmap='viridis', ec='k',
                                norm = mcolors.Normalize(vmin=0, vmax=1000))
                axes[7].bar(x=hit.channelNumber, height=hit.charge, width=1,
                            color="green", alpha=0.5)
            elif channelNumber in event0.induction2ChannelsID:
                nInd2Hits += 1
                ind2Charge += hit.charge
                deltaI2.append(hit.channelNumber)
                axes[4].plot(hit.Peak.waveform.waveform, color="black", linewidth=1.0)
                axes[4].fill_between(hit.Peak.waveform.xPoints,
                     hit.Peak.waveform.waveform,
                     hit.Peak.waveform.baseline,
                     where = [hit.Peak.tini < x < hit.Peak.tend for x in hit.Peak.waveform.xPoints],
                     alpha=0.5)
                axes[5].scatter(hit.channelNumber, hit.hitTime, s=hit.charge, c=hit.charge, cmap='viridis', ec='k',
                                norm = mcolors.Normalize(vmin=0, vmax=1000))
                axes[8].bar(x=hit.channelNumber, height=hit.charge, width=1,
                            color="green", alpha=0.5)
    collectionWaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nCollHits}")
    axes[0].legend(handles=[collectionWaveformPatch])
    induction1WaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nInd1Hits}")
    axes[2].legend(handles=[induction1WaveformPatch])
    induction2WaveformPatch = mpatches.Patch(color="grey", label=fr"#Hits = {nInd2Hits}")
    axes[4].legend(handles=[induction2WaveformPatch])

    deltaT = np.round(np.max(deltaT) - np.min(deltaT), 2)*0.5
    deltaC = np.max(deltaC) - np.min(deltaC) + 1 #starts from 1
    deltaI1 = np.max(deltaI1) - np.min(deltaI1)
    deltaI2 = np.max(deltaI2) - np.min(deltaI2)
    deltaL = np.round(np.sqrt((deltaT*(1.55))**2 + (4/3)*((deltaC*5)**2 + (deltaI1*7.5)**2 - (deltaC*5*7.5*deltaI1))),2)
    deltaZ = np.round(deltaT/deltaC,2)
    collectionWaveformPatch = mpatches.Patch(color="dimgray", label=fr"$\Delta$C = {deltaC} channels")
    timeWaveformPatch = mpatches.Patch(color="darkgray", label=fr"$\Delta$T = {deltaT} $\mu$s")
    lengthWaveformPatch = mpatches.Patch(color="silver", label=fr"$\Delta$L = {deltaL} mm")
    projectionWaveformPatch = mpatches.Patch(color="gainsboro", label=fr"$\Delta$z = {deltaZ} mm")
    axes[1].legend(handles=[collectionWaveformPatch, timeWaveformPatch, lengthWaveformPatch, projectionWaveformPatch])
    induction1WaveformPatch = mpatches.Patch(color="gray", label=fr"$\Delta$Ind1 = {deltaI1} channels")
    axes[3].legend(handles=[induction1WaveformPatch])
    induction2WaveformPatch = mpatches.Patch(color="gray", label=fr"$\Delta$Ind2 = {deltaI2} channels")
    axes[5].legend(handles=[induction2WaveformPatch])

    collectionChargePatch = mpatches.Patch(color="green", label=fr"$\sum$ = {int(np.round(collCharge, -1))} ADC")
    axes[6].legend(handles=[collectionChargePatch])
    induction1ChargePatch = mpatches.Patch(color="green", label=fr"$\sum$ = {int(np.round(ind1Charge, -1))} ADC")
    axes[7].legend(handles=[induction1ChargePatch])
    induction2ChargePatch = mpatches.Patch(color="green", label=fr"$\sum$ = {int(np.round(ind2Charge, -1))} ADC")
    axes[8].legend(handles=[induction2ChargePatch])
    fig.savefig(f'{pathToData}/Plots/{fName.split(".")[0]}/peakFinderPerformance/reco{event0.eventID}.pdf', dpi=300, format="pdf")
    # clear the axes that will be used for the ADC superposition
    for ax in axes:
        ax.cla()

    axes[1].imshow(np.matrix(event0.ADC[0:48]).T, cmap='viridis', aspect='auto', origin="lower")
    axes[3].imshow(np.matrix(event0.ADC[48:98]).T, cmap='viridis', aspect='auto', origin="lower")
    axes[5].imshow(np.matrix(event0.ADC[98:]).T, cmap='viridis', aspect='auto', origin="lower")
    for channelNumber, ADCrow in enumerate(event0.ADC):
        if f"chn{channelNumber}" in event0.collectionChannelsID:
            axes[0].plot(ADCrow, color="black")
        elif f"chn{channelNumber}" in event0.induction1ChannelsID:
            axes[2].plot(ADCrow, color="black")
        elif f"chn{channelNumber}" in event0.induction2ChannelsID:
            axes[4].plot(ADCrow, color="black")

    fig.savefig(f'{pathToData}/Plots/{fName.split(".")[0]}/peakFinderPerformance/raw{event0.eventID}.pdf', dpi=300, format="pdf")



