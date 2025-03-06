from eventHandling.cluster import Cluster
from eventHandling.voxel import Voxel
from eventHandling.hit import Hit
from eventHandling.peak import Peak
from eventHandling.waveform import Waveform
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
import json
class Event():
    def __init__(self, event):
        self.event = event
        self.settings =self.settings = json.load(open("/Users/jcapo/cernbox/DUNE-IFIC/Software/bi207/analysis/eventHandling/settings.json"))
        self.ADC = np.zeros((128, len(self.event.chn0)))
        self.hits = {}
        self.nCollChannels = 48
        self.nInd1Channels = 40
        self.nInd2Channels = 40
        self.matchingTime = self.settings["hit"]["minPeakSeparation"]
        self.voxelSeparation = self.settings["hit"]["voxelSeparation"]
        self.nChannels = self.nCollChannels + self.nInd1Channels + self.nInd2Channels
        self.set_Eventparameters()
        self.set_ADC()
        self.find_hits()
        self.reshape_channels()
        # self.calculate_track_angle()

    def set_Eventparameters(self):
        self.runUnixTime = self.event.runUnixTime
        self.runTime = self.event.runTime
        self.runDate = self.event.runDate
        self.eventTime = self.event.eventTime
        self.collectionBiasVoltage = self.event.COLLECTION_BIAS_VOLTAGE
        self.cathodeHV = self.event.CATHODE_HV
        self.binaryFileID = self.event.binaryFileID
        self.convertedFileID = self.event.convertedFileID
        self.eventID = self.event.eventId
        self.binaryEventID = self.event.binaryEventID
        self.convertedEventID = self.event.convertedEventID
        self.channelsID = ["chn"+str(i) for i in range(0,self.nChannels)]
        self.collectionChannelsID = ["chn"+str(i) for i in range(0, self.nCollChannels)]
        self.induction1ChannelsID = ["chn"+str(i+self.nCollChannels) for i in range(0, self.nInd1Channels)]
        self.induction2ChannelsID = ["chn"+str(i+self.nCollChannels+self.nInd1Channels) for i in range(0, self.nInd2Channels)]

    def set_ADC(self):
        for nchan, channel in enumerate(self.channelsID):
            self.ADC[nchan] = Waveform(getattr(self.event, channel)).waveform  # Waveform.waveform returns an array

    def find_hits(self):
        for nchan, channel in enumerate(self.channelsID):
            self.hits[channel] = {}
            waveform = Waveform(getattr(self.event, channel))
            for index, peak in enumerate(waveform.peaks):
                if len(self.hits[channel]) == 0:
                    self.hits[channel][index] = Hit(channel, Peak(waveform.waveform, peak))
                    peak_index = index
                else:
                    if abs(Hit(channel, Peak(waveform.waveform, peak)).hitTime - self.hits[channel][peak_index].hitTime) > self.matchingTime:
                        self.hits[channel][index] = Hit(channel, Peak(waveform.waveform, peak))
                        peak_index = index
                    else:
                        continue
        return self.hits

    def find_nhits_per_channel(self):
        if len(self.hits) == 0:
            self.find_hits()
        self.nhits = {}
        self.nHits = 0
        for nchan, channel in enumerate(self.channelsID):
            self.nhits[channel] = (len(self.hits[channel]))
            self.nHits += len(self.hits[channel])
        return self.nhits, self.nHits

    def reshape_channels(self):
        self.reshapedHits = {}
        for nchan, channel in enumerate(self.collectionChannelsID):
            self.reshapedHits[f"chn{self.nCollChannels-1-nchan}"] = self.hits[channel]
        for index, channel in enumerate(self.induction1ChannelsID):
            self.reshapedHits[self.induction1ChannelsID[index]] = self.hits[self.induction2ChannelsID[index]]
            self.reshapedHits[self.induction2ChannelsID[index]] = self.hits[self.induction1ChannelsID[index]]
        return self

    def make_voxels(self):
        self.voxels = []
        ind1HitRef = None
        ind2HitRef = None
        for collChannel, collHits in self.hits.items():
            if len(collHits) == 0:
                continue
            elif len(collHits) != 0:
                if collChannel not in self.collectionChannelsID:
                    continue
                for hitCollId, hitColl in collHits.items():
                    hitCollTime = hitColl.hitTime
                    minimumTimeDiff = 999
                    collHitRef = None
                    for indChannel in self.induction1ChannelsID:
                        if len(self.hits[indChannel]) == 0:
                            continue
                        elif len(self.hits[indChannel]) != 0:
                            hitsInd1 = self.hits[indChannel]
                            for hitInd1Id, hitInd1 in hitsInd1.items():
                                hitIndTime = hitInd1.hitTime
                                if abs(hitCollTime-hitIndTime) < minimumTimeDiff:
                                    minimumTimeDiff = abs(hitCollTime-hitIndTime)
                                    if minimumTimeDiff < self.voxelSeparation:
                                        collHitRef, ind1HitRef = hitColl, hitInd1
                    minimumTimeDiff = 999
                    for indChannel in self.induction2ChannelsID:
                        if len(self.hits[indChannel]) == 0:
                            continue
                        elif len(self.hits[indChannel]) != 0:
                            hitsInd2 = self.hits[indChannel]
                            for hitInd2Id, hitInd2 in hitsInd2.items():
                                hitIndTime = hitInd2.hitTime
                                if abs(hitCollTime-hitIndTime) < minimumTimeDiff:
                                    minimumTimeDiff = abs(hitCollTime-hitIndTime)
                                    if minimumTimeDiff < self.voxelSeparation:
                                        collHitRef, ind2HitRef = hitColl, hitInd2
                    if (collHitRef is not None and ind1HitRef is not None) or (collHitRef is not None and ind2HitRef is not None):
                        self.voxels.append(Voxel(collHitRef, ind1HitRef, ind2HitRef))
        return self

    def make_clusters(self):
        self.clusters = None
        if len(self.voxels) == 0:
            self.clusters = {}
        else:
            eps = self.settings["cluster"]["eps"]
            min_samples = self.settings["cluster"]["minNumberOfVoxels"]
            voxel_coordinates = []
            voxels = []
            for vox in self.voxels:
                if vox.xCord is not None and vox.yCord is not None and vox.zCord is not None:
                    # voxel_coordinates.append((vox.xCord, vox.yCord, vox.zCord))
                    voxel_coordinates.append((vox.xCord, vox.zCord))
                    voxels.append(vox)
            if len(voxels) > 0:
                dist_matrix = pairwise_distances(voxel_coordinates, metric="euclidean")
                dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric="euclidean")
                cluster_labels = dbscan.fit_predict(dist_matrix)
                clusters = {i: [] for i in set(cluster_labels)}
                for i, voxel in enumerate(voxels):
                    cluster_id = cluster_labels[i]
                    clusters[cluster_id].append(voxel)
                self.clusters = {i: [] for i in set(cluster_labels)}
                for clusterId, clusterVoxels in clusters.items():
                    self.clusters[clusterId] = Cluster(clusterVoxels, clusterId)
        if self.clusters is not None:
            self.nClusters = len(self.clusters)
        else:
            self.nClusters = None
        return self


