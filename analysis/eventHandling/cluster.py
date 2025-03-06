import matplotlib.pyplot as plt
import matplotlib.cm as cm
cmap = cm.get_cmap('viridis')
import numpy as np
import math

class Cluster():
    def __init__(self, voxels, clusterId):
        self.voxels = voxels
        self.clusterId = clusterId
        self.nVoxels = len(self.voxels)
        self.get_crossed_channels()
        self.get_track_time()
        self.find_cluster_hits()

    def find_cluster_hits(self):
        self.clusterHits = []
        for voxel in self.voxels:
            if voxel.hitColl is not None:
                self.clusterHits.append(voxel.hitColl)
            if voxel.hitInd1 is not None:
                self.clusterHits.append(voxel.hitInd1)
            if voxel.hitInd2 is not None:
                self.clusterHits.append(voxel.hitInd2)
        return self


    def get_crossed_channels(self):
        self.deltaC, self.deltaI1 = 0, 0
        self.collChannelsCrossed, self.ind1ChannelsCrossed = [], []
        for voxel in self.voxels:
            try:
                self.collChannelsCrossed.append(voxel.hitColl.channelNumber)
                self.ind1ChannelsCrossed.append(voxel.hitInd1.channelNumber)
            except:
                continue
        try:
            self.deltaC = abs(max(self.collChannelsCrossed) - min(self.collChannelsCrossed))
            self.deltaI1 = max(self.ind1ChannelsCrossed) - min(self.ind1ChannelsCrossed)
        except: pass
        return self

    def find_first_voxel(self):
        minVoxelTime = 999
        self.firstVoxel = None
        for voxel in self.voxels:
            if abs(voxel.voxelTime < minVoxelTime):
                self.firstVoxel = voxel
                minVoxelTime = voxel.voxelTime
        return self

    def find_last_voxel(self):
        maxVoxelTime = 0
        self.lastVoxel = None
        for voxel in self.voxels:
            if abs(voxel.voxelTime > maxVoxelTime):
                self.lastVoxel = voxel
                maxVoxelTime = voxel.voxelTime
        return self

    def get_strip_quantities(self):
        self.deltaC = self.lastVoxel.hitColl.channelNumber - self.firstVoxel.hitColl.channelNumber
        self.deltaI = self.lastVoxel.hitInd1.channelNumber - self.firstVoxel.hitInd1.channelNumber
        return self

    def get_track_time(self):
        self.find_first_voxel()
        self.find_last_voxel()
        self.deltaT = self.lastVoxel.voxelTime - self.firstVoxel.voxelTime
        return self

    def find_cluster_angle(self):
        self.find_first_voxel()
        self.find_last_voxel()
        self.deltaX = self.lastVoxel.xCord - self.firstVoxel.xCord
        self.deltaY = self.lastVoxel.yCord - self.firstVoxel.yCord
        self.deltaZ = self.lastVoxel.zCord - self.firstVoxel.zCord
        self.theta = math.atan(self.deltaY/self.deltaX)*180/np.pi
        self.phi = math.atan(self.deltaZ/self.deltaX)*180/np.pi
        self.psi = math.atan(self.deltaZ/self.deltaY)*180/np.pi
        return self

    def find_cluster_centroid(self):
        self.centroid = None
        if self.voxels is not None:
            x, y, z = 0, 0, 0
            for voxel in self.voxels:
                x += voxel.xCord
                y += voxel.yCord
                z += voxel.zCord
            self.centroid = np.array([x/len(self.voxels), y/len(self.voxels), z/len(self.voxels)])
        return self