from eventHandling.planes import Mesh
import numpy as np
import json

class Voxel():
    def __init__(self, hitColl, hitInd1, hitInd2):
        self.settings =self.settings = json.load(open("/Users/jcapo/cernbox/DUNE-IFIC/Software/bi207/analysis/eventHandling/settings.json"))
        self.driftVelocity = self.settings["voxel"]["driftVelocity"]
        self.timeTick = self.settings["voxel"]["timeTick"]
        self.hitColl = hitColl
        self.hitInd1 = hitInd1
        self.hitInd2 = hitInd2
        self.hits = [self.hitColl, self.hitInd1, self.hitInd2]
        self.xCord = None
        self.yCord = None
        self.zCord = None
        self.xCordErr = self.settings["strip"]["xCordErr"]
        self.yCordErr = self.settings["strip"]["yCordErr"]
        self.zCordErr = self.settings["strip"]["zCordErr"]
        self.voxelTime = self.hitColl.hitTime*self.timeTick
        self.positionVector = {}
        self.make_cords()
        self.get_position_vector()
        self.make_energy()

    def make_cords(self):
        self.mesh = Mesh()
        self.mesh.find_intersection(self.hitColl, self.hitInd1, self.hitInd2)
        self.xCord, self.yCord = self.mesh.xCord, self.mesh.yCord
        self.xCordErr, self.yCordErr = self.mesh.xCordErr, self.mesh.yCordErr
        self.zCord = round(self.timeTick*self.hitColl.hitTime*self.driftVelocity, 1) #drift_velocitiy
        return self

    def get_position_vector(self):
        self.positionVector = np.array([self.xCord, self.yCord, self.zCord], dtype=float)
        self.positionVectorErr = np.array([self.xCordErr, self.yCordErr, self.zCordErr], dtype=float)
        return self

    def make_energy(self):
        self.charge = self.hitColl.charge
        return self

