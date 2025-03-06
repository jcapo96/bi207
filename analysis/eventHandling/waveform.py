import numpy as np
from scipy.signal import find_peaks
import json

class Waveform():
    def __init__(self, waveform):
        self.settings = json.load(open("/Users/jcapo/cernbox/DUNE-IFIC/Software/bi207/analysis/eventHandling/settings.json"))["waveform"]
        self.timeTick = self.settings["timeTick"]
        self.waveform = np.array(waveform)
        self.baseline = np.mean(self.waveform)
        self.nPoints = len(self.waveform)
        self.xPoints = range(self.nPoints)
        self.threshdoldADC = self.settings["thresholdADC"]
        self.threshold = self.baseline + self.threshdoldADC
        self.initialize_waveform()

    def initialize_waveform(self):
        self.nPoints = len(self.waveform)
        #self.normalize_waveform()
        self.find_peaks()
        self.is_bipolar()
        return self.nPoints

    def normalize_waveform(self):
        self.waveform = self.waveform/(np.mean(self.waveform))
        return self.waveform

    def find_peaks(self):
        self.peaks, _ = find_peaks(self.waveform, height=self.threshold)
        self.nPeaks = len(self.peaks)
        return self.peaks, self.nPeaks

    def is_bipolar(self):
        if min(self.waveform) < self.baseline - self.threshdoldADC:
            self.bipolar = True
        else:
            self.bipolar = False
        return self.bipolar
