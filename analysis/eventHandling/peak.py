import numpy as np
from eventHandling.waveform import Waveform

class Peak:
    def __init__(self, waveform, peak):
        self.threshold = 1.00
        self.waveform = Waveform(waveform)
        self.peak = peak
        self.bipolar = self.waveform.bipolar
        self.tini = self.waveform.nPoints
        self.tend = 0
        self.width_value = None
        self.integral_value = None
        self.height_value = None
        self.SNR_value = None
        self.initialize()

    def initialize(self):
        self.find_peak_range()
        self.find_peak_center()
        self.calculate_integral()
        self.calculate_height()
        self.calculate_width()
        self.calculate_SNR()
        #self.remove_coherent_noise()

    def find_peak_range(self):
        for timetick in range(0, self.peak):
            if self.waveform.waveform[self.peak - timetick] < self.waveform.baseline * self.threshold:
                self.tini = self.peak - timetick
                break
        if self.bipolar == False:
            for timetick in range(0, self.waveform.nPoints - self.peak):
                if self.waveform.waveform[self.peak + timetick] < self.waveform.baseline * self.threshold:
                    self.tend = self.peak + timetick
                    break
        if self.bipolar == True:
            nBaselineCrosses = 0
            for timetick in range(0, self.waveform.nPoints - self.peak):
                if nBaselineCrosses == 0:
                    if self.waveform.waveform[self.peak + timetick] < self.waveform.baseline * self.threshold:
                        nBaselineCrosses += 1
                        continue
                if nBaselineCrosses > 0:
                    if self.waveform.waveform[self.peak + timetick] > self.waveform.baseline * self.threshold:
                        self.tend = self.peak + timetick
                        break
        return self.tini, self.tend

    def calculate_integral(self):
        if self.tini is None or self.tend is None:
            self.find_peak_range()
        if self.width_value is None:
            self.calculate_width()
        if self.bipolar == False:
            integral_range = abs(self.waveform.waveform[self.tini:self.tend] - self.waveform.baseline)
            self.integral_value = np.trapz(integral_range) * self.waveform.timeTick
        if self.bipolar == True:
            integral_range1 = abs(self.waveform.waveform[self.tini:self.peak_center] - self.waveform.baseline)
            integral_value1 = np.trapz(integral_range1) * self.waveform.timeTick
            integral_range2 = abs(self.waveform.waveform[self.peak_center:self.tend] - self.waveform.baseline)
            integral_value2 = np.trapz(integral_range2) * self.waveform.timeTick
            self.integral_value = integral_value1 + integral_value2
        return self.integral_value

    def calculate_width(self):
        if self.tini is None or self.tend is None:
            self.find_peak_range()
        self.width_value = self.tend - self.tini
        return self.width_value

    def calculate_height(self):
        self.height_value = self.waveform.waveform[self.peak]
        return self.height_value

    def calculate_SNR(self):
        if self.height_value is None:
            self.calculate_height()
        self.SNR_value = self.height_value / self.waveform.baseline
        return self.SNR_value

    def find_peak_center(self):
        self.peak_center = int(np.round(np.mean([self.tini, self.tend])))
        return self.peak_center

    def remove_coherent_noise(self):
        self.waveform.waveform[:self.tini] = [self.waveform.baseline for element in range(0, self.tini)]
        self.waveform.waveform[self.tend:] = [self.waveform.baseline for element in range(self.tend, self.waveform.nPoints)]
        return self.waveform