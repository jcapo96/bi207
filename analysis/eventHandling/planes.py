import numpy as np
import sympy as sp
from sympy import symbols, lambdify
import matplotlib.pyplot as plt
class Strip():
    def __init__(self, channel):
        self.channel = channel
        self.nCollChannels = 48
        self.nInd1Channels = 40
        self.nInd2Channels = 40
        self.collPitch = 5
        self.indPitch = 7.5
        self.step = 5E-1
        self.xLength = self.collPitch*self.nCollChannels
        self.yLength = self.indPitch*(self.nInd1Channels-8)
        self.collectionChannelsID = ["chn"+str(i) for i in range(0, self.nCollChannels)]
        self.induction1ChannelsID = ["chn"+str(i+self.nCollChannels) for i in range(0, self.nInd1Channels)]
        self.induction2ChannelsID = ["chn"+str(i+self.nCollChannels+self.nInd1Channels) for i in range(0, self.nInd2Channels)]
        self.classify_channel()
        self.set_channel_angle()
        self.set_channel_equation()

    def classify_channel(self):
        if self.channel in self.collectionChannelsID:
            self.channelType = "coll"
        elif self.channel in self.induction1ChannelsID:
            self.channelType = "ind1"
        elif self.channel in self.induction2ChannelsID:
            self.channelType = "ind2"
        else:
            print(self.channel)
        return self

    def set_channel_angle(self):
        if self.channelType == "coll":
            self.channelAngle = 0
        if self.channelType == "ind1":
            self.channelAngle = 30*np.pi/180
        if self.channelType == "ind2":
            self.channelAngle = -30*np.pi/180
        return self

    def set_channel_equation(self):
        x = sp.symbols('x')
        channelNumber = int(self.channel.split("chn")[1])
        if self.channelType == "coll":
            self.function = channelNumber*self.collPitch
        elif self.channelType == "ind1":
            channelNumber = channelNumber - 48
            m = np.tan(self.channelAngle)
            x0 = (channelNumber-31)*(7.5/np.sin(abs(self.channelAngle)))
            self.function = m * (x - x0)
        elif self.channelType == "ind2":
            channelNumber = channelNumber - 88
            m = np.tan(self.channelAngle)
            x0 = (8+channelNumber)*(7.5/np.sin(abs(self.channelAngle)))
            self.function = m * (x - x0)
        return self

class Mesh():
    def __init__(self):
        self.nCollChannels = 48
        self.nInd1Channels = 40
        self.nInd2Channels = 40
        self.xCord = None
        self.yCord = None
        self.xCordErr = 1.5
        self.yCordErr = 1.5
        self.channels = ["chn"+str(i) for i in range(128)]
        self.collectionChannelsID = ["chn"+str(i) for i in range(0, self.nCollChannels)]
        self.induction1ChannelsID = ["chn"+str(i+self.nCollChannels) for i in range(0, self.nInd1Channels)]
        self.induction2ChannelsID = ["chn"+str(i+self.nCollChannels+self.nInd1Channels) for i in range(0, self.nInd2Channels)]

    def draw_mesh(self):
        fig, axes = plt.subplots(1, 1)
        channels = [f"chn{chan}" for chan in range(0, 128)]
        for chan in channels:
            strip = Strip(chan)
            if (type(strip.function)) == int:
                axes.axvline(strip.function, linewidth=1.0)
            else:
                x = symbols('x')
                xValues = np.linspace(0, 47*5, 1000)
                func = lambdify(x, strip.function, 'numpy')
                if chan in self.induction1ChannelsID:
                    axes.plot(xValues, func(xValues), linewidth=1.0, label=chan)
                else:
                    axes.plot(xValues, func(xValues), linewidth=1.0, color="black")
        axes.set_xlim(0, 48*5)
        axes.set_ylim(0, 40*7.5)
        axes.legend(ncol=2, loc='upper left', fontsize=6)
        fig.savefig("mesh.png", dpi=300)

    def find_intersection(self, hit1=None, hit2=None, hit3=None):
        x = sp.symbols('x')
        wires = {}
        self.xCord, self.yCord = None, None
        if hit1 is not None:
            wires[hit1.channel] = Strip(hit1.channel)
        if hit2 is not None:
            wires[hit2.channel] = Strip(hit2.channel)
        if hit3 is not None:
            wires[hit3.channel] = Strip(hit3.channel)

        if len(wires) <= 1:
            pass
            #print("You can not find the intersection")
        elif len(wires) == 2:
            for channel, wire in wires.items():
                if wire.channelAngle == 0:
                    self.xCord = wire.function
                    collChannel = channel
                    break

            if self.xCord is None:
                self.xCord = sp.solve(wires[list(wires.keys())[0]].function - wires[list(wires.keys())[1]].function, x)[0]
                self.yCord = wires[list(wires.keys())[0]].function.subs(x, self.xCord)
            else:
                wires.pop(collChannel)
                self.yCord = wires[list(wires.keys())[0]].function.subs(x, self.xCord)
                if self.yCord < 0 or self.yCord > 32*7.5/np.sin(abs(wires[list(wires.keys())[0]].channelAngle)):
                    self.xCord = None
                    self.yCord = None
                    # print(f"Length = 2: The intersection does not exists between channels {[collChannel] + list(wires.keys())}")
                else:
                    self.xCord = round(self.xCord, 1)
                    self.yCord = round(self.yCord, 1)
        elif len(wires) == 3:
            collChannel = list(wires.keys())[wires.values() == 0]
            collXCord = wires[collChannel].function
            self.xCord = round(collXCord, 0)
            yCord = []
            yCord.append(round(wires[list(wires.keys())[1]].function.subs(x, self.xCord), 1))
            yCord.append(round(wires[list(wires.keys())[2]].function.subs(x, self.xCord), 1))
            yCord = np.array(yCord, dtype=float)
            self.yCordErr = round(abs(yCord[0] - yCord[1])/2,1)
            self.yCord = round(np.mean(yCord),1)
            # self.yCord = round((yCord[0]),1)

        self.wires = wires
        return self