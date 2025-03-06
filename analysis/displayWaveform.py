import ROOT
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')
import matplotlib.pyplot as plt
from eventHandling.event import Event

import ROOT
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')
import matplotlib.pyplot as plt
from eventHandling.event import Event

pathToData = "/eos/user/j/jcapotor/bi207Data/"
fName = "20230623.root"
rootFile = ROOT.TFile(f"{pathToData}{fName}", "READ")
for event in rootFile.Get("events"):
    event0 = Event(event)
    print(f"EventId: {event0.eventID}")
    if event0.eventID > 1:
        break
    for channel, hit in event0.hits.items():
        if len(hit) != 0:
            print(channel)
            plt.plot(hit[0].Peak.waveform.waveform)
            plt.savefig("displayedWaveform.png")
