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
    if event0.eventID > 47:
        break
    plt.imshow(event0.ADC, cmap='viridis', aspect='auto')
    plt.savefig("ADC.png")