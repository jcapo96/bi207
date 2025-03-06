import time
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
from array import array
import sys
sys.path.append('./OopFramework/')
sys.path.append('./OopFramework/eventHandling')

from eventHandling.event import Event

class Equalisation():
    def __init__(self):
        # self.pathToRawData = "/eos/project/f/flic-bi207/bi207/DATA/20220511/jsonData/"
        self.pathToRawData = "/eos/user/j/jcapotor/bi207Data/"
        self.pathToSaveSelection = "/"
        self.rootSaveName = "equalisation_2023.root"

    def line(self, x, a, deg):
        return a + np.tan(deg*np.pi/180)*x

    def fill_selection(self, date):
        startTime = time.time()
        rootFile = ROOT.TFile(f"{self.pathToRawData}{date}.root", "READ")
        print(f"Total number of entries: {rootFile.Get('events').GetEntries()}")
        self.saveFile = ROOT.TFile(self.rootSaveName, "RECREATE")

        ### Create arrays that will contain histogram variables
        rootTree = ROOT.TTree(f"selection_{date}", "Selected events and results from selection")
        AccRootTree = ROOT.TTree(f"selection", "Selected events and results from selection")
        deltaT = np.array([0.0], dtype=np.float32)
        deltaC = np.array([0.0], dtype=np.float32)
        deltaI = np.array([0.0], dtype=np.float32)
        deltaL = np.array([0.0], dtype=np.float32)
        deltaZ = np.array([0.0], dtype=np.float32)
        theta = np.array([0.0], dtype=np.float32)
        phi = np.array([0.0], dtype=np.float32)
        psi = np.array([0.0], dtype=np.float32)
        eventId = np.array([0], dtype=np.int32)
        chargePerChannel = {f"chn{chan}":np.array([0.0], dtype=np.float32) for chan in range(0,128)}
        # Create and save the histogram as an image
        selection = ROOT.TH2F(f"selectedEvents_{date}", f"Selected Events {date}", 127, 0, 128, 300, 0, 644)
        selection.GetXaxis().SetTitle("Channel Number")
        selection.GetYaxis().SetTitle("Time Tick")

        selectionXY = ROOT.TH2F(f"XY_selectedEvents_{date}", f"Selected Events {date}", 47, 0, 240, 40, 0,  33*7.5/np.sin(60*np.pi/180))
        selectionXY.GetXaxis().SetTitle("X axis (mm)")
        selectionXY.GetYaxis().SetTitle("Y axis (mm)")

        selectionXZ = ROOT.TH2F(f"XZ_selectedEvents_{date}", f"Selected Events {date}", 47, 0, 240, 300, 0, 644*0.5*1.55)
        selectionXZ.GetXaxis().SetTitle("X axis (mm)")
        selectionXZ.GetYaxis().SetTitle("Z axis (mm)")

        selectionYZ = ROOT.TH2F(f"YZ_selectedEvents_{date}", f"Selected Events {date}", 40, 0, 33*7.5/np.sin(60*np.pi/180), 300, 0, 644*0.5*1.55)
        selectionYZ.GetXaxis().SetTitle("Y axis (mm)")
        selectionYZ.GetYaxis().SetTitle("Z axis (mm)")

        accSelectionXY = ROOT.TH2F(f"XY_selectedEvents", f"Selected Events", 47, 0, 240, 40, 0,  33*7.5/np.sin(60*np.pi/180))
        accSelectionXY.GetXaxis().SetTitle("X axis (mm)")
        selectionXY.GetYaxis().SetTitle("Y axis (mm)")

        accSelectionXZ = ROOT.TH2F(f"XZ_selectedEvents", f"Selected Events", 47, 0, 240, 300, 0, 644*0.5*1.55)
        accSelectionXZ.GetXaxis().SetTitle("X axis (mm)")
        accSelectionXZ.GetYaxis().SetTitle("Z axis (mm)")

        accSelectionYZ = ROOT.TH2F(f"YZ_selectedEvents", f"Selected Events", 40, 0, 33*7.5/np.sin(60*np.pi/180), 300, 0, 644*0.5*1.55)
        accSelectionYZ.GetXaxis().SetTitle("Y axis (mm)")
        accSelectionYZ.GetYaxis().SetTitle("Z axis (mm)")

        # Initialize branches according to specified array types
        rootTree.Branch("deltaT", deltaT, f"deltaT/F")
        rootTree.Branch("deltaC", deltaC, f"deltaC/F")
        rootTree.Branch("deltaI", deltaI, f"deltaI/F")
        rootTree.Branch("deltaL", deltaL, f"deltaL/F")
        rootTree.Branch("deltaZ", deltaZ, f"deltaZ/F")
        rootTree.Branch("theta", theta, f"theta/F")
        rootTree.Branch("phi", phi, f"phi/F")
        rootTree.Branch("psi", psi, f"psi/F")
        rootTree.Branch("eventId", eventId, f"eventId/I")

        AccRootTree.Branch("deltaT", deltaT, f"deltaT/F")
        AccRootTree.Branch("deltaC", deltaC, f"deltaC/F")
        AccRootTree.Branch("deltaI", deltaI, f"deltaI/F")
        AccRootTree.Branch("deltaL", deltaL, f"deltaL/F")
        AccRootTree.Branch("deltaZ", deltaZ, f"deltaZ/F")
        AccRootTree.Branch("theta", theta, f"theta/F")
        AccRootTree.Branch("phi", phi, f"phi/F")
        AccRootTree.Branch("psi", psi, f"psi/F")
        AccRootTree.Branch("eventId", eventId, f"eventId/I")
        for chan, charge in chargePerChannel.items():
            rootTree.Branch(f"charge_{chan}", charge, f"charge_{chan}/F")
            AccRootTree.Branch(f"charge_{chan}", charge, f"charge_{chan}/F")

        ### Starts the event loop where branches are filled
        for event in rootFile.Get("events"):
            event0 = Event(event)
            if (event0.eventID) % 1000 == True:
                print(f"Current event {event0.eventID} -> Elapsed time {time.time() - startTime} s")
            event0.make_voxels()
            event0.make_clusters()
            if event0.clusters is None:
                continue
            for clusterId, cluster in event0.clusters.items():
                if clusterId == -1:
                    continue
                if cluster.nVoxels < 45:
                    continue
                if abs(cluster.deltaC) < 46:
                    continue
                print(f"Found event {event0.eventID}")
                eventId[0] = event0.eventID
                cluster.get_strip_quantities()
                deltaT[0] = cluster.deltaT
                deltaC[0] = abs(cluster.deltaC)
                deltaI[0] = cluster.deltaI
                deltaL[0] = np.sqrt((deltaT[0]*(1.55))**2 + (4/3)*((deltaC[0]*5)**2 + (deltaI[0]*7.5)**2 - (deltaC[0]*5*7.5*deltaI[0])))
                deltaZ[0] = deltaL[0]/abs(deltaC[0])
                hits = cluster.clusterHits
                for hit in hits:
                    channelNumber = hit.channelNumber
                    hitTime = hit.hitTime
                    hitCharge = hit.charge
                    selection.Fill(channelNumber, hitTime, hitCharge/deltaZ[0])
                try:
                    cluster.find_cluster_angle()
                    theta[0] = cluster.theta
                    phi[0] = cluster.phi
                    psi[0] = cluster.psi
                except:
                    print("Not able to fill the cluster")
                    continue
                for voxel in cluster.voxels:
                    charge = voxel.charge
                    xCord = voxel.xCord
                    yCord = voxel.yCord
                    zCord = voxel.zCord
                    chan = voxel.hitColl.channel
                    chargePerChannel[chan][0] = charge/deltaZ[0]
                    selectionXY.Fill(xCord, yCord, charge)
                    selectionXZ.Fill(xCord, zCord, charge)
                    selectionYZ.Fill(yCord, zCord, charge)
                    accSelectionXY.Fill(xCord, yCord, charge)
                    accSelectionXZ.Fill(xCord, zCord, charge)
                    accSelectionYZ.Fill(yCord, zCord, charge)
                rootTree.Fill()
                AccRootTree.Fill()

        # Save the histogram in the ROOT file
        rootFile.Close()

        self.saveFile.cd()
        selectionXY.Write()
        selectionXZ.Write()
        selectionYZ.Write()
        accSelectionXY.Write()
        accSelectionXZ.Write()
        accSelectionYZ.Write()
        selection.Write()
        rootTree.Write()
        AccRootTree.Write()
        self.saveFile.Close()

        endTime = time.time()
        print(f"Elapsed time: {endTime-startTime}")

    def add_data(self, date):
        startTime = time.time()
        rootFile = ROOT.TFile(f"{self.pathToRawData}{date}.root", "READ")
        print(f"Total number of entries: {rootFile.Get('events').GetEntries()}")
        self.saveFile = ROOT.TFile("equalisation.root", "UPDATE")
        rootTree = ROOT.TTree(f"selection_{date}", "Selected events and results from selection")
        accRootTree = self.saveFile.Get("selection")
        # Create and save the histogram as an image
        selectionXY = ROOT.TH2F(f"XY_selectedEvents_{date}", f"Selected Events {date}", 47, 0, 240, 40, 0,  33*7.5/np.sin(60*np.pi/180))
        selectionXY.GetXaxis().SetTitle("X axis (mm)")
        selectionXY.GetYaxis().SetTitle("Y axis (mm)")

        selectionXZ = ROOT.TH2F(f"XZ_selectedEvents_{date}", f"Selected Events {date}", 47, 0, 240, 490, 0, 644*0.5*1.55)
        selectionXZ.GetXaxis().SetTitle("X axis (mm)")
        selectionXZ.GetYaxis().SetTitle("Z axis (mm)")

        selectionYZ = ROOT.TH2F(f"YZ_selectedEvents_{date}", f"Selected Events {date}", 40, 0, 33*7.5/np.sin(60*np.pi/180), 590, 0, 644*0.5*1.55)
        selectionYZ.GetXaxis().SetTitle("Y axis (mm)")
        selectionYZ.GetYaxis().SetTitle("Z axis (mm)")

        accSelectionXY = self.saveFile.Get("XY_selectedEvents")
        accSelectionXZ = self.saveFile.Get("XZ_selectedEvents")
        accSelectionYZ = self.saveFile.Get("YZ_selectedEvents")

        deltaT = np.array([0.0], dtype=np.float32)
        deltaC = np.array([0.0], dtype=np.float32)
        deltaI = np.array([0.0], dtype=np.float32)
        deltaL = np.array([0.0], dtype=np.float32)
        deltaZ = np.array([0.0], dtype=np.float32)
        theta = np.array([0.0], dtype=np.float32)
        phi = np.array([0.0], dtype=np.float32)
        psi = np.array([0.0], dtype=np.float32)
        eventId = np.array([0], dtype=np.int32)

        rootTree.Branch("deltaT", deltaT, f"deltaT/F")
        rootTree.Branch("deltaC", deltaC, f"deltaC/F")
        rootTree.Branch("deltaI", deltaI, f"deltaI/F")
        rootTree.Branch("deltaL", deltaL, f"deltaL/F")
        rootTree.Branch("deltaZ", deltaZ, f"deltaZ/F")
        rootTree.Branch("theta", theta, f"theta/F")
        rootTree.Branch("phi", phi, f"phi/F")
        rootTree.Branch("psi", psi, f"psi/F")
        rootTree.Branch("eventId", eventId, f"eventId/I")
        chargePerChannel = {f"chn{chan}":np.array([0.0], dtype=np.float32) for chan in range(0,128)}
        for chan, charge in chargePerChannel.items():
            rootTree.Branch(f"charge_{chan}", chargePerChannel[chan], f"charge_{chan}/F")
            accRootTree.SetBranchAddress(f"charge_{chan}", chargePerChannel[chan])

        # Set the branch addresses to the numpy arrays
        accRootTree.SetBranchAddress("deltaT", deltaT)
        accRootTree.SetBranchAddress("deltaC", deltaC)
        accRootTree.SetBranchAddress("deltaI", deltaI)
        accRootTree.SetBranchAddress("deltaL", deltaL)
        accRootTree.SetBranchAddress("deltaZ", deltaZ)
        accRootTree.SetBranchAddress("theta", theta)
        accRootTree.SetBranchAddress("phi", phi)
        accRootTree.SetBranchAddress("psi", psi)
        accRootTree.SetBranchAddress("eventId", eventId)
        for event in rootFile.Get("events"):
            event0 = Event(event)
            # if (event0.eventID) > 100:
            #     break
            if (event0.eventID) % 1000 == True:
                print(f"Current event {event0.eventID} -> Elapsed time {time.time() - startTime} s")
            event0.make_voxels()
            event0.make_clusters()
            if event0.clusters is None:
                continue
            for clusterId, cluster in event0.clusters.items():
                if clusterId == -1:
                    continue
                if cluster.nVoxels < 45:
                    continue
                if abs(cluster.deltaC) < 46:
                    continue
                eventId[0] = event0.eventID
                cluster.get_strip_quantities()
                deltaT[0] = cluster.deltaT
                deltaC[0] = abs(cluster.deltaC)
                deltaI[0] = cluster.deltaI
                deltaL[0] = np.sqrt((deltaT[0]*(1.55))**2 + (4/3)*((deltaC[0]*5)**2 + (deltaI[0]*7.5)**2 - (deltaC[0]*5*7.5*deltaI[0])))
                deltaZ[0] = deltaL[0]/abs(deltaC[0])
                try:
                    cluster.find_cluster_angle()
                    theta[0] = cluster.theta
                    phi[0] = cluster.phi
                    psi[0] = cluster.psi
                except:
                    continue
                for voxel in cluster.voxels:
                    charge = voxel.charge
                    xCord = voxel.xCord
                    yCord = voxel.yCord
                    zCord = voxel.zCord
                    chan = voxel.hitColl.channel
                    chargePerChannel[chan][0] = charge/deltaZ[0]
                    #chargePerChannelBranch = accRootTree.GetBranch(f"charge_{chan}")
                    #chargePerChannelBranch.SetAddress(np.array([chargePerChannel[chan][0]], dtype=np.float32))
                    selectionXY.Fill(xCord, yCord, charge)
                    selectionXZ.Fill(xCord, zCord, charge)
                    selectionYZ.Fill(yCord, zCord, charge)
                    accSelectionXY.Fill(xCord, yCord, charge)
                    accSelectionXZ.Fill(xCord, zCord, charge)
                    accSelectionYZ.Fill(yCord, zCord, charge)
                rootTree.Fill()
                accRootTree.Fill()

        # Save the histogram in the ROOT file
        selectionXY.Write()
        selectionXZ.Write()
        selectionYZ.Write()
        accSelectionXY.Write()
        accSelectionXZ.Write()
        accSelectionYZ.Write()
        rootTree.Write()
        accRootTree.Write()
        self.saveFile.Close()
        rootFile.Close()
        endTime = time.time()
        print(f"Elapsed time: {endTime-startTime}")

    def do_equalisation(self):
        rootFile = ROOT.TFile("/afs/cern.ch/user/j/jcapotor/bi207/equalisation_2023.root", "UPDATE")
        collChannels = [f"chn{index}" for index in range(0,48)]
        listOfEqualisationSets = []
        cnt_selection = 1
        for key in rootFile.GetListOfKeys():
            if "selection" not in key.GetName():
                continue
            if "selection" == key.GetName():
                listOfEqualisationSets.append(f"{key.GetName()};{cnt_selection}")
                cnt_selection += 1
            else:
                listOfEqualisationSets.append(key.GetName())
        print(listOfEqualisationSets)
        for selection in listOfEqualisationSets:
            print(selection)
            rootFile.cd()
            if rootFile.GetListOfKeys().Contains(f"EqualisationPlot_{selection}"):
                rootFile.Delete(f"EqualisationPlot_{selection};*")
            if rootFile.GetListOfKeys().Contains(f"landauFits_{selection}"):
                rootFile.Delete(f"landauFits_{selection};*")
            if ";" in selection:
                selectionName = selection.split(";")
            else:
                selectionName = (selection,"")
            rootFile.mkdir(f"landauFits_{selectionName[0]}_{selectionName[1]}")
            landau_func = ROOT.TF1("landau_func", "landau", 0, 150)
            composite_func = ROOT.TF1("composite_func", "[0]*TMath::Landau(x, [1], [2])*[3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -10, 200)
            composite_func.SetParameters(20.0, 40.0, 2.0, 20, 0.0, 5.0, 20., 50.0, 10.0)
            landau_func.SetParameters(10.0, 40.0, 2.0)
            mpv = np.array([0.0], dtype=np.float32)
            sigma = np.array([0.0], dtype=np.float32)
            rootTree = ROOT.TTree(f"landaus_{selection}", f"Landau fits for Collection Channels - {selection}")
            rootTree.Branch("mpv", mpv, "mpv/F")
            rootTree.Branch("sigma", sigma, "sigma/F")
            mpvs, channel_numbers = [], []
            for nchan, channel in enumerate(collChannels):
                canvas = ROOT.TCanvas("canvas", "Landau-Gaussian Convolution Fit", 800, 600)
                data_hist = ROOT.TH1F("data_hist", "Landau Fit Example", 100, -10, 200)
                for entry in rootFile.Get(f"{selection}"):
                    data_hist.Fill(getattr(entry, f"charge_{channel}"))
                data_hist.Fit(landau_func, "R")
                fit_parameters = landau_func.GetParameters()
                data_hist.Draw()
                landau_func.Draw("SAME")
                rootFile.cd(f"landauFits_{selectionName[0]}_{selectionName[1]}")
                canvas.Write(f"landau_fit_{channel}.png")
                mpv[0] = fit_parameters[1]
                mpvs.append(fit_parameters[1])
                channel_numbers.append(nchan)
                rootTree.Fill()
            rootTree.Write()
            mpvs = array("f", mpvs)
            channel_numbers = array("f", channel_numbers)
            graph = ROOT.TGraph(len(channel_numbers), channel_numbers, mpvs)
            graph.SetMarkerStyle(20)
            graph.SetTitle(rf"Equalisation {selection}: Collection channels -> mu={0:.2f}".format(np.mean(mpvs[23:-2]/np.mean(mpvs[1:23]))))
            graph.GetXaxis().SetTitle("Channel Number")
            graph.GetYaxis().SetTitle("dQ/dx (ADC/mm)")
            rootFile.cd()
            graph.Write(f"EqualisationPlot_{selection}")
        rootFile.Close()

# Equalisation().fill_selection("20230722")
# Equalisation().add_data("20220512")
# Equalisation().add_data("20220513")
# Equalisation().add_data("20220514")
# Equalisation().add_data("20220515")
Equalisation().do_equalisation()