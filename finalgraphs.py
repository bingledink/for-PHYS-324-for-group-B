import uproot
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
from array import array

plt.style.use([hep.style.ROOT, hep.style.firamath])

#all elec and muon variables
fileptr = uproot.open("TT_Dilept_13.root")



e_pt = fileptr["Delphes_Ntuples"]['elec_pt'].array()
e_eta = fileptr["Delphes_Ntuples"]['elec_eta'].array()
e_phi = fileptr["Delphes_Ntuples"]['elec_phi'].array()
e_pt_flat = []
e_eta_flat = []
e_phi_flat = []
for evtidx in range(len(e_pt)):
  for i in range(len(e_pt[evtidx])):
    e_pt_flat.append(e_pt[evtidx][i])
    e_eta_flat.append(e_eta[evtidx][i])
    e_phi_flat.append(e_phi[evtidx][i])
  #for i in range(len(mu_pt[evtidx])):

plt.hist(e_eta_flat, density = True, bins = 300,)
plt.xlim(-2.4,2.4)
plt.title('Elec Eta')
plt.ylabel('Frequency')
plt.xlabel('Eta (GeV)')
'''
muon_eta = fileptr["CutTree"]['muon_eta'].array()
plt.hist(muon_eta, density = True, bins = 300, log = True)
plt.ylabel('Frequency')
plt.xlabel('Transverse Momentum (GeV)')

muon_phi = fileptr["CutTree"]['muon_phi'].array()
plt.hist(muon_phi, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

elec_pt = fileptr["CutTree"]['elec_pt'].array()
plt.hist(elec_pt, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

elec_eta = fileptr["CutTree"]['elec_eta'].array()
plt.hist(elec_eta, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

elec_phi = fileptr["CutTree"]['elec_phi'].array()
plt.hist(elec_phi, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

jet_pt = fileptr["CutTree"]['jet_pt'].array()
plt.hist(jet_pt, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

jet_eta = fileptr["CutTree"]['jet_eta'].array()
plt.hist(jet_eta, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')

jet_phi = fileptr["CutTree"]['jet_phi'].array()
plt.hist(jet_phi, density = True, bins = 300, log = True)
plt.ylabel('Momentum')
plt.xlabel('Energy Level ')
'''
'''
muon_eta = fileptr["CutTree"]['muon_eta'].array()
plt.hist(muon_eta, density = True, bins = 300)

muon_phi = fileptr["CutTree"]['muon_phi'].array()
plt.hist(muon_phi, density = True, bins = 300)

muon_charge = fileptr["CutTree"]['muon_charge'].array()
plt.hist(muon_charge, density = True, bins = 300)


elec_pt = fileptr["CutTree"]['elec_pt'].array()
plt.hist(elec_pt, density = True, bins = 300)

elec_eta = fileptr["CutTree"]['elec_eta'].array()
plt.hist(elec_eta, density = True, bins = 300)

elec_phi = fileptr["CutTree"]['elec_phi'].array()
plt.hist(elec_phi, density = True, bins = 300)

elec_charge = fileptr["CutTree"]['elec_charge'].array()
plt.hist(elec_charge, density = True, bins = 300)
'''




plt.show()