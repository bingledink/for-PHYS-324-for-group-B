import ROOT
import uproot
import numpy as np
from   array import array

fileptr = uproot.open("TT_Dilept_13.root")

outputfile = ROOT.TFile("NewRoot.root", 'recreate')
tree = ROOT.TTree("CutTree", "CutTree")

#Electron arrays
elec_pt_arr = array('f', [0.])
elec_eta_arr = array('f', [0.])
elec_phi_arr = array('f', [0.])
elec_charge_arr = array('f', [0.])

#Muon arrays
muon_pt_arr = array('f', [0.])
muon_eta_arr = array('f', [0.])
muon_phi_arr = array('f', [0.])
muon_charge_arr = array('f', [0.])

#Jet arrays
ljet_pt_arr = array('f',[.0])
ljet_eta_arr = array('f',[.0])
ljet_phi_arr = array('f',[.0])
ljet_mass_arr = array('f',[.0])

sljet_pt_arr = array('f',[.0])
sljet_eta_arr = array('f',[.0])
sljet_phi_arr = array('f',[.0])
sljet_mass_arr = array('f',[.0])
#Lepton arrays
l_pt_arr = array('f',[.0])
l_eta_arr = array('f',[.0])
l_phi_arr = array('f',[.0])
l_mass_arr = array('f',[.0])

sl_pt_arr = array('f',[.0])
sl_eta_arr = array('f',[.0])
sl_phi_arr = array('f',[.0])
sl_mass_arr = array('f',[.0])

#Elec branches
tree.Branch("elec_pt", elec_pt_arr, "elec_pt/F")
tree.Branch("elec_eta", elec_arr, "elec_eta/F")
tree.Branch("elec_phi", elec_phi_arr, "elec_phi/F")
tree.Branch("elec_charge", elec_charge_arr, "elec_charge/F")

#Muon branches
tree.Branch("muon_pt", muon_pt_arr, "muon_pt/F")
tree.Branch("muon_eta", muon_eta_arr, "muon_eta/F")
tree.Branch("muon_phi", muon_phi_arr, "muon_phi/F")
tree.Branch("muon_charge", muon_charge_arr, "muon_charge/F")

#Jets branches
tree.Branch("ljet_pt", ljet_pt_arr, "ljet_pt/F")
tree.Branch("ljet_eta", ljet_eta_arr, "ljet_eta/F")
tree.Branch("ljet_phi", ljet_phi_arr, "ljet_phi/F")
tree.Branch("ljet_mass", ljet_mass_arr, "ljet_mass/F")
tree.Branch("sljet_pt", sljet_pt_arr, "sljet_pt/F")
tree.Branch("sljet_eta", sljet_eta_arr, "sljet_eta/F")
tree.Branch("sljet_phi", sljet_phi_arr, "sljet_phi/F")
tree.Branch("sljet_mass", sljet_mass_arr, "sljet_mass/F")


#Lepton Branches
tree.Branch("l_pt", l_pt_arr, "l_pt/F")
tree.Branch("l_eta", l_eta_arr, "l_eta/F")
tree.Branch("l_phi", l_phi_arr, "l_phi/F")
tree.Branch("l_mass", l_mass_arr, "l_mass/F")
tree.Branch("sl_pt", sl_pt_arr, "sl_pt/F")
tree.Branch("sl_eta", sl_eta_arr, "sl_eta/F")
tree.Branch("sl_phi", sl_phi_arr, "sl_phi/F")
tree.Branch("sl_mass", sl_mass_arr, "sl_mass/F")


#originaltags = [mu_pt, mu_eta, mu_phi, mu_charge]
#flattags = [mu_pt_flat, mu_eta_flat, mu_phi_flat, mu_charge_flat]

for i in range(len(e_pt)):
    e_pt_arr[0] = e_pt[i]
    e_phi_arr[0] = e_phi[i]
    e_eta_arr[0] = e_eta[i]
    e_charge_arr[0] = e_charge[i]
    
    muon_pt_arr[0] = mu_pt[i]
    muon_eta_arr[0] = mu_eta[i]
    muon_phi_arr[0] = mu_phi[i]
    muon_charge_arr[0] = mu_charge[i]
    
    ljet_pt_arr[0] = ljet_pt[i]
    ljet_eta_arr[0] = ljet_eta[i]
    ljet_phi_arr[0] = ljet_phi[i]
    ljet_mass_arr[0] = ljet_mass[i]
    
    sljet_pt_arr[0] = sljet_pt[i]
    sljet_eta_arr[0] = sljet_eta[i]
    sljet_phi_arr[0] = sljet_phi[i]
    sljet_mass_arr[0] = sljet_mass[i]
    
    l_pt_arr[0] = l_pt[i]
    l_eta_arr[0] = l_eta[i]
    l_phi_arr[0] = l_phi[i]
    l_mass_arr[0] = l_mass[i]
    
    sl_pt_arr[0] = sl_pt[i]
    sl_eta_arr[0] = sl_eta[i]
    sl_phi_arr[0] = sl_phi[i]
    sl_pt_mass[0] = sl_mass[i]
    
    tree.Fill()

outputfile.Write()
outputfile.Close()