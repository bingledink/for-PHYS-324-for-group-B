import ROOT
import uproot
import numpy as np
from   array import array

fileptr = uproot.open("TT_Dilept_13.root")
elec_pt = fileptr['Delphes_Ntuples']['elec_pt'].array()

jet_pt = fileptr['Delphes_Ntuples']['ljet_pt'].array()
jet_phi = fileptr['Delphes_Ntuples']['ljet_phi'].array()
jet_eta = fileptr['Delphes_Ntuples']['sljet_eta'].array()
jet_mass = fileptr['Delphes_Ntuples']['sljet_mass'].array()

outputfile = ROOT.TFile("NewRoot.root", 'recreate')
tree = ROOT.TTree("CutTree", "CutTree")

elec_pt_arr = array('f', [0.])


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


tree.Branch("elec_pt", elec_pt_arr, "elec_pt/F")



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



elec_pt_flat = []

for i in elec_pt:
    for v in i:
        if v > 20:
            elec_pt_flat.append(v)

for i in elec_pt_flat:
    elec_pt_arr[0] = i
    tree.Fill()

outputfile.Write()
outputfile.Close()