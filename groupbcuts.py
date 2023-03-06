import ROOT
import uproot
import numpy as np
from   array import array
#import matplotlib.pypolt as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input')
parser.add_argument('-o', '--output', help='Output')
args = parser.parse_args()

fileptr = uproot.open(args.input)
outputfile = ROOT.TFile(args.output, 'recreate')
tree = ROOT.TTree("CutTree", "CutTree")

elec_pt = fileptr['Delphes_Ntuples']['elec_pt'].array()
elec_eta = fileptr['Delphes_Ntuples']['elec_eta'].array()
elec_charge = fileptr['Delphes_Ntuples']['elec_charge'].array()
elec_phi = fileptr['Delphes_Ntuples']['elec_phi'].array()
elec_mass = fileptr['Delphes_Ntuples']['elec_mass'].array()

muon_pt = fileptr['Delphes_Ntuples']['muon_pt'].array()
muon_eta = fileptr['Delphes_Ntuples']['muon_eta'].array()
muon_reliso = fileptr['Delphes_Ntuples']['muon_reliso'].array()
muon_charge = fileptr['Delphes_Ntuples']['muon_charge'].array()
muon_phi = fileptr['Delphes_Ntuples']['muon_phi'].array()
muon_mass = fileptr['Delphes_Ntuples']['muon_mass'].array()

jet_pt = fileptr['Delphes_Ntuples']['jet_pt'].array()
jet_eta = fileptr['Delphes_Ntuples']['jet_eta'].array()
jet_btag = fileptr['Delphes_Ntuples']['jet_btag'].array()
jet_phi = fileptr['Delphes_Ntuples']['jet_phi'].array()
jet_mass = fileptr['Delphes_Ntuples']['jet_phi'].array()

#new
met_pt = fileptr['Delphes_Ntuples']['met_pt'].array()
met_phi = fileptr['Delphes_Ntuples']['met_phi'].array()

weight = fileptr['Delphes_Ntuples']['weight'].array()
scalar_ht = fileptr['Delphes_Ntuples']['scalar_ht'].array()

genjet_pt = fileptr['Delphes_Ntuples']['genjet_pt'].array()
genjet_eta = fileptr['Delphes_Ntuples']['genjet_eta'].array()
genjet_phi = fileptr['Delphes_Ntuples']['genjet_phi'].array()
genjet_mass = fileptr['Delphes_Ntuples']['genjet_mass'].array()

genpart_pt = fileptr['Delphes_Ntuples']['genpart_pt'].array()
genpart_eta = fileptr['Delphes_Ntuples']['genpart_eta'].array()
genpart_phi = fileptr['Delphes_Ntuples']['genpart_phi'].array()
genpart_mass = fileptr['Delphes_Ntuples']['genpart_mass'].array()
genpart_pid = fileptr['Delphes_Ntuples']['genpart_pid'].array()
genpart_status = fileptr['Delphes_Ntuples']['genpart_status'].array()
genpart_charge = fileptr['Delphes_Ntuples']['genpart_charge'].array()

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
tree.Branch("elec_eta", elec_eta_arr, "elec_eta/F")
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



#array for events that passed cuts
onesnzeroes = [0] * len(elec_pt)

e4vector = ROOT.TLorentzVector()
mu4vector = ROOT.TLorentzVector()

e_pt = []
e_eta = []
e_phi = []
e_charge = []

mu_pt = []
mu_eta = []
mu_phi = []
mu_charge = []

ljet_pt = []
ljet_eta = []
ljet_phi = []
ljet_btag = []
ljet_mass = []

sljet_pt = []
sljet_eta = []
sljet_phi = []
sljet_btag = []
sljet_mass = []


l_pt = []
l_eta = []
l_phi = []
l_mass = []

sl_pt = []
sl_eta = []
sl_phi = []
sl_mass = []


def deltaphi(e_phi, m_phi):
    d_phi = e_phi - m_phi
    if (d_phi > np.pi):
        d_phi -= 2*np.pi
    if (d_phi < -np.pi):
        d_phi += 2*np.pi
    return d_phi


def dR(e_phi, e_eta, m_phi, m_eta):
    d_eta = abs(e_eta - m_eta)
    d_phi = deltaphi(e_phi, m_phi)
    return np.sqrt(d_phi**2 + d_eta**2)

for event_idx in range(len(elec_pt)):
    e_idx = []
    mu_idx = []

    ef_idx = []
    muf_idx = []

    counter = 0
    j_idx = []
    jf_idx = []
    
    for i in range(len(muon_pt[event_idx])):
        if muon_pt[event_idx][i] < 20:
            continue
        if abs(muon_eta[event_idx][i]) > 2.4:
            continue
        if muon_reliso[event_idx][i] > 0.15:
            continue
        mu_idx.append(i)
    
    
    for i in range(len(elec_pt[event_idx])):
        if elec_pt[event_idx][i] < 20:
            continue
        if abs(elec_eta[event_idx][i]) > 2.4 or (1.4442<abs(elec_eta[event_idx][i])<1.5660):
            continue
        e_idx.append(i)
    
    if (len(e_idx) == 0 or len(mu_idx) == 0):
        continue
    
    for i in range(len(e_idx)):
        for j in range(len(mu_idx)):

            tmp_e_idx = e_idx[i]
            tmp_mu_idx = mu_idx[j]

            if (elec_charge[event_idx][tmp_e_idx] * muon_charge[event_idx][tmp_mu_idx] == -1):
                    ef_idx.append(tmp_e_idx)
                    muf_idx.append(tmp_mu_idx)

        # Ensure such a pairing exists
    if (len(ef_idx) == 0 or len(muf_idx) == 0):
        continue

    e_index = ef_idx[0]
    mu_index = muf_idx[0]

    for i in range(len(jet_pt[event_idx])):
        if jet_pt[event_idx][i] < 30:
            continue
        if abs(jet_eta[event_idx][i]) > 2.4 or (1.4442<abs(jet_eta[event_idx][i])<1.5660):
            continue
        if (dR(elec_phi[event_idx][e_index], elec_eta[event_idx][e_index], jet_phi[event_idx][i], jet_eta[event_idx][i]) < 0.4) or dR(muon_phi[event_idx][mu_index],muon_eta[event_idx][mu_index],jet_phi[event_idx][i],jet_eta[event_idx][i]) < 0.4:
            continue

        j_idx.append(i)

    for i in range(len(jet_btag[event_idx])):
        if jet_btag[event_idx][i] > 0:
            counter+=1
    if counter == 0:
        continue

    if len(j_idx) < 2:
        continue



    ljet_idx = j_idx[0]
    sljet_idx = j_idx[1]



    if elec_pt[event_idx][e_index] > muon_pt[event_idx][mu_index] and elec_pt[event_idx][e_index] > 25 :
        l_pt.append(elec_pt[event_idx][e_index])
        sl_pt.append(muon_pt[event_idx][mu_index])
        l_eta.append(elec_eta[event_idx][e_index])
        sl_eta.append(muon_eta[event_idx][mu_index])
        l_phi.append(elec_phi[event_idx][e_index])
        sl_phi.append(muon_phi[event_idx][mu_index])
        l_mass.append(elec_mass[event_idx][e_index])
        sl_mass.append(muon_mass[event_idx][mu_index])

    elif muon_pt[event_idx][mu_index] > elec_pt[event_idx][e_index] and muon_pt[event_idx][mu_index] > 25:
        l_pt.append(muon_pt[event_idx][mu_index])
        sl_pt.append(elec_pt[event_idx][e_index])
        sl_eta.append(elec_eta[event_idx][e_index])
        l_eta.append(muon_eta[event_idx][mu_index])
        sl_phi.append(elec_phi[event_idx][e_index])
        l_phi.append(muon_phi[event_idx][mu_index])
        sl_mass.append(elec_mass[event_idx][e_index])
        l_mass.append(muon_mass[event_idx][mu_index])
    else:
        continue
        
    e4vector.SetPtEtaPhiM((elec_pt[event_idx][e_index]),(elec_eta[event_idx][e_index]),(elec_phi[event_idx][e_index]),(elec_mass[event_idx][e_index]))
    mu4vector.SetPtEtaPhiM((muon_pt[event_idx][mu_index]),(muon_eta[event_idx][mu_index]),(muon_phi[event_idx][mu_index]),(muon_mass[event_idx][mu_index]))
    if (e4vector + mu4vector).M() < 20:
        continue







    e_pt.append(elec_pt[event_idx][e_index])
    e_eta.append(elec_eta[event_idx][e_index])
    e_phi.append(elec_phi[event_idx][e_index])
    e_charge.append(elec_charge[event_idx][e_index])
    
    mu_pt.append(muon_pt[event_idx][mu_index])
    mu_eta.append(muon_eta[event_idx][mu_index])
    mu_phi.append(muon_phi[event_idx][mu_index])
    mu_charge.append(muon_charge[event_idx][mu_index])

    ljet_pt.append(jet_pt[event_idx][ljet_idx])
    ljet_phi.append(jet_phi[event_idx][ljet_idx])
    ljet_eta.append(jet_eta[event_idx][ljet_idx])
    ljet_mass.append(jet_mass[event_idx][ljet_idx])

    sljet_pt.append(jet_pt[event_idx][sljet_idx])
    sljet_phi.append(jet_phi[event_idx][sljet_idx])
    sljet_eta.append(jet_eta[event_idx][sljet_idx])
    sljet_mass.append(jet_mass[event_idx][sljet_idx])
    
    onesnzeroes = [0] * len(elec_pt)

for i in range(len(e_pt)):
    elec_pt_arr[0] = e_pt[i]
    elec_phi_arr[0] = e_phi[i]
    elec_eta_arr[0] = e_eta[i]
    elec_charge_arr[0] = e_charge[i]
    
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
    sl_mass_arr[0] = sl_mass[i]
    
    tree.Fill()

outputfile.Write()
outputfile.Close()

print(len(l_pt))
print(counter)
#print(j_pt)
