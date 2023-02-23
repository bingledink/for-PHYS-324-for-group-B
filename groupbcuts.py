import ROOT
import uproot
import numpy as np
from   array import array
#import matplotlib.pypolt as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input')
args = parser.parse_args()

fileptr = uproot.open(args.input)

#Pull da data

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




#Setting arrays for tings
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

#some function setting
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

#Electron and muons
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
#jets cuts
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


#Lepton cuts
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






#append all the arrays
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


print(len(l_pt))
print(counter)
#print(j_pt)



