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

met_pt_arr = []
met_phi_arr = []

scalar_ht_arr = []
total_pt = []
total_jet_pt = []

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

    e4vector.SetPtEtaPhiM((elec_pt[event_idx][e_index]),(elec_eta[event_idx][e_index]),(elec_phi[event_idx][e_index]),(elec_mass[event_idx][e_index]))
    mu4vector.SetPtEtaPhiM((muon_pt[event_idx][mu_index]),(muon_eta[event_idx][mu_index]),(muon_phi[event_idx][mu_index]),(muon_mass[event_idx][mu_index]))
    if (e4vector + mu4vector).M() < 20:
        continue


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
    
    met_pt_arr.append(met_pt[event_idx][0])
    met_phi_arr.append(met_pt[event_idx][0])
    
    scalar_ht_arr.append(scalar_ht[event_idx][0])
    
    onesnzeroes[event_idx] = 1

    x_total_jet_pt = 0
    for i in range(len(jet_pt[event_idx])):
        x_total_jet_pt += jet_pt[event_idx][i]
    total_jet_pt.append(x_total_jet_pt)
    x_elec_pt = 0
    x_muon_pt = 0
    for i in range(len(elec_pt[event_idx])):
        x_elec_pt += elec_pt[event_idx][i]
    for i in range(len(muon_pt[event_idx])):
        x_muon_pt += muon_pt[event_idx][i]
    x_total_pt = x_elec_pt + x_muon_pt + x_total_jet_pt
    total_pt.append(x_total_pt)
   

#np arrays
onesnzeroes = np.array(onesnzeroes)
'''    
np_weight = np.array(weight)
    
np_jet_pt = np.array(jet_pt)
np_jet_eta = np.array(jet_eta)
np_jet_phi = np.array(jet_phi)
np_jet_mass = np.array(jet_mass)
np_jet_btag = np.array(jet_btag)
    
np_genjet_pt = np.array(genjet_pt)
np_genjet_eta = np.array(genjet_eta)
np_genjet_phi = np.array(genjet_phi)
np_genjet_mass = np.array(genjet_mass)
    
np_genpart_pt = np.array(genpart_pt)
np_genpart_eta = np.array(genpart_eta)
np_genpart_phi = np.array(genpart_phi)
np_genpart_mass = np.array(genpart_mass)
np_genpart_charge = np.array(genpart_charge)
np_genpart_status = np.array(genpart_status)
np_genpart_pid = np.array(genpart_pid)
'''
#np array setting

weight_ones = weight[onesnzeroes == 1]

jet_pt_ones = jet_pt[onesnzeroes == 1]
jet_eta_ones = jet_eta[onesnzeroes == 1]
jet_phi_ones = jet_phi[onesnzeroes == 1]
jet_mass_ones = jet_mass[onesnzeroes == 1]
jet_btag_ones = jet_btag[onesnzeroes == 1]

genjet_pt_ones = genjet_pt[onesnzeroes == 1]
genjet_eta_ones = genjet_eta[onesnzeroes == 1]
genjet_phi_ones = genjet_phi[onesnzeroes == 1]
genjet_mass_ones = genjet_mass[onesnzeroes == 1]

genpart_pt_ones = genpart_pt[onesnzeroes == 1]
genpart_eta_ones = genpart_eta[onesnzeroes == 1]
genpart_phi_ones = genpart_phi[onesnzeroes == 1]
genpart_mass_ones = genpart_mass[onesnzeroes == 1]
genpart_charge_ones = genpart_charge[onesnzeroes == 1]
genpart_status_ones = genpart_status[onesnzeroes == 1]
genpart_pid_ones = genpart_pid[onesnzeroes == 1]
    
#np_ones array defining

weight_ones_arr = array('f',10000*[0.])

jet_pt_ones_arr = array('f',10000*[0.])
jet_eta_ones_arr = array('f',10000*[0.])
jet_phi_ones_arr = array('f',10000*[0.])
jet_mass_ones_arr = array('f',10000*[0.])
jet_btag_ones_arr = array('f',10000*[0.])

genjet_pt_ones_arr = array('f',10000*[0.])
genjet_eta_ones_arr = array('f',10000*[0.])
genjet_phi_ones_arr = array('f',10000*[0.])
genjet_mass_ones_arr = array('f',10000*[0.])
genpart_pt_ones_arr = array('f',10000*[0.])
genpart_eta_ones_arr = array('f',10000*[0.])
genpart_phi_ones_arr = array('f',10000*[0.])
genpart_mass_ones_arr = array('f',10000*[0.])
genpart_charge_ones_arr = array('f',10000*[0.])
genpart_status_ones_arr = array('f',10000*[0.])
genpart_pid_ones_arr = array('f',10000*[0.])

#np_ones array branches

tree.Branch("weight_ones", weight_ones_arr, "weight_ones/F")

tree.Branch("jet_pt_ones", jet_pt_ones_arr, "jet_pt_ones/F")
tree.Branch("jet_eta_ones", jet_eta_ones_arr, "jet_eta_ones/F")
tree.Branch("jet_phi_ones", jet_phi_ones_arr, "jet_phi_ones/F")
tree.Branch("jet_mass_ones", jet_mass_ones_arr, "jet_mass_ones/F")
tree.Branch("jet_btag_ones", jet_btag_ones_arr, "jet_btag_ones/F")

tree.Branch("genjet_pt_ones", genjet_pt_ones_arr, "genjet_pt_ones/F")
tree.Branch("genjet_eta_ones", genjet_eta_ones_arr, "genjet_eta_ones/F")
tree.Branch("genjet_phi_ones", genjet_phi_ones_arr, "genjet_phi_ones/F")
tree.Branch("genjet_mass_ones", genjet_mass_ones_arr, "genjet_mass_ones/F")

tree.Branch("genpart_pt_ones", genpart_pt_ones_arr, "genpart_pt_ones/F")
tree.Branch("genpart_eta_ones", genpart_eta_ones_arr, "genpart_eta_ones/F")
tree.Branch("genpart_phi_ones", genpart_phi_ones_arr, "genpart_phi_ones/F")
tree.Branch("genpart_mass_ones", genpart_mass_ones_arr, "genpart_mass_ones/F")
tree.Branch("genpart_charge_ones", genpart_charge_ones_arr, "genpart_charge_ones/F")
tree.Branch("genpart_status_ones", genpart_status_ones_arr, "genpart_status_ones/F")
tree.Branch("genpart_pid_ones", genpart_pid_ones_arr, "genpart_pid_ones/F")

# Lbar creation
ljet_eta = np.array(ljet_eta)
sljet_eta = np.array(sljet_eta)

ljet_phi = np.array(ljet_phi)
sljet_phi = np.array(sljet_phi)

 
l_eta = np.array(l_eta)
sl_eta = np.array(sl_eta)

l_phi = np.array(l_phi)
sl_phi = np.array(sl_phi)

bbbar_dphi = abs(abs(abs(ljet_phi-sljet_phi)-np.pi)-np.pi)
bbbar_deta = abs(ljet_eta-sljet_eta)

llbar_dphi = abs(abs(abs(l_phi-sl_phi)-np.pi)-np.pi)
llbar_deta = abs(l_eta-sl_eta)

#New arrays & branches for bar
bbbar_dphi_arr = array('f', [0.])
bbbar_deta_arr = array('f', [0.])
llbar_dphi_arr = array('f', [0.])
llbar_deta_arr = array('f', [0.])

tree.Branch("bbbar_dphi", bbbar_dphi_arr, "bbbar_dphi/F")
tree.Branch("bbbar_deta", bbbar_deta_arr, "bbbar_deta/F")
tree.Branch("llbar_dphi", llbar_dphi_arr, "llbar_dphi/F")
tree.Branch("llbar_deta", llbar_deta_arr, "llbar_deta/F")


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
    
    bbbar_dphi_arr[0] = bbbar_dphi[i]
    bbbar_deta_arr[0] = bbbar_deta[i]

    llbar_dphi_arr[0] = llbar_dphi[i]
    llbar_deta_arr[0] = llbar_deta[i]

    for v in range(len(weight[i])):
        weight_ones_arr[0][v] = weight[i][v]

    for v in range(len(genjet_pt[i])):
        genjet_pt_ones_arr[0][v] = genjet_pt[v][i]
        genjet_eta_ones_arr[0][v] = genjet_eta[v][i]
        genjet_phi_ones_arr[0][v] = genjet_phi[v][i]
        genjet_mass_ones_arr[0][v] = genjet_mass[v][i]
    for v in range(len(genpart_pt[i])):
        genpart_pt_ones_arr[0][v] = genpart_pt[v][i]
        genpart_eta_ones_arr[0][v] = genpart_eta[v][i]
        genpart_phi_ones_arr[0][v] = genpart_phi[v][i]
        genpart_mass_ones_arr[0][v] = genpart_mass[v][i]
        genpart_charge_ones_arr[0][v] = genpart_charge[v][i]
        genpart_pid_ones_arr[0][v] = genpart_pid[v][i]
        genpart_status_ones_arr[0][v] = genpart_status[v][i]
    for v in range(len(jet_pt[i])):
        jet_pt_ones_arr[0][v] = jet_pt[v][i]
        jet_eta_ones_arr[0][v] = jet_eta[v][i]
        jet_phi_ones_arr[0][v] = jet_phi[v][i]
        jet_mass_ones_arr[0][v] = jet_mass[v][i]
        jet_btag_ones_arr[0][v] = jet_btag[v][i]
    
    tree.Fill()

outputfile.Write()
outputfile.Close()

print(len(l_pt))
print(counter)
#print(j_pt)
