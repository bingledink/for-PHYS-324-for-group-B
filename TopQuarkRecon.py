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
#tree = ROOT.TTree("CutTree", "CutTree")

genpart_pt = fileptr['CutTree']['genpart_pt'].array()
genpart_eta = fileptr['CutTree']['genpart_eta'].array()
genpart_phi = fileptr['CutTree']['genpart_phi'].array()
genpart_mass = fileptr['CutTree']['genpart_mass'].array()
genpart_pid = fileptr['CutTree']['genpart_pid'].array()
genpart_status = fileptr['CutTree']['genpart_status'].array()
genpart_charge = fileptr['CutTree']['genpart_charge'].array()

met_pt = fileptr['CutTree']['met_pt'].array()
met_phi = fileptr['CutTree']['met_phi'].array()

elec_pt = fileptr['CutTree']['elec_pt'].array()
elec_eta = fileptr['CutTree']['elec_eta'].array()
elec_charge = fileptr['CutTree']['elec_charge'].array()
elec_phi = fileptr['CutTree']['elec_phi'].array()
elec_mass = fileptr['CutTree']['elec_mass'].array()

muon_pt = fileptr['CutTree']['muon_pt'].array()
muon_eta = fileptr['CutTree']['muon_eta'].array()
muon_charge = fileptr['CutTree']['muon_charge'].array()
muon_phi = fileptr['CutTree']['muon_phi'].array()
muon_mass = fileptr['CutTree']['muon_mass'].array()

jet_pt = fileptr['CutTree']['jet_pt'].array()
jet_eta = fileptr['CutTree']['jet_eta'].array()
jet_btag = fileptr['CutTree']['jet_btag'].array()
jet_phi = fileptr['CutTree']['jet_phi'].array()
jet_mass = fileptr['CutTree']['jet_phi'].array()

weight = fileptr['CutTree']['weight'].array()

onesnzeroes = [0] * len(elec_pt)

e4vector = ROOT.TLorentzVector()
mu4vector = ROOT.TLorentzVector()

lep_pt = []
lep_eta = []
lep_phi = []
lep_mass = []

alep_pt = []
alep_eta = []
alep_phi = []
alep_mass = []

top_pt = []
top_eta = []
top_phi = []
top_mass = []

atop_pt = []
atop_eta = []
atop_phi = []
atop_mass = []

b_pt = []
b_eta = []
b_phi = []
b_mass = []

bbar_pt = []
bbar_eta = []
bbar_phi = []
bbar_mass = []

nu_pt = []
nu_eta = []
nu_phi = []

anu_pt = []
anu_eta = []
anu_phi = []

tt_mass = []
gentt_mass = []

gen_top_pt = []
gen_top_eta = []
gen_top_phi = []
gen_top_rap = []

gen_atop_pt = []
gen_atop_eta = []
gen_atop_phi = []
gen_atop_rap = []

for event_idx in range(len(elec_pt)):
    lep4vector = ROOT.TLorentzVector
    alep4vector = ROOT.TLorentzVector
    if elec_charge[event_idx] == -1 and muon_charge[event_idx] == 1:
        lep4vector.SetPtEtaPhiM((elec_pt[event_idx]), (elec_eta[event_idx]),
                              (elec_phi[event_idx]), (elec_mass[event_idx]))
        alep4vector.SetPtEtaPhiM((muon_pt[event_idx]), (muon_eta[event_idx]),
                               (muon_phi[event_idx]), (muon_mass[event_idx]))

    elif elec_charge[event_idx] == 1 and muon_charge[event_idx] == -1:
        alep4vector.SetPtEtaPhiM((elec_pt[event_idx]), (elec_eta[event_idx]),
                              (elec_phi[event_idx]), (elec_mass[event_idx]))
        lep4vector.SetPtEtaPhiM((muon_pt[event_idx]), (muon_eta[event_idx]),
                               (muon_phi[event_idx]), (muon_mass[event_idx]))
    met_x = met_pt[event_idx] * cos(met_phi[event_idx])
    met_y = met_pt[event_idx] * sin(met_phi[event_idx])
    
    btag_counter = 0
    high_w = 0
    tt_mass_final = 0
    
    for i_1 in range(len(jet_pt[event_idx])):
        
        for i_2 in range(len(jet_pt[event_idx])):
            if i_1 >= i_2:
                continue
            if jet_pt[event_idx][i_1] < 30:
                continue
            if jet_pt[event_idx][i_2] < 30:
                continue
            if abs(jet_eta[event_idx][i_1]) > 2.4:
                continue
            if abs(jet_eta[event_idx][i_2]) > 2.4:
                continue
            if jet_btag[event_idx][i_1] ==0 and jet_btag[event_idx][i_2] == 0:
                continue
            if(getptag[event_idx][i_1] != 0 and getptag[event_idx][i_2]):
                tt_mass_1, top_p4_1, atop_p4_1, new_p4_1, newbar_p4_1, sw_1 = try_smear(jet_1, jet_2, alep4vector, lep4vector, met_x, met_y, event_idx)
                tt_mass_2, top_p4_2, atop_p4_2, new_p4_2, newbar_p4_2, sw_2 = try_smear(jet_2, jet_1, alep4vector, lep4vector, met_x, met_y, event_idx)
    
                if(tt_mass_1 == -999 and tt_mass_2 == -999):
                    continue

                btag_counter = 2
    
                if(tt_mass_1 == -999):
                    tt_mass_final = tt_mass_2
                    top_p4_final = top_p4_2
                    atop_p4_final = atop_p4_2
                    new_p4_final = new_p4_2
                    anew_p4_final = anew_p4_2
                    b_v4_final = jet_2
                    bbar_v4_final = jet_1
        
                if(tt_mass_2 == -999):
                    tt_mass_final = tt_mass_1
                    top_p4_final = top_p4_1
                    atop_p4_final = atop_p4_1
                    new_p4_final = new_p4_1
                    anew_p4_final = anew_p4_1
                    b_v4_final = jet_1
                    bbar_v4_final = jet_2
    
                if(tt_mass_1 != -999 and tt_mass_2 != -999 and sw_2 <= sw_1):
                    tt_mass_final = tt_mass_1
                    top_p4_final = top_p4_1
                    atop_p4_final = atop_p4_1
                    new_p4_final = new_p4_1
                    anew_p4_final = anew_p4_1
                    b_v4_final = jet_1
                    bbar_v4_final = jet_2
        
                if(tt_mass_1 != -999 and tt_mass_2 != -999 and sw_1 <= sw_2):
                    tt_mass_final = tt_mass_2
                    top_p4_final = top_p4_2
                    atop_p4_final = atop_p4_2
                    new_p4_final = new_p4_2
                    anew_p4_final = anew_p4_2
                    b_v4_final = jet_2
                    bbar_v4_final = jet_1
                continue
            if(jet_btag[event_idx][i_1] + jet_btag[event_idx][i_2] == 1):
                tt_mass_1, top_p4_1, atop_p4_1, new_p4_1, newbar_p4_1, sw_1 = try_smear(jet_1, jet_2, alep4vector, lep4vector, met_x, met_y, event_idx)
                tt_mass_2, top_p4_2, atop_p4_2, new_p4_2, newbar_p4_2, sw_2 = try_smear(jet_2, jet_1, alep4vector, lep4vector, met_x, met_y, event_idx)
                if(tt_mass_1 == -999 and tt_mass_2 == -999):
                    continue
                if(tt_mass_2 == -999 and high_w <= sw_1):
                    tt_mass_final = tt_mass_1
                    top_p4_final = top_p4_1
                    atop_p4_final = atop_p4_1
                    new_p4_final = new_p4_1
                    anew_p4_final = anew_p4_1
                    b_v4_final = jet_1
                    bbar_v4_final = jet_2
                if(tt_mass_1 == -999 and high_w <= sw_2):
                    tt_mass_final = tt_mass_2
                    top_p4_final = top_p4_2
                    atop_p4_final = atop_p4_2
                    new_p4_final = new_p4_2
                    anew_p4_final = anew_p4_2
                    b_v4_final = jet_2
                    bbar_v4_final = jet_1
                if(tt_mass_1 != -999 and tt_mass_2 != -999 and sw_2 <= sw_1 and high_w < sw_1):
                    tt_mass_final = tt_mass_1
                    top_p4_final = top_p4_1
                    atop_p4_final = atop_p4_1
                    new_p4_final = new_p4_1
                    anew_p4_final = anew_p4_1
                    b_v4_final = jet_1
                    bbar_v4_final = jet_2
                if(tt_mass_1 != -999 and tt_mass_2 != -999 and sw_1 <= sw_2 and high_w < sw_2):
                    tt_mass_final = tt_mass_2
                    top_p4_final = top_p4_2
                    atop_p4_final = atop_p4_2
                    new_p4_final = new_p4_2
                    anew_p4_final = anew_p4_2
                    b_v4_final = jet_2
                    bbar_v4_final = jet_1
                continue
    if(tt_mass_final == 0):
        continue
    for i in range(len(genpart_pt[event_idx])):
        if(genpart_pid[event_idx][i] == 6 and genpart_status[event_idx][i] == 62):
            gen_top_4vec == ROOT.TLorentzVector()
            gen_top_4vec.SetPtEtaPhiM(genpart_pt[event_idx][i], genpart_eta[event_idx][i], genpart_phi[event_idx][i], genpart_mass[event_idx][i])
        if(genpart_pid[event_idx][i] == -6 and genpart_status[event_idx][i] == 62):
            gen_atop_4vec == ROOT.TLorentzVector()
            gen_atop_4vec.SetPtEtaPhiM(genpart_pt[event_idx][i], genpart_eta[event_idx][i], genpart_phi[event_idx][i], genpart_mass[event_idx][i])
    com_4vec = gen_top_4vec + gen_atop_4vec
    tt_mass.append(tt_mass_final)
    top_pt.append(top_p4_final.Pt())

            
            
    
        
    
    


