import ROOT
import uproot
import numpy as np
from   array import array
#import matplotlib.pypolt as plt
import argparse

gen_lep_pt = []
gen_lep_eta = []
gen_lep_phi = []
gen_lep_mass = []

gen_alep_pt = []
gen_alep_eta = []
gen_alep_phi = []
gen_alep_mass = []

gen_top_pt = []
gen_top_eta = []
gen_top_phi = []
gen_top_mass = []

gen_atop_pt = []
gen_atop_eta = []
gen_atop_phi = []
gen_atop_mass = []

num_leps = 0
num_antileps = 0
lep_idx = []
antilep_idx = []

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input')
args = parser.parse_args()

fileptr = uproot.open(args.input)

genpart_pt = fileptr['CutTree']['genpart_pt_ones'].array()
genpart_eta = fileptr['CutTree']['genpart_eta_ones'].array()
genpart_phi = fileptr['CutTree']['genpart_phi_ones'].array()
genpart_mass = fileptr['CutTree']['genpart_mass_ones'].array()
genpart_pid = fileptr['CutTree']['genpart_pid_ones'].array()
genpart_status = fileptr['CutTree']['genpart_status_ones'].array()
genpart_charge = fileptr['CutTree']['genpart_charge_ones'].array()




for i in range(len(genpart_pt)):
    num_leps = 0
    num_antileps = 0
    lep_idx = 0
    antilep_idx = 0

    for j in range(len(genpart_pt(i-1))):
        if (genpart_pid[i][j] == 11 or 13 or 15 ) and (genpart_pid[i][j+1] == -12 or -14 or -16) and num_leps == 0:
            lep_idx = j
            num_leps += 1

        if (genpart_pid[i][j] == -11 or -13 or -15) and (genpart_pid[i][j + 1] == 12 or 14 or 16) and num_leps == 0:
            antilep_idx = j
            num_antileps += 1

    if num_lep == 0 or num_antileps == 0:
        continue

    gen_lep_pt.append(genpart_pt[i][lep_idx])
    gen_lep_eta.append(genpart_pt[i][lep_idx])
    gen_lep_phi.append(genpart_pt[i][lep_idx])
    gen_lep_mass.append(genpart_pt[i][lep_idx])

    gen_alep_pt.append(genpart_pt[i][lep_idx])
    gen_alep_eta.append(genpart_pt[i][lep_idx])
    gen_alep_phi.append(genpart_pt[i][lep_idx])
    gen_alep_mass.append(genpart_pt[i][lep_idx])

    gen_top_pt.append(genpart_pt[i][2])
    gen_top_eta.append(genpart_pt[i][2])
    gen_top_phi.append(genpart_pt[i][2])
    gen_top_mass.append(genpart_pt[i][2])

    gen_atop_pt.append(genpart_pt[i][3])
    gen_atop_eta.append(genpart_pt[i][3])
    gen_top_phi.append(genpart_pt[i][3])
    gen_top_mass.append(genpart_pt[i][3])

outputfile = ROOT.TFile(args.imput, "recreate")    

tree = ROOT.TTree("GenSliceTree", "GenSliceTree")

# arrays and branches
gen_lep_pt_arr = array('f', [0.])
gen_lep_eta_arr = array('f', [0.])
gen_lep_phi_arr = array('f', [0.])
gen_lep_mass_arr = array('f', [0.])

gen_alep_pt_arr = array('f', [0.])
gen_alep_eta_arr = array('f', [0.])
gen_alep_phi_arr = array('f', [0.])
gen_alep_mass_arr = array('f', [0.])

gen_top_pt_arr = array('f', [0.])
gen_top_eta_arr = array('f', [0.])
gen_top_phi_arr = array('f', [0.])
gen_top_mass_arr = array('f', [0.])

gen_atop_pt_arr = array('f', [0.])
gen_atop_eta_arr = array('f', [0.])
gen_atop_phi_arr = array('f', [0.])
gen_atop_mass_arr = array('f', [0.])

tree.Branch("gen_lep_pt", gen_lep_pt_arr, "gen_lep_pt/F")
tree.Branch("gen_lep_eta", gen_lep_eta_arr, "gen_lep_eta/F")
tree.Branch("gen_lep_phi", gen_lep_phi_arr, "gen_lep_phi/F")
tree.Branch("gen_lep_mass", gen_lep_mass_arr, "gen_lep_mass/F")

tree.Branch("gen_alep_pt", gen_alep_pt_arr, "gen_alep_pt/F")
tree.Branch("gen_alep_eta", gen_alep_eta_arr, "gen_alep_eta/F")
tree.Branch("gen_alep_phi", gen_alep_phi_arr, "gen_alep_phi/F")
tree.Branch("gen_alep_mass", gen_alep_mass_arr, "gen_alep_mass/F")

tree.Branch("gen_top_pt", gen_top_pt_arr, "gen_top_pt/F")
tree.Branch("gen_top_eta", gen_top_eta_arr, "gen_top_eta/F")
tree.Branch("gen_top_phi", gen_top_phi_arr, "gen_top_phi/F")
tree.Branch("gen_top_mass", gen_top_mass_arr, "gen_top_mass/F")

tree.Branch("gen_atop_pt", gen_atop_pt_arr, "gen_atop_pt/F")
tree.Branch("gen_atop_eta", gen_atop_eta_arr, "gen_atop_eta/F")
tree.Branch("gen_atop_phi", gen_atop_phi_arr, "gen_atop_phi/F")
tree.Branch("gen_atop_mass", gen_atop_mass_arr, "gen_atop_mass/F")

for i in range(len(gen_lep_pt)):
    gen_lep_pt_arr[0] = gen_lep_pt[i]
    gen_lep_eta_arr[0] = gen_lep_eta[i]
    gen_lep_phi_arr[0] = gen_lep_phi[i]
    gen_lep_mass_arr[0] = gen_lep_mass[i]
    
    gen_alep_pt_arr[0] = gen_alep_pt[i]
    gen_alep_eta_arr[0] = gen_alep_eta[i]
    gen_alep_phi_arr[0] = gen_alep_phi[i]
    gen_alep_mass_arr[0] = gen_alep_mass[i]
    
    gen_top_pt_arr[0] = gen_top_pt[i]
    gen_top_eta_arr[0] = gen_top_eta[i]
    gen_top_phi_arr[0] = gen_top_phi[i]
    gen_top_mass_arr[0] = gen_top_mass[i]
    
    gen_atop_pt_arr[0] = gen_atop_pt[i]
    gen_atop_eta_arr[0] = gen_atop_eta[i]
    gen_atop_phi_arr[0] = gen_atop_phi[i]
    gen_atop_mass_arr[0] = gen_atop_mass[i]
    
    tree.Fill()
    
outputfile.Write()
outputfile.Close()
    

    
