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
parser.add_argument('-o', '--output', help='Output')
args = parser.parse_args()

fileptr = uproot.open(args.input)
outputfile = ROOT.TFile(args.output, 'recreate')
tree = ROOT.TTree("CutTree", "CutTree")

genpart_pt = fileptr['CutTree']['genpart_pt'].array()
genpart_eta = fileptr['CutTree']['genpart_eta'].array()
genpart_phi = fileptr['CutTree']['genpart_phi'].array()
genpart_mass = fileptr['CutTree']['genpart_mass'].array()
genpart_pid = fileptr['CutTree']['genpart_pid'].array()
genpart_status = fileptr['CutTree']['genpart_status'].array()
genpart_charge = fileptr['CutTree']['genpart_charge'].array()




for i in range(len(genpart_pt)):
    num_leps = 0
    num_antileps = 0
    lep_idx = 0
    antilep_idx = 0

    for j in range(len(genpart_pt(i-1))):
        if (genpart_pid[i][j] == 11 or 13 or 15 ) and (genpart_pid[i][j+1] == -12 or -14 or -16) and num_leps = 0
            lep_idx = j
            num_leps += 1

        if (genpart_pid[i][j] == -11 or -13 or -15) and (genpart_pid[i][j + 1] == 12 or 14 or 16) and num_leps = 0
            antilep_idx = j
            num_antileps += 1

    if num_lep = 0 or num_antileps = 0
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

# arrays and branches
genpart_pt_arr = array('f', [0.])
genpart_eta_arr = array('f', [0.])
genpart_phi_arr = array('f', [0.])
genpart_mass_arr = array('f', [0.])
genpart_pid_arr = array('f', [0.])
genpart_status_arr = array('f', [0.])
genpart_charge_arr = array('f', [0.])

tree.Branch("genpart_pt", genpart_pt_arr, "genpart_pt/F")
tree.Branch("genpart_eta", genpart_eta_arr, "genpart_eta/F")
tree.Branch("genpart_phi", genpart_phi_arr, "genpart_phi/F")
tree.Branch("genpart_mass", genpart_mass_arr, "genpart_mass/F")
tree.Branch("genpart_pid", genpart_pid_arr, "genpart_pid/F")
tree.Branch("genpart_status", genpart_status_arr, "genpart_status/F")
tree.Branch("genpart_charge", genpart_charge_arr, "genpart_charge/F")

