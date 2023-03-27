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