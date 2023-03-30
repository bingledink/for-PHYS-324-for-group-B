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

for i in range (len(elec_pt)):
    lep4vector = ROOT.TLorentzVector
    alep4vector = ROOT.TLorentzVector
    


