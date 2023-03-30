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
