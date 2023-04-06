import uproot
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt

plt.style.use([hep.style.ROOT, hep.style.firamath])

from groupbcuts.py import mu_pt

fileptr = uproot.open("TT_Dilept_13.root")
plt.hist(mu_pt)

