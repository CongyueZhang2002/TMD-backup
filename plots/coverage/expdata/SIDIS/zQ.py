import pandas as pd
import numpy as np
import os
import pylab as py
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib
#matplotlib.rcParams.update({'errorbar.capsize': 2})
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
#from scipy.interpolate import spline,interp1d
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

dftot = pd.DataFrame(columns=["Nu","Z","Q2","pt2","MULT-RATIO","STAT","SYS","fuu","fuua","DIS"])

for file in os.listdir("./HERMES_DIS/"):
    if "z.dat" in file:
        df = pd.read_csv("HERMES_DIS/"+file, delim_whitespace=True)
        df = df[df.Z<0.7]
        df = df[df.pt2<0.3]
        dftot = pd.concat([dftot,df])

ax = py.subplot(111)
h = ax.hist2d(dftot['Z'],dftot['Q2'],bins = 15,cmap=py.cm.Reds)
cbar = py.colorbar(h[3])

cbar.ax.tick_params(labelsize=30)
ax.tick_params(labelsize = 30)
ax.set_ylabel(r'\rm $Q^2$',fontsize = 30)
ax.set_xlabel(r'\rm $z$'  ,fontsize = 30)
ax.legend(frameon = False,fontsize = 30,loc='upper center', bbox_to_anchor=(1.4, 0.8))
#ax.set_ylim(0.8,3e4)
#ax.semilogy()

py.tight_layout()

py.savefig('filter_SIDIS.pdf', bbox_inches='tight')
