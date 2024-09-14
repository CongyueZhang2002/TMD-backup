import os
import pandas as pd
import numpy as np
import pylab as py
import warnings
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
#from scipy.interpolate import spline,interp1d
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
warnings.filterwarnings('ignore')

src = './fit'

for filename in os.listdir(src+'/expdata/SIDIS/HERMES_DIS/'):
    path2src = src+'/expdata/SIDIS/HERMES_DIS/'
    df = pd.read_csv(path2src+filename, delim_whitespace=True)
    if 'MULT-RATIO' not in df.columns:
        print("{} does not contain 'MULT-RATIO'".format(filename))

