import os
import pandas as pd
import numpy as np
import pylab as py
import warnings
import matplotlib.pyplot as plt
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline as spline
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from matplotlib.ticker import ScalarFormatter,MaxNLocator,LogLocator,NullFormatter,FuncFormatter
warnings.filterwarnings('ignore')
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
from matplotlib.pyplot import gca

results = pd.read_csv('params.dat', delim_whitespace=True)

for i in range(len(results)-1):
    results['gamma'][i+1] = results['gamma'][0]
    results['g3f'][i+1] = results['g3f'][0]
    results['g3D'][i+1] = results['g3D'][0]
    results['Nq1'][i+1] = results['Nq1'][0]
    results['gq1'][i+1] = results['gq1'][0]
    results['dq1'][i+1] = results['dq1'][0]
    results['Nq2'][i+1] = results['Nq2'][0]
    results['gq2'][i+1] = results['gq2'][0]
    results['dq2'][i+1] = results['dq2'][0]
    results['p_10'][i+1] = results['p_10'][0]
    results['p_11'][i+1] = results['p_11'][0]
    results['p_12'][i+1] = results['p_12'][0]
    results['rep'][i+1] = i+1

df =[]
for i in range(0,176): 
    df.append(results[i:i+1])

for j in range(0,176): 
    if len(str(j)) == 1:
        numb = '00'+str(j)
    elif len(str(j)) == 2:
        numb = '0'+str(j)
    else:
        numb = str(j)
    dest = './replica_'+numb + '/expdata/'
    df[j].to_csv(dest+"/collected_results.csv", sep='\t', index = False)
