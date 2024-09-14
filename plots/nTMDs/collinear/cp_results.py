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

src = '../../'
os.listdir(src)
results = pd.read_csv(src+'collected_results.txt', delim_whitespace=True)

for i in range(len(results)-1):
    results['aN'][i+1] = results['aN'][0]
    results['bN'][i+1] = results['bN'][0]


df =[]
for i in range(97):
    df.append(results[i:i+1])


for j in range(97):
    if len(str(j)) == 1:
        numb = '00'+str(j)
    elif len(str(j)) == 2:
        numb = '0'+str(j)
    else:
        numb = str(j)
    dest = './replica_'+numb
    df[j].to_csv(dest+"/collected_results.csv", sep='\t', index = False)
