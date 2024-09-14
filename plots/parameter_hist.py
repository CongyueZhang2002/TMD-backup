#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:41:01 2019

@author: john
"""

import pandas as pd
import numpy as np
import pylab as py
import warnings
from  matplotlib import rc
import matplotlib
#matplotlib.rcParams.update({'errorbar.cB2size': 2})
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
warnings.filterwarnings('ignore')

xlabels = [r'\rm $a_N$',r'\rm $b_N$', r'\rm $\alpha$',  r'\rm $\beta$',  r'\rm $g_{2B}$',  r'\rm $b_2$']
dfinit = pd.read_csv('./collected_results.txt', delimiter = ',',delim_whitespace=True)

print ('aN')
print (np.sort(dfinit['aN'].tolist())[16] -dfinit['aN'].tolist()[0])
print (dfinit['aN'].tolist()[0])
print (np.sort(dfinit['aN'].tolist())[84]-dfinit['aN'].tolist()[0])

print ('bN')
print (np.sort(dfinit['bN'].tolist())[16] -dfinit['bN'].tolist()[0])
print (dfinit['bN'].tolist()[0])
print (np.sort(dfinit['bN'].tolist())[84]-dfinit['bN'].tolist()[0])


print ('g2A')
print (np.sort(dfinit['g2A'].tolist())[16] -dfinit['g2A'].tolist()[0])
print (dfinit['g2A'].tolist()[0])
print (np.sort(dfinit['g2A'].tolist())[84]-dfinit['g2A'].tolist()[0])


spacing = {}
spacing['aN']    = 0.01
spacing['bN']    = 0.01
spacing['g2A']    = 0.01
spacing['a2']    = 0.01

fig = py.figure(figsize = (20,20))
i = 0
keys = ['aN','bN']
for i in range(len(keys)):
  key = keys[i]
  print (i, key)
  data = dfinit[key]
  print (dfinit[key])
  print("standard deviation for", key, ":",  np.std(data) )
  #binss = np.linspace(min(dfinit[key]),max(dfinit[key]),int(round((max(dfinit[key])-min(dfinit[key]))/spacing[key],0)))
  binss = np.linspace(min(data),max(data),20)
  ax = py.subplot(4,4,i+1)
  ax.hist(data, bins = binss)
  ax.axvline(dfinit[key].iloc[0],color = 'k')
  xlabel = xlabels[i]
  ax.set_xlabel(str(xlabel), fontsize = 40)
  if i == 3:
    ax.set_ylabel(r"$\rm frequency$", fontsize = 40,y=-0.)
  py.tick_params(labelsize=30)
  i+= 1
  axs = py.gca()
py.tight_layout(pad = 0.)
fig.savefig('parameter_hist.pdf',bbox_inches='tight')
