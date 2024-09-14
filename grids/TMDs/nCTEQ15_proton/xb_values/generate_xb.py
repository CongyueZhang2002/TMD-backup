#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#from scipy.interpolate import spline,interp1d
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from matplotlib.ticker import ScalarFormatter,MaxNLocator,LogLocator,NullFormatter,FuncFormatter
warnings.filterwarnings('ignore')
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
from matplotlib.pyplot import gca


# In[15]:


# get list of data files
pwd = os.getcwd()
path2dat = pwd + '/HERMES_DIS/'
datalist = os.listdir(path2dat)


# In[42]:


he_datalist = []
ne_datalist = []
kr_datalist = []
xe_datalist = []

for filename in datalist:
    if ('pi' in filename):
        if('he' in filename):
            he_datalist.append(filename)
        if('ne' in filename):
            ne_datalist.append(filename)
        if('kr' in filename):
            kr_datalist.append(filename)
        if('xe' in filename):
            xe_datalist.append(filename)


# In[43]:


# Mass of the proton (in GeV) 0.938272046d0 
M = 0.938272046

def find_xb(nu, Qs):
    xb = Qs/(2*M*nu)
    return xb


# In[61]:


he_xblist = []
ne_xblist = []
kr_xblist = []
xe_xblist = []

for filename in he_datalist:
    df = pd.read_csv(path2dat+filename,delim_whitespace=True)
    nulist = df['Nu']
    Qslist = df['Q2']  
    for i in range(len(nulist)):
        nu = nulist[i]
        Qs = Qslist[i]
        xb = find_xb(nu,Qs)
        he_xblist.append(xb)
        
for filename in ne_datalist:
    df = pd.read_csv(path2dat+filename,delim_whitespace=True)
    nulist = df['Nu']
    Qslist = df['Q2']  
    for i in range(len(nulist)):
        nu = nulist[i]
        Qs = Qslist[i]
        xb = find_xb(nu,Qs)
        ne_xblist.append(xb)
        
for filename in kr_datalist:
    df = pd.read_csv(path2dat+filename,delim_whitespace=True)
    nulist = df['Nu']
    Qslist = df['Q2']  
    for i in range(len(nulist)):
        nu = nulist[i]
        Qs = Qslist[i]
        xb = find_xb(nu,Qs)
        kr_xblist.append(xb)

for filename in xe_datalist:
    df = pd.read_csv(path2dat+filename,delim_whitespace=True)
    nulist = df['Nu']
    Qslist = df['Q2']  
    for i in range(len(nulist)):
        nu = nulist[i]
        Qs = Qslist[i]
        xb = find_xb(nu,Qs)
        xe_xblist.append(xb)


# In[71]:


#sort from smallest to biggest xb
he_xblist.sort()
ne_xblist.sort()
kr_xblist.sort()
xe_xblist.sort()


# In[84]:


# rm duplicates
he_xblist_rm = list(set(he_xblist))
ne_xblist_rm = list(set(ne_xblist))
kr_xblist_rm = list(set(kr_xblist))
xe_xblist_rm = list(set(xe_xblist))
he_xblist_rm.sort()
ne_xblist_rm.sort()
kr_xblist_rm.sort()
xe_xblist_rm.sort()


# In[85]:


#output as csv files. 
he_df = pd.DataFrame(he_xblist_rm, columns = ['xb'])
ne_df = pd.DataFrame(ne_xblist_rm, columns = ['xb'])
kr_df = pd.DataFrame(kr_xblist_rm, columns = ['xb'])
xe_df = pd.DataFrame(xe_xblist_rm, columns = ['xb'])
he_df.to_csv('he_xb.csv', index=False)
ne_df.to_csv('ne_xb.csv', index=False)
kr_df.to_csv('kr_xb.csv', index=False)
xe_df.to_csv('xe_xb.csv', index=False)


# In[ ]:




