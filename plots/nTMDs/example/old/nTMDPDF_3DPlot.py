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
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline as spline
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from matplotlib.ticker import ScalarFormatter,MaxNLocator,LogLocator,NullFormatter,FuncFormatter
warnings.filterwarnings('ignore')
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
from matplotlib.pyplot import gca


# In[2]:


df3d = pd.read_csv('fort.3000', delimiter = ',',delim_whitespace=True)


# In[3]:


from mpl_toolkits import mplot3d
from numpy import nan


# In[4]:


kt = df3d['kt']
x = df3d['x']
PDF_P = x*df3d['up']
PDF_1 = x*df3d['u1']
PDF_2 = x*df3d['u2']

xhigh = 0.4
for i in range(len(x)):
    if (x[i] > xhigh):
        x[i]     = nan
        kt[i]    = nan
        PDF_P[i] = nan
        PDF_1[i] = nan
        PDF_2[i] = nan


# In[5]:


fig = plt.figure()
ax = plt.axes(projection ='3d')

ax.scatter(kt,x,PDF_P,color= "blue")
ax.set_xlabel(r"\rm $k_T$ (GeV)", fontsize = 25, labelpad = 20)
ax.set_ylabel(r"\rm $x$", fontsize = 25, labelpad = 20)
ax.zaxis.set_rotate_label(False)  # disable automatic rotation
ax.set_zticks([0.0,0.1,0.2,0.3,0.4])
ax.set_zlabel(r"\rm $x f_{u/p}$", fontsize = 25, labelpad = 25)

ax.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig.set_size_inches(24,14)


# In[16]:


fig.savefig('TMD3D_nucleon.pdf', bbox_inches = "tight")


# In[9]:


fig2 = plt.figure()
ax2 = plt.axes(projection ='3d')

ax2.scatter(kt,x,PDF_1,color= "green")
ax2.set_xlabel(r"\rm $k_T$ (GeV)", fontsize = 25, labelpad = 20)
ax2.set_ylabel(r"\rm $x$", fontsize = 25, labelpad = 20)
ax2.zaxis.set_rotate_label(False)  # disable automatic rotation
ax2.set_zlabel(r"\rm $x f_{u/A}$", fontsize = 25, labelpad = 25)

ax2.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax2.set_yticks([0.0,0.1,0.2,0.3,0.4])
ax2.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax2.set_zticks([0.0,0.1,0.2,0.3,0.4])
ax2.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax2.invert_yax2is()
fig2.set_size_inches(24,14)


# In[10]:


fig2.savefig('TMD3D_HE.pdf', bbox_inches = "tight")


# In[12]:


fig3 = plt.figure()
ax3 = plt.axes(projection ='3d')

ax3.scatter(kt,x,PDF_2,color= "red")
ax3.set_xlabel(r"\rm $k_T$ (GeV)", fontsize = 25, labelpad = 20)
ax3.set_ylabel(r"\rm $x$", fontsize = 25, labelpad = 20)
ax3.zaxis.set_rotate_label(False)  # disable automatic rotation
ax3.set_zlabel(r"\rm $x f_{u/A}$", fontsize = 25, labelpad = 25)

ax3.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax3.set_yticks([0.0,0.1,0.2,0.3,0.4])
ax3.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax3.set_zticks([0.0,0.1,0.2,0.3,0.4])
#ax3.set_zlim(0,2.2)
ax3.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax3.invert_yax3is()
fig3.set_size_inches(24,14)


# In[13]:


fig3.savefig('TMD3D_XE.pdf', bbox_inches = "tight")


# In[16]:


fig4 = plt.figure()
ax4 = plt.axes(projection ='3d')

ax4.scatter(kt,x,PDF_1/PDF_P,color= "green")
ax4.set_xlabel(r"\rm $k_T$ (GeV)", fontsize = 25, labelpad = 20)
ax4.set_ylabel(r"\rm $x$", fontsize = 25, labelpad = 20)
ax4.zaxis.set_rotate_label(False)  # disable automatic rotation
ax4.set_yticks([0.0,0.1,0.2,0.3,0.4])
ax4.set_zlabel(r"\rm $\frac{f_A}{f_{p/n}}$", fontsize = 40, labelpad = 20)
ax4.set_zticks([0.6,0.8,1.0,1.2,1.4])
ax4.set_zlim(0.6,1.4)
ax4.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax4.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax4.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig4.set_size_inches(24,14)


# In[17]:


fig4.savefig('TMD3D_HeRatio.pdf', bbox_inches = "tight")


# In[24]:


fig5 = plt.figure()
ax5 = plt.axes(projection ='3d')

ax5.scatter(kt,x,PDF_2/PDF_P,color= "red")
ax5.set_xlabel(r"\rm $k_T$ (GeV)", fontsize = 25, labelpad = 20)
ax5.set_ylabel(r"\rm $x$", fontsize = 25, labelpad = 20)
ax5.zaxis.set_rotate_label(False)  # disable automatic rotation
ax5.set_yticks([0.0,0.1,0.2,0.3,0.4])
ax5.set_zlabel(r"\rm $\frac{f_A}{f_{p/n}}$", fontsize = 40, labelpad = 20)
ax5.set_zticks([0.6,0.8,1.0,1.2,1.4])
ax5.set_zlim(0.4,1.3)
ax5.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax5.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax5.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig5.set_size_inches(24,14)


# In[25]:


fig5.savefig('TMD3D_XeRatio.pdf', bbox_inches = "tight")


# In[ ]:
