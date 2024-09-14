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


df3d = pd.read_csv('fort.2000', delimiter = ',',delim_whitespace=True)


# In[3]:


from mpl_toolkits import mplot3d
from numpy import nan


# In[4]:


kt = df3d['kt']
Q = df3d['Q']
PDF_P = df3d['up']
PDF_1 = df3d['u1']
PDF_2 = df3d['u2']


# In[35]:


fig = plt.figure()
ax = plt.axes(projection ='3d')

ax.scatter(kt,Q,PDF_P,color= "blue")
ax.set_xlabel(r"\rm $p_T$ (GeV)", fontsize = 25, labelpad = 20)
ax.set_ylabel(r"\rm $Q$ (GeV)", fontsize = 25, labelpad = 20)
ax.zaxis.set_rotate_label(False)  # disable automatic rotation
ax.set_zlabel(r"\rm $D_{\pi^0/p}$", fontsize = 25, labelpad = 20)

ax.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig.set_size_inches(25,14)


# In[36]:


fig.savefig('TMDFF3D_nucleon.pdf', bbox_inches = "tight")


# In[37]:


fig2 = plt.figure()
ax2 = plt.axes(projection ='3d')

ax2.scatter(kt,Q,PDF_1,color= "green")
ax2.set_xlabel(r"\rm $p_T$ (GeV)", fontsize = 25, labelpad = 20)
ax2.set_ylabel(r"\rm $Q$ (GeV)", fontsize = 25, labelpad = 20)
ax2.zaxis.set_rotate_label(False)  # disable automatic rotation
ax2.set_zlabel(r"\rm $D_{\pi^{0}/A}$", fontsize = 25, labelpad = 25)

ax2.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax2.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax2.set_zticks([0.08,0.10,0.12,0.14,0.16])
ax2.set_zlim(0.06,0.17)
ax2.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax2.invert_yax2is()
fig2.set_size_inches(24,14)


# In[38]:


fig2.savefig('TMDFF3D_HE.pdf', bbox_inches = "tight")


# In[39]:


fig3 = plt.figure()
ax3 = plt.axes(projection ='3d')

ax3.scatter(kt,Q,PDF_2,color= "red")
ax3.set_xlabel(r"\rm $p_T$ (GeV)", fontsize = 25, labelpad = 20)
ax3.set_ylabel(r"\rm $Q$ (GeV)", fontsize = 25, labelpad = 20)
ax3.zaxis.set_rotate_label(False)  # disable automatic rotation
ax3.set_zlabel(r"\rm $D_{\pi^0/A}$", fontsize = 25, labelpad = 25)

ax3.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax3.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax3.set_zticks([0.06,0.08,0.10,0.12,0.14,0.16,0.18])
ax3.set_zlim(0.06,0.17)
ax3.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax3.invert_yax3is()
fig3.set_size_inches(24,14)


# In[40]:


fig3.savefig('TMDFF3D_XE.pdf', bbox_inches = "tight")


# In[50]:


fig4 = plt.figure()
ax4 = plt.axes(projection ='3d')

ax4.scatter(kt,Q,PDF_1/PDF_P,color= "green")
ax4.set_xlabel(r"\rm $p_T$ (GeV)", fontsize = 25, labelpad = 20)
ax4.set_ylabel(r"\rm $Q$ (GeV)", fontsize = 25, labelpad = 20)
ax4.zaxis.set_rotate_label(False)  # disable automatic rotation
ax4.set_zlabel(r"\rm $\frac{f_A}{f_{p/n}}$", fontsize = 40, labelpad = 20)
ax4.set_zticks([0.6,0.8,1.0,1.2,1.4])
ax4.set_zlim(0.6,1.2)
ax4.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax4.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax4.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig4.set_size_inches(24,14)


# In[42]:


fig4.savefig('TMDFF3D_HeRatio.pdf', bbox_inches = "tight")


# In[48]:


fig5 = plt.figure()
ax5 = plt.axes(projection ='3d')

ax5.scatter(kt,Q,PDF_2/PDF_P,color= "red")
ax5.set_xlabel(r"\rm $p_T$ (GeV)", fontsize = 25, labelpad = 20)
ax5.set_ylabel(r"\rm $Q$ (GeV)", fontsize = 25, labelpad = 20)
ax5.zaxis.set_rotate_label(False)  # disable automatic rotation
ax5.set_zlabel(r"\rm $\frac{f_A}{f_{p/n}}$", fontsize = 40, labelpad = 20)
ax5.set_zticks([0.6,0.8,1.0,1.2,1.4])
ax5.set_zlim(0.6,1.2)
ax5.tick_params(axis = 'x', direction ='in',labelsize = 25)
ax5.tick_params(axis = 'y', direction ='in',labelsize = 25)
ax5.tick_params(axis = 'z', direction ='in',labelsize = 25)

#ax.invert_yaxis()
fig5.set_size_inches(24,14)


# In[49]:


fig5.savefig('TMDFF3D_XeRatio.pdf', bbox_inches = "tight")


# In[ ]:




