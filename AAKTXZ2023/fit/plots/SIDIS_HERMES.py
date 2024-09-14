#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex = True)


# In[13]:
pion_data_dir = '../plot_data/HERMES/'
kaon_data_dir = '../plot_data/HERMES/'

#Get Data
#pip
pip_he_nu_df = pd.read_csv(pion_data_dir + 'pip_he_nu.dat', sep ='\s+')
pip_ne_nu_df = pd.read_csv(pion_data_dir + 'pip_ne_nu.dat', sep ='\s+')
pip_kr_nu_df = pd.read_csv(pion_data_dir + 'pip_kr_nu.dat', sep ='\s+')
pip_xe_nu_df = pd.read_csv(pion_data_dir + 'pip_xe_nu.dat', sep ='\s+')

pip_he_z_df = pd.read_csv(pion_data_dir + 'pip_he_z.dat', sep ='\s+')
pip_ne_z_df = pd.read_csv(pion_data_dir + 'pip_ne_z.dat', sep ='\s+')
pip_kr_z_df = pd.read_csv(pion_data_dir + 'pip_kr_z.dat', sep ='\s+')
pip_xe_z_df = pd.read_csv(pion_data_dir + 'pip_xe_z.dat', sep ='\s+')

pip_he_q2_df = pd.read_csv(pion_data_dir + 'pip_he_q2.dat', sep ='\s+')
pip_ne_q2_df = pd.read_csv(pion_data_dir + 'pip_ne_q2.dat', sep ='\s+')
pip_kr_q2_df = pd.read_csv(pion_data_dir + 'pip_kr_q2.dat', sep ='\s+')
pip_xe_q2_df = pd.read_csv(pion_data_dir + 'pip_xe_q2.dat', sep ='\s+')

pip_he_pt2_df = pd.read_csv(pion_data_dir + 'pip_he_pt2.dat', sep ='\s+')
pip_ne_pt2_df = pd.read_csv(pion_data_dir + 'pip_ne_pt2.dat', sep ='\s+')
pip_kr_pt2_df = pd.read_csv(pion_data_dir + 'pip_kr_pt2.dat', sep ='\s+')
pip_xe_pt2_df = pd.read_csv(pion_data_dir + 'pip_xe_pt2.dat', sep ='\s+')

#kp
kp_he_nu_df = pd.read_csv(kaon_data_dir + 'kp_he_nu.dat', sep ='\s+')
kp_ne_nu_df = pd.read_csv(kaon_data_dir + 'kp_ne_nu.dat', sep ='\s+')
kp_kr_nu_df = pd.read_csv(kaon_data_dir + 'kp_kr_nu.dat', sep ='\s+')
kp_xe_nu_df = pd.read_csv(kaon_data_dir + 'kp_xe_nu.dat', sep ='\s+')

kp_he_z_df = pd.read_csv(kaon_data_dir + 'kp_he_z.dat', sep ='\s+')
kp_ne_z_df = pd.read_csv(kaon_data_dir + 'kp_ne_z.dat', sep ='\s+')
kp_kr_z_df = pd.read_csv(kaon_data_dir + 'kp_kr_z.dat', sep ='\s+')
kp_xe_z_df = pd.read_csv(kaon_data_dir + 'kp_xe_z.dat', sep ='\s+')

kp_he_q2_df = pd.read_csv(kaon_data_dir + 'kp_he_q2.dat', sep ='\s+')
kp_ne_q2_df = pd.read_csv(kaon_data_dir + 'kp_ne_q2.dat', sep ='\s+')
kp_kr_q2_df = pd.read_csv(kaon_data_dir + 'kp_kr_q2.dat', sep ='\s+')
kp_xe_q2_df = pd.read_csv(kaon_data_dir + 'kp_xe_q2.dat', sep ='\s+')

kp_he_pt2_df = pd.read_csv(kaon_data_dir + 'kp_he_pt2.dat', sep ='\s+')
kp_ne_pt2_df = pd.read_csv(kaon_data_dir + 'kp_ne_pt2.dat', sep ='\s+')
kp_kr_pt2_df = pd.read_csv(kaon_data_dir + 'kp_kr_pt2.dat', sep ='\s+')
kp_xe_pt2_df = pd.read_csv(kaon_data_dir + 'kp_xe_pt2.dat', sep ='\s+')

#pim
pim_he_nu_df = pd.read_csv(pion_data_dir + 'pim_he_nu.dat', sep ='\s+')
pim_ne_nu_df = pd.read_csv(pion_data_dir + 'pim_ne_nu.dat', sep ='\s+')
pim_kr_nu_df = pd.read_csv(pion_data_dir + 'pim_kr_nu.dat', sep ='\s+')
pim_xe_nu_df = pd.read_csv(pion_data_dir + 'pim_xe_nu.dat', sep ='\s+')

pim_he_z_df = pd.read_csv(pion_data_dir + 'pim_he_z.dat', sep ='\s+')
pim_ne_z_df = pd.read_csv(pion_data_dir + 'pim_ne_z.dat', sep ='\s+')
pim_kr_z_df = pd.read_csv(pion_data_dir + 'pim_kr_z.dat', sep ='\s+')
pim_xe_z_df = pd.read_csv(pion_data_dir + 'pim_xe_z.dat', sep ='\s+')

pim_he_q2_df = pd.read_csv(pion_data_dir + 'pim_he_q2.dat', sep ='\s+')
pim_ne_q2_df = pd.read_csv(pion_data_dir + 'pim_ne_q2.dat', sep ='\s+')
pim_kr_q2_df = pd.read_csv(pion_data_dir + 'pim_kr_q2.dat', sep ='\s+')
pim_xe_q2_df = pd.read_csv(pion_data_dir + 'pim_xe_q2.dat', sep ='\s+')

pim_he_pt2_df = pd.read_csv(pion_data_dir + 'pim_he_pt2.dat', sep ='\s+')
pim_ne_pt2_df = pd.read_csv(pion_data_dir + 'pim_ne_pt2.dat', sep ='\s+')
pim_kr_pt2_df = pd.read_csv(pion_data_dir + 'pim_kr_pt2.dat', sep ='\s+')
pim_xe_pt2_df = pd.read_csv(pion_data_dir + 'pim_xe_pt2.dat', sep ='\s+')

#km
km_he_nu_df = pd.read_csv(kaon_data_dir + 'km_he_nu.dat', sep ='\s+')
km_ne_nu_df = pd.read_csv(kaon_data_dir + 'km_ne_nu.dat', sep ='\s+')
km_kr_nu_df = pd.read_csv(kaon_data_dir + 'km_kr_nu.dat', sep ='\s+')
km_xe_nu_df = pd.read_csv(kaon_data_dir + 'km_xe_nu.dat', sep ='\s+')

km_he_z_df = pd.read_csv(kaon_data_dir + 'km_he_z.dat', sep ='\s+')
km_ne_z_df = pd.read_csv(kaon_data_dir + 'km_ne_z.dat', sep ='\s+')
km_kr_z_df = pd.read_csv(kaon_data_dir + 'km_kr_z.dat', sep ='\s+')
km_xe_z_df = pd.read_csv(kaon_data_dir + 'km_xe_z.dat', sep ='\s+')

km_he_q2_df = pd.read_csv(kaon_data_dir + 'km_he_q2.dat', sep ='\s+')
km_ne_q2_df = pd.read_csv(kaon_data_dir + 'km_ne_q2.dat', sep ='\s+')
km_kr_q2_df = pd.read_csv(kaon_data_dir + 'km_kr_q2.dat', sep ='\s+')
km_xe_q2_df = pd.read_csv(kaon_data_dir + 'km_xe_q2.dat', sep ='\s+')

km_he_pt2_df = pd.read_csv(kaon_data_dir + 'km_he_pt2.dat', sep ='\s+')
km_ne_pt2_df = pd.read_csv(kaon_data_dir + 'km_ne_pt2.dat', sep ='\s+')
km_kr_pt2_df = pd.read_csv(kaon_data_dir + 'km_kr_pt2.dat', sep ='\s+')
km_xe_pt2_df = pd.read_csv(kaon_data_dir + 'km_xe_pt2.dat', sep ='\s+')


#pi0
pi0_he_nu_df = pd.read_csv(pion_data_dir + 'pi0_he_nu.dat', sep ='\s+')
pi0_ne_nu_df = pd.read_csv(pion_data_dir + 'pi0_ne_nu.dat', sep ='\s+')
pi0_kr_nu_df = pd.read_csv(pion_data_dir + 'pi0_kr_nu.dat', sep ='\s+')
pi0_xe_nu_df = pd.read_csv(pion_data_dir + 'pi0_xe_nu.dat', sep ='\s+')

pi0_he_z_df = pd.read_csv(pion_data_dir + 'pi0_he_z.dat', sep ='\s+')
pi0_ne_z_df = pd.read_csv(pion_data_dir + 'pi0_ne_z.dat', sep ='\s+')
pi0_kr_z_df = pd.read_csv(pion_data_dir + 'pi0_kr_z.dat', sep ='\s+')
pi0_xe_z_df = pd.read_csv(pion_data_dir + 'pi0_xe_z.dat', sep ='\s+')

pi0_he_q2_df = pd.read_csv(pion_data_dir + 'pi0_he_q2.dat', sep ='\s+')
pi0_ne_q2_df = pd.read_csv(pion_data_dir + 'pi0_ne_q2.dat', sep ='\s+')
pi0_kr_q2_df = pd.read_csv(pion_data_dir + 'pi0_kr_q2.dat', sep ='\s+')
pi0_xe_q2_df = pd.read_csv(pion_data_dir + 'pi0_xe_q2.dat', sep ='\s+')

pi0_he_pt2_df = pd.read_csv(pion_data_dir + 'pi0_he_pt2.dat', sep ='\s+')
pi0_ne_pt2_df = pd.read_csv(pion_data_dir + 'pi0_ne_pt2.dat', sep ='\s+')
pi0_kr_pt2_df = pd.read_csv(pion_data_dir + 'pi0_kr_pt2.dat', sep ='\s+')
pi0_xe_pt2_df = pd.read_csv(pion_data_dir + 'pi0_xe_pt2.dat', sep ='\s+')


#pip_ne_nu_df.head(10)


# In[14]:


class multdata:

    def __init__(self, file, var):
        file = file
        self.xdata  = file[var]
        self.var    = var
        self.nu     = file['Nu']
        self.z      = file['Z']
        self.Q2     = file['Q2']
        self.pt2    = file['pt2']
        self.mult   = file['MULT-RATIO']
        self.stat   = file['STAT']
        self.sys    = file['SYS']
        self.Ra     = file['DIS']
        self.error  = np.sqrt((file['STAT'])**2+ (file['SYS'])**2)


# In[35]:


class multplot:

    def __init__(self, data, nucleus, plotcolor):
        self.data = data
        self.nucleus = nucleus
        self.plotcolor = plotcolor

    def newplot(self,label):
        plt.title(label, fontsize = 'x-large')
        plt.axhline(1,color ='black',linewidth = 0.5)
        plt.ylabel(r"\rm $R^h_A$", fontsize = 'x-large',rotation = 0,labelpad = 15)
        plt.xlabel(self.data.var,fontsize = 'x-large', labelpad = 10)

    def plot(self):
        plt.errorbar(self.data.xdata, self.data.mult, self.data.error, fmt = "o",color = self.plotcolor, ecolor = self.plotcolor, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False,label = self.nucleus)
        plt.plot(self.data.xdata, self.data.Ra, color = self.plotcolor)
        plt.legend(loc = 'best')

    def endplot(self,ylow,yhigh):
        plt.ylim(ylow,yhigh)
        fig = plt.gcf()
        ax = plt.gca()
        ax.tick_params(axis = 'both', direction ='in',labelsize = 30)
        return fig, ax

    def newaxesplot(self,axes,ylow,yhigh):
        axes.axhline(1,color = 'black', linewidth = 1)
        axes.tick_params(axis = 'both', direction ='in',labelsize = 30)
        axes.set_ylim(ylow,yhigh)

    def axesplot(self,axes):
        axes.errorbar(self.data.xdata, self.data.mult, self.data.error, fmt = "o",color = self.plotcolor, ecolor = self.plotcolor, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False,label = self.nucleus)
        axes.plot(self.data.xdata, self.data.Ra, color = self.plotcolor)


# In[36]:


# Setup Arrays
#pip
pip_he_nu = multdata(pip_he_nu_df,'Nu')
pip_ne_nu = multdata(pip_ne_nu_df,'Nu')
pip_kr_nu = multdata(pip_kr_nu_df,'Nu')
pip_xe_nu = multdata(pip_xe_nu_df,'Nu')

pip_he_z = multdata(pip_he_z_df,'Z')
pip_ne_z = multdata(pip_ne_z_df,'Z')
pip_kr_z = multdata(pip_kr_z_df,'Z')
pip_xe_z = multdata(pip_xe_z_df,'Z')

pip_he_q2 = multdata(pip_he_q2_df,'Q2')
pip_ne_q2 = multdata(pip_ne_q2_df,'Q2')
pip_kr_q2 = multdata(pip_kr_q2_df,'Q2')
pip_xe_q2 = multdata(pip_xe_q2_df,'Q2')

pip_he_pt2 = multdata(pip_he_pt2_df,'pt2')
pip_ne_pt2 = multdata(pip_ne_pt2_df,'pt2')
pip_kr_pt2 = multdata(pip_kr_pt2_df,'pt2')
pip_xe_pt2 = multdata(pip_xe_pt2_df,'pt2')

#kp
kp_he_nu = multdata(kp_he_nu_df,'Nu')
kp_ne_nu = multdata(kp_ne_nu_df,'Nu')
kp_kr_nu = multdata(kp_kr_nu_df,'Nu')
kp_xe_nu = multdata(kp_xe_nu_df,'Nu')

kp_he_z = multdata(kp_he_z_df,'Z')
kp_ne_z = multdata(kp_ne_z_df,'Z')
kp_kr_z = multdata(kp_kr_z_df,'Z')
kp_xe_z = multdata(kp_xe_z_df,'Z')

kp_he_q2 = multdata(kp_he_q2_df,'Q2')
kp_ne_q2 = multdata(kp_ne_q2_df,'Q2')
kp_kr_q2 = multdata(kp_kr_q2_df,'Q2')
kp_xe_q2 = multdata(kp_xe_q2_df,'Q2')

kp_he_pt2 = multdata(kp_he_pt2_df,'pt2')
kp_ne_pt2 = multdata(kp_ne_pt2_df,'pt2')
kp_kr_pt2 = multdata(kp_kr_pt2_df,'pt2')
kp_xe_pt2 = multdata(kp_xe_pt2_df,'pt2')

#pim
pim_he_nu = multdata(pim_he_nu_df,'Nu')
pim_ne_nu = multdata(pim_ne_nu_df,'Nu')
pim_kr_nu = multdata(pim_kr_nu_df,'Nu')
pim_xe_nu = multdata(pim_xe_nu_df,'Nu')

pim_he_z = multdata(pim_he_z_df,'Z')
pim_ne_z = multdata(pim_ne_z_df,'Z')
pim_kr_z = multdata(pim_kr_z_df,'Z')
pim_xe_z = multdata(pim_xe_z_df,'Z')

pim_he_q2 = multdata(pim_he_q2_df,'Q2')
pim_ne_q2 = multdata(pim_ne_q2_df,'Q2')
pim_kr_q2 = multdata(pim_kr_q2_df,'Q2')
pim_xe_q2 = multdata(pim_xe_q2_df,'Q2')

pim_he_pt2 = multdata(pim_he_pt2_df,'pt2')
pim_ne_pt2 = multdata(pim_ne_pt2_df,'pt2')
pim_kr_pt2 = multdata(pim_kr_pt2_df,'pt2')
pim_xe_pt2 = multdata(pim_xe_pt2_df,'pt2')

#km
km_he_nu = multdata(km_he_nu_df,'Nu')
km_ne_nu = multdata(km_ne_nu_df,'Nu')
km_kr_nu = multdata(km_kr_nu_df,'Nu')
km_xe_nu = multdata(km_xe_nu_df,'Nu')

km_he_z = multdata(km_he_z_df,'Z')
km_ne_z = multdata(km_ne_z_df,'Z')
km_kr_z = multdata(km_kr_z_df,'Z')
km_xe_z = multdata(km_xe_z_df,'Z')

km_he_q2 = multdata(km_he_q2_df,'Q2')
km_ne_q2 = multdata(km_ne_q2_df,'Q2')
km_kr_q2 = multdata(km_kr_q2_df,'Q2')
km_xe_q2 = multdata(km_xe_q2_df,'Q2')

km_he_pt2 = multdata(km_he_pt2_df,'pt2')
km_ne_pt2 = multdata(km_ne_pt2_df,'pt2')
km_kr_pt2 = multdata(km_kr_pt2_df,'pt2')
km_xe_pt2 = multdata(km_xe_pt2_df,'pt2')

#pi0
pi0_he_nu = multdata(pi0_he_nu_df,'Nu')
pi0_ne_nu = multdata(pi0_ne_nu_df,'Nu')
pi0_kr_nu = multdata(pi0_kr_nu_df,'Nu')
pi0_xe_nu = multdata(pi0_xe_nu_df,'Nu')

pi0_he_z = multdata(pi0_he_z_df,'Z')
pi0_ne_z = multdata(pi0_ne_z_df,'Z')
pi0_kr_z = multdata(pi0_kr_z_df,'Z')
pi0_xe_z = multdata(pi0_xe_z_df,'Z')

pi0_he_q2 = multdata(pi0_he_q2_df,'Q2')
pi0_ne_q2 = multdata(pi0_ne_q2_df,'Q2')
pi0_kr_q2 = multdata(pi0_kr_q2_df,'Q2')
pi0_xe_q2 = multdata(pi0_xe_q2_df,'Q2')

pi0_he_pt2 = multdata(pi0_he_pt2_df,'pt2')
pi0_ne_pt2 = multdata(pi0_ne_pt2_df,'pt2')
pi0_kr_pt2 = multdata(pi0_kr_pt2_df,'pt2')
pi0_xe_pt2 = multdata(pi0_xe_pt2_df,'pt2')


# In[37]:


# Setup Plot Data
#pip
multplot_pip_he_nu = multplot(pip_he_nu, r'\rm He', 'red')
multplot_pip_ne_nu = multplot(pip_ne_nu, r'\rm Ne', 'black')
multplot_pip_kr_nu = multplot(pip_kr_nu, r'\rm Kr', 'green')
multplot_pip_xe_nu = multplot(pip_xe_nu, r'\rm Xe', 'blue')

multplot_pip_he_z = multplot(pip_he_z, r'\rm He', 'red')
multplot_pip_ne_z = multplot(pip_ne_z, r'\rm Ne', 'black')
multplot_pip_kr_z = multplot(pip_kr_z, r'\rm Kr', 'green')
multplot_pip_xe_z = multplot(pip_xe_z, r'\rm Xe', 'blue')

multplot_pip_he_q2 = multplot(pip_he_q2, r'\rm He', 'red')
multplot_pip_ne_q2 = multplot(pip_ne_q2, r'\rm Ne', 'black')
multplot_pip_kr_q2 = multplot(pip_kr_q2, r'\rm Kr', 'green')
multplot_pip_xe_q2 = multplot(pip_xe_q2, r'\rm Xe', 'blue')

multplot_pip_he_pt2 = multplot(pip_he_pt2, r'\rm He', 'red')
multplot_pip_ne_pt2 = multplot(pip_ne_pt2, r'\rm Ne', 'black')
multplot_pip_kr_pt2 = multplot(pip_kr_pt2, r'\rm Kr', 'green')
multplot_pip_xe_pt2 = multplot(pip_xe_pt2, r'\rm Xe', 'blue')

#kp
multplot_kp_he_nu = multplot(kp_he_nu, r'\rm He', 'red')
multplot_kp_ne_nu = multplot(kp_ne_nu, r'\rm Ne', 'black')
multplot_kp_kr_nu = multplot(kp_kr_nu, r'\rm Kr', 'green')
multplot_kp_xe_nu = multplot(kp_xe_nu, r'\rm Xe', 'blue')


multplot_kp_he_z = multplot(kp_he_z, r'\rm He', 'red')
multplot_kp_ne_z = multplot(kp_ne_z, r'\rm Ne', 'black')
multplot_kp_kr_z = multplot(kp_kr_z, r'\rm Kr', 'green')
multplot_kp_xe_z = multplot(kp_xe_z, r'\rm Xe', 'blue')

multplot_kp_he_q2 = multplot(kp_he_q2, r'\rm He', 'red')
multplot_kp_ne_q2 = multplot(kp_ne_q2, r'\rm Ne', 'black')
multplot_kp_kr_q2 = multplot(kp_kr_q2, r'\rm Kr', 'green')
multplot_kp_xe_q2 = multplot(kp_xe_q2, r'\rm Xe', 'blue')

multplot_kp_he_pt2 = multplot(kp_he_pt2, r'\rm He', 'red')
multplot_kp_ne_pt2 = multplot(kp_ne_pt2, r'\rm Ne', 'black')
multplot_kp_kr_pt2 = multplot(kp_kr_pt2, r'\rm Kr', 'green')
multplot_kp_xe_pt2 = multplot(kp_xe_pt2, r'\rm Xe', 'blue')

# Setup Plot Data
#pim
multplot_pim_he_nu = multplot(pim_he_nu, r'\rm He', 'red')
multplot_pim_ne_nu = multplot(pim_ne_nu, r'\rm Ne', 'black')
multplot_pim_kr_nu = multplot(pim_kr_nu, r'\rm Kr', 'green')
multplot_pim_xe_nu = multplot(pim_xe_nu, r'\rm Xe', 'blue')


multplot_pim_he_z = multplot(pim_he_z, r'\rm He', 'red')
multplot_pim_ne_z = multplot(pim_ne_z, r'\rm Ne', 'black')
multplot_pim_kr_z = multplot(pim_kr_z, r'\rm Kr', 'green')
multplot_pim_xe_z = multplot(pim_xe_z, r'\rm Xe', 'blue')

multplot_pim_he_q2 = multplot(pim_he_q2, r'\rm He', 'red')
multplot_pim_ne_q2 = multplot(pim_ne_q2, r'\rm Ne', 'black')
multplot_pim_kr_q2 = multplot(pim_kr_q2, r'\rm Kr', 'green')
multplot_pim_xe_q2 = multplot(pim_xe_q2, r'\rm Xe', 'blue')

multplot_pim_he_pt2 = multplot(pim_he_pt2, r'\rm He', 'red')
multplot_pim_ne_pt2 = multplot(pim_ne_pt2, r'\rm Ne', 'black')
multplot_pim_kr_pt2 = multplot(pim_kr_pt2, r'\rm Kr', 'green')
multplot_pim_xe_pt2 = multplot(pim_xe_pt2, r'\rm Xe', 'blue')

#km
multplot_km_he_nu = multplot(km_he_nu, r'\rm He', 'red')
multplot_km_ne_nu = multplot(km_ne_nu, r'\rm Ne', 'black')
multplot_km_kr_nu = multplot(km_kr_nu, r'\rm Kr', 'green')
multplot_km_xe_nu = multplot(km_xe_nu, r'\rm Xe', 'blue')


multplot_km_he_z = multplot(km_he_z, r'\rm He', 'red')
multplot_km_ne_z = multplot(km_ne_z, r'\rm Ne', 'black')
multplot_km_kr_z = multplot(km_kr_z, r'\rm Kr', 'green')
multplot_km_xe_z = multplot(km_xe_z, r'\rm Xe', 'blue')

multplot_km_he_q2 = multplot(km_he_q2, r'\rm He', 'red')
multplot_km_ne_q2 = multplot(km_ne_q2, r'\rm Ne', 'black')
multplot_km_kr_q2 = multplot(km_kr_q2, r'\rm Kr', 'green')
multplot_km_xe_q2 = multplot(km_xe_q2, r'\rm Xe', 'blue')

multplot_km_he_pt2 = multplot(km_he_pt2, r'\rm He', 'red')
multplot_km_ne_pt2 = multplot(km_ne_pt2, r'\rm Ne', 'black')
multplot_km_kr_pt2 = multplot(km_kr_pt2, r'\rm Kr', 'green')
multplot_km_xe_pt2 = multplot(km_xe_pt2, r'\rm Xe', 'blue')

#pi0
multplot_pi0_he_nu = multplot(pi0_he_nu, r'\rm He', 'red')
multplot_pi0_ne_nu = multplot(pi0_ne_nu, r'\rm Ne', 'black')
multplot_pi0_kr_nu = multplot(pi0_kr_nu, r'\rm Kr', 'green')
multplot_pi0_xe_nu = multplot(pi0_xe_nu, r'\rm Xe', 'blue')

multplot_pi0_he_z = multplot(pi0_he_z, r'\rm He', 'red')
multplot_pi0_ne_z = multplot(pi0_ne_z, r'\rm Ne', 'black')
multplot_pi0_kr_z = multplot(pi0_kr_z, r'\rm Kr', 'green')
multplot_pi0_xe_z = multplot(pi0_xe_z, r'\rm Xe', 'blue')

multplot_pi0_he_q2 = multplot(pi0_he_q2, r'\rm He', 'red')
multplot_pi0_ne_q2 = multplot(pi0_ne_q2, r'\rm Ne', 'black')
multplot_pi0_kr_q2 = multplot(pi0_kr_q2, r'\rm Kr', 'green')
multplot_pi0_xe_q2 = multplot(pi0_xe_q2, r'\rm Xe', 'blue')

multplot_pi0_he_pt2 = multplot(pi0_he_pt2, r'\rm He', 'red')
multplot_pi0_ne_pt2 = multplot(pi0_ne_pt2, r'\rm Ne', 'black')
multplot_pi0_kr_pt2 = multplot(pi0_kr_pt2, r'\rm Kr', 'green')
multplot_pi0_xe_pt2 = multplot(pi0_xe_pt2, r'\rm Xe', 'blue')


# In[38]:


fig, axs = plt.subplots(3, 4, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})
#pip
multplot_pip_he_nu.newaxesplot(axs[0,0],0.4,1.2)
multplot_pip_he_nu.axesplot(axs[0,0])
multplot_pip_ne_nu.axesplot(axs[0,0])
multplot_pip_kr_nu.axesplot(axs[0,0])
multplot_pip_xe_nu.axesplot(axs[0,0])


multplot_pip_he_z.newaxesplot(axs[0,1],0.4,1.2)
multplot_pip_he_z.axesplot(axs[0,1])
multplot_pip_ne_z.axesplot(axs[0,1])
multplot_pip_kr_z.axesplot(axs[0,1])
multplot_pip_xe_z.axesplot(axs[0,1])

multplot_pip_he_q2.newaxesplot(axs[0,2],0.4,1.4)
multplot_pip_he_q2.axesplot(axs[0,2])
multplot_pip_ne_q2.axesplot(axs[0,2])
multplot_pip_kr_q2.axesplot(axs[0,2])
multplot_pip_xe_q2.axesplot(axs[0,2])

multplot_pip_he_pt2.newaxesplot(axs[0,3],0.4,1.4)
multplot_pip_he_pt2.axesplot(axs[0,3])
multplot_pip_ne_pt2.axesplot(axs[0,3])
multplot_pip_xe_pt2.axesplot(axs[0,3])

#pi0
multplot_pi0_he_nu.newaxesplot(axs[1,0],0.4,1.2)
multplot_pi0_he_nu.axesplot(axs[1,0])
multplot_pi0_ne_nu.axesplot(axs[1,0])
multplot_pi0_kr_nu.axesplot(axs[1,0])
multplot_pi0_xe_nu.axesplot(axs[1,0])


multplot_pi0_he_z.newaxesplot(axs[1,1],0.4,1.2)
multplot_pi0_he_z.axesplot(axs[1,1])
multplot_pi0_ne_z.axesplot(axs[1,1])
multplot_pi0_kr_z.axesplot(axs[1,1])
multplot_pi0_xe_z.axesplot(axs[1,1])

multplot_pi0_he_q2.newaxesplot(axs[1,2],0.4,1.4)
multplot_pi0_he_q2.axesplot(axs[1,2])
multplot_pi0_ne_q2.axesplot(axs[1,2])
multplot_pi0_kr_q2.axesplot(axs[1,2])
multplot_pi0_xe_q2.axesplot(axs[1,2])

multplot_pi0_he_pt2.newaxesplot(axs[1,3],0.4,1.4)
multplot_pi0_he_pt2.axesplot(axs[1,3])
multplot_pi0_ne_pt2.axesplot(axs[1,3])
multplot_pi0_kr_pt2.axesplot(axs[1,3])
multplot_pi0_xe_pt2.axesplot(axs[1,3])

#pim
multplot_pim_he_nu.newaxesplot(axs[2,0],0.4,1.2)
multplot_pim_he_nu.axesplot(axs[2,0])
multplot_pim_ne_nu.axesplot(axs[2,0])
multplot_pim_kr_nu.axesplot(axs[2,0])
multplot_pim_xe_nu.axesplot(axs[2,0])


multplot_pim_he_z.newaxesplot(axs[2,1],0.4,1.2)
multplot_pim_he_z.axesplot(axs[2,1])
multplot_pim_ne_z.axesplot(axs[2,1])
multplot_pim_kr_z.axesplot(axs[2,1])
multplot_pim_xe_z.axesplot(axs[2,1])


multplot_pim_he_q2.newaxesplot(axs[2,2],0.4,1.4)
multplot_pim_he_q2.axesplot(axs[2,2])
multplot_pim_ne_q2.axesplot(axs[2,2])
multplot_pim_kr_q2.axesplot(axs[2,2])
multplot_pim_xe_q2.axesplot(axs[2,2])

multplot_pim_he_pt2.newaxesplot(axs[2,3],0.4,1.4)
multplot_pim_he_pt2.axesplot(axs[2,3])
multplot_pim_ne_pt2.axesplot(axs[2,3])
multplot_pim_kr_pt2.axesplot(axs[2,3])
multplot_pim_xe_pt2.axesplot(axs[2,3])


#axs[2,0].set_xlim(0,20)
#axs[2,1].set_xlim(0,1)
#axs[2,2].set_xlim(0,20)
#axs[2,3].set_xlim(0,10)

axs[0,0].set_ylabel(r"\rm $R^h_A$", fontsize = 30,rotation = 0,labelpad = 30)
axs[2,0].set_xlabel(r"\rm $\nu$ (GeV)", fontsize = 30, labelpad = 10)
axs[2,1].set_xlabel(r"\rm $z$", fontsize = 30, labelpad = 10)
axs[2,2].set_xlabel(r"\rm $Q^2$ (Ge$V^2$)", fontsize = 30, labelpad = 10)
axs[2,3].set_xlabel(r"\rm $p_t^2$ (Ge$V^2$)", fontsize = 30, labelpad = 10)

fig.set_size_inches(18,14)

# Text Boxes ##
props = dict(boxstyle='round', facecolor='white', alpha=0)

axs[0,0].text(0.7, 0.9,  r"\rm $\pi^+$", transform=axs[0,0].transAxes, fontsize=30,
    verticalalignment='top', bbox=props)

axs[1,0].text(0.7, 0.9,  r"\rm $\pi^0$", transform=axs[1,0].transAxes, fontsize=30,
    verticalalignment='top', bbox=props)

axs[2,0].text(0.7, 0.9,  r"\rm $\pi^-$", transform=axs[2,0].transAxes, fontsize=30,
    verticalalignment='top', bbox=props)

axs[2,3].legend(frameon = False,fontsize = 25,loc='upper center', bbox_to_anchor=(1.2, 0.9))


# In[34]:


fig.savefig('SIDIS_HERMES_PIONS.pdf', bbox_inches = "tight")


# In[39]:


fig2, axs2 = plt.subplots(2, 4, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})

#kp
multplot_kp_he_nu.newaxesplot(axs2[0,0],0.2,1.2)
multplot_kp_he_nu.axesplot(axs2[0,0])
multplot_kp_ne_nu.axesplot(axs2[0,0])
multplot_kp_kr_nu.axesplot(axs2[0,0])
multplot_kp_xe_nu.axesplot(axs2[0,0])

multplot_kp_he_z.newaxesplot(axs2[0,1],0.2,1.2)
multplot_kp_he_z.axesplot(axs2[0,1])
multplot_kp_ne_z.axesplot(axs2[0,1])
multplot_kp_kr_z.axesplot(axs2[0,1])
multplot_kp_xe_z.axesplot(axs2[0,1])


multplot_kp_he_q2.newaxesplot(axs2[0,2],0.2,1.4)
multplot_kp_he_q2.axesplot(axs2[0,2])
multplot_kp_ne_q2.axesplot(axs2[0,2])
multplot_kp_kr_q2.axesplot(axs2[0,2])
multplot_kp_xe_q2.axesplot(axs2[0,2])

multplot_kp_he_pt2.newaxesplot(axs2[0,3],0.2,1.4)
multplot_kp_he_pt2.axesplot(axs2[0,3])
multplot_kp_ne_pt2.axesplot(axs2[0,3])
multplot_kp_kr_pt2.axesplot(axs2[0,3])
multplot_kp_xe_pt2.axesplot(axs2[0,3])

#km
multplot_km_he_nu.newaxesplot(axs2[1,0],0.2,1.8)
multplot_km_he_nu.axesplot(axs2[1,0])
multplot_km_ne_nu.axesplot(axs2[1,0])
multplot_km_kr_nu.axesplot(axs2[1,0])
multplot_km_xe_nu.axesplot(axs2[1,0])

multplot_km_he_z.newaxesplot(axs2[1,1],0.4,1.8)
multplot_km_he_z.axesplot(axs2[1,1])
multplot_km_ne_z.axesplot(axs2[1,1])
multplot_km_kr_z.axesplot(axs2[1,1])
multplot_km_xe_z.axesplot(axs2[1,1])

multplot_km_he_q2.newaxesplot(axs2[1,2],0.2,1.8)
multplot_km_he_q2.axesplot(axs2[1,2])
multplot_km_ne_q2.axesplot(axs2[1,2])
multplot_km_kr_q2.axesplot(axs2[1,2])
multplot_km_xe_q2.axesplot(axs2[1,2])

multplot_km_he_pt2.newaxesplot(axs2[1,3],0.2,1.8)
multplot_km_he_pt2.axesplot(axs2[1,3])
multplot_km_ne_pt2.axesplot(axs2[1,3])
multplot_km_kr_pt2.axesplot(axs2[1,3])
multplot_km_xe_pt2.axesplot(axs2[1,3])

axs2[0,0].set_ylabel(r"\rm $R^h_A$", fontsize = 30,rotation = 0,labelpad = 30)
axs2[1,0].set_xlabel(r"\rm $\nu$ (GeV)", fontsize = 30, labelpad = 10)
axs2[1,1].set_xlabel(r"\rm $z$", fontsize = 30, labelpad = 10)
axs2[1,2].set_xlabel(r"\rm $Q^2$ (Ge$V^2$)", fontsize = 30, labelpad = 10)
axs2[1,3].set_xlabel(r"\rm $p_t^2$ (Ge$V^2$)", fontsize = 30, labelpad = 10)

fig2.set_size_inches(18,10)

# Text Boxes ##
props = dict(boxstyle='round', facecolor='white', alpha=0)

axs2[0,0].text(0.7, 0.9,  r"\rm $K^+$", transform=axs2[0,0].transAxes, fontsize=30,
    verticalalignment='top', bbox=props)

axs2[1,0].text(0.7, 0.9,  r"\rm $K^-$", transform=axs2[1,0].transAxes, fontsize=30,
    verticalalignment='top', bbox=props)

axs2[0,3].legend(frameon = False,fontsize = 25,loc='upper center', bbox_to_anchor=(1.2, 0.9))


# In[33]:


fig2.savefig('SIDIS_HERMES_KAONS.pdf', bbox_inches = "tight")


# In[ ]:
