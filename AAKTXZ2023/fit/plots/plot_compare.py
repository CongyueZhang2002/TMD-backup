import numpy as np
import pylab as py
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib import rc
import matplotlib.font_manager
from matplotlib.pyplot import gca
from matplotlib.ticker import AutoMinorLocator
a = gca()
sizeOfFont = 40
#fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],'weight' : 'normal', 'size' : sizeOfFont}
fontProperties = {'sans-serif':['Helvetica'],'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = matplotlib.font_manager.FontProperties(family='Helvetica', style='normal',size=sizeOfFont, weight='normal', stretch='normal')
rc('text',usetex=True)
rc('font',**fontProperties)

normalize = True

#------------------------------------------------------------------
#      Load Theory
#------------------------------------------------------------------
df200_off = pd.read_csv("../plot_data/E288_200_off.dat", delim_whitespace=True)
df200_off_5 = df200_off
df200_off_5 = df200_off_5[df200_off_5.Q>4]
df200_off_5 = df200_off_5[df200_off_5.Q<5]
df200_off_5 = df200_off_5[df200_off_5.ds>0.]

df200_off_6 = df200_off
df200_off_6 = df200_off_6[df200_off_6.Q>5]
df200_off_6 = df200_off_6[df200_off_6.Q<6]
df200_off_6 = df200_off_6[df200_off_6.ds>0.]

df200_off_7 = df200_off
df200_off_7 = df200_off_7[df200_off_7.Q>6]
df200_off_7 = df200_off_7[df200_off_7.Q<7]
df200_off_7 = df200_off_7[df200_off_7.ds>0.]

df200_off_8 = df200_off
df200_off_8 = df200_off_8[df200_off_8.Q>7]
df200_off_8 = df200_off_8[df200_off_8.Q<8]
df200_off_8 = df200_off_8[df200_off_8.ds>0.]

df200_off_9 = df200_off
df200_off_9 = df200_off_9[df200_off_9.Q>8]
df200_off_9 = df200_off_9[df200_off_9.Q<9]
df200_off_9 = df200_off_9[df200_off_9.ds>0.]

#------------------------------------------------------------------
df300_off = pd.read_csv("../plot_data/E288_300_off.dat", delim_whitespace=True)
df300_off_5 = df300_off
df300_off_5 = df300_off_5[df300_off_5.Q>4]
df300_off_5 = df300_off_5[df300_off_5.Q<5]
df300_off_5 = df300_off_5[df300_off_5.ds>0.]

df300_off_6 = df300_off
df300_off_6 = df300_off_6[df300_off_6.Q>5]
df300_off_6 = df300_off_6[df300_off_6.Q<6]
df300_off_6 = df300_off_6[df300_off_6.ds>0.]

df300_off_7 = df300_off
df300_off_7 = df300_off_7[df300_off_7.Q>6]
df300_off_7 = df300_off_7[df300_off_7.Q<7]
df300_off_7 = df300_off_7[df300_off_7.ds>0.]

df300_off_8 = df300_off
df300_off_8 = df300_off_8[df300_off_8.Q>7]
df300_off_8 = df300_off_8[df300_off_8.Q<8]
df300_off_8 = df300_off_8[df300_off_8.ds>0.]

df300_off_9 = df300_off
df300_off_9 = df300_off_9[df300_off_9.Q>8]
df300_off_9 = df300_off_9[df300_off_9.Q<9]
df300_off_9 = df300_off_9[df300_off_9.ds>0.]

df300_off_12 = df300_off
df300_off_12 = df300_off_12[df300_off_12.Q>11]
df300_off_12 = df300_off_12[df300_off_12.Q<12]
df300_off_12 = df300_off_12[df300_off_12.ds>0.]

#------------------------------------------------------------------
df400_off = pd.read_csv("../plot_data/E288_400_off.dat", delim_whitespace=True)
df400_off_6 = df400_off
df400_off_6 = df400_off_6[df400_off_6.Q>5]
df400_off_6 = df400_off_6[df400_off_6.Q<6]
df400_off_6 = df400_off_6[df400_off_6.ds>0.]

df400_off_7 = df400_off
df400_off_7 = df400_off_7[df400_off_7.Q>6]
df400_off_7 = df400_off_7[df400_off_7.Q<7]
df400_off_7 = df400_off_7[df400_off_7.ds>0.]

df400_off_8 = df400_off
df400_off_8 = df400_off_8[df400_off_8.Q>7]
df400_off_8 = df400_off_8[df400_off_8.Q<8]
df400_off_8 = df400_off_8[df400_off_8.ds>0.]

df400_off_9 = df400_off
df400_off_9 = df400_off_9[df400_off_9.Q>8]
df400_off_9 = df400_off_9[df400_off_9.Q<9]
df400_off_9 = df400_off_9[df400_off_9.ds>0.]

df400_off_12 = df400_off
df400_off_12 = df400_off_12[df400_off_12.Q>11]
df400_off_12 = df400_off_12[df400_off_12.Q<12]
df400_off_12 = df400_off_12[df400_off_12.ds>0.]

df400_off_13 = df400_off
df400_off_13 = df400_off_13[df400_off_13.Q>12]
df400_off_13 = df400_off_13[df400_off_13.Q<13]
df400_off_13 = df400_off_13[df400_off_13.ds>0.]

df400_off_14 = df400_off
df400_off_14 = df400_off_14[df400_off_14.Q>13]
df400_off_14 = df400_off_14[df400_off_14.Q<14]
df400_off_14 = df400_off_14[df400_off_14.ds>0.]



#------------------------------------------------------------------
#      Load Theory
#------------------------------------------------------------------
df200 = pd.read_csv("../plot_data/E288_200.dat", delim_whitespace=True)
df200_5 = df200
df200_5 = df200_5[df200_5.Q>4]
df200_5 = df200_5[df200_5.Q<5]
df200_5 = df200_5[df200_5.ds>0.]

df200_6 = df200
df200_6 = df200_6[df200_6.Q>5]
df200_6 = df200_6[df200_6.Q<6]
df200_6 = df200_6[df200_6.ds>0.]

df200_7 = df200
df200_7 = df200_7[df200_7.Q>6]
df200_7 = df200_7[df200_7.Q<7]
df200_7 = df200_7[df200_7.ds>0.]

df200_8 = df200
df200_8 = df200_8[df200_8.Q>7]
df200_8 = df200_8[df200_8.Q<8]
df200_8 = df200_8[df200_8.ds>0.]

df200_9 = df200
df200_9 = df200_9[df200_9.Q>8]
df200_9 = df200_9[df200_9.Q<9]
df200_9 = df200_9[df200_9.ds>0.]

#------------------------------------------------------------------
df300 = pd.read_csv("../plot_data/E288_300.dat", delim_whitespace=True)
df300_5 = df300
df300_5 = df300_5[df300_5.Q>4]
df300_5 = df300_5[df300_5.Q<5]
df300_5 = df300_5[df300_5.ds>0.]

df300_6 = df300
df300_6 = df300_6[df300_6.Q>5]
df300_6 = df300_6[df300_6.Q<6]
df300_6 = df300_6[df300_6.ds>0.]

df300_7 = df300
df300_7 = df300_7[df300_7.Q>6]
df300_7 = df300_7[df300_7.Q<7]
df300_7 = df300_7[df300_7.ds>0.]

df300_8 = df300
df300_8 = df300_8[df300_8.Q>7]
df300_8 = df300_8[df300_8.Q<8]
df300_8 = df300_8[df300_8.ds>0.]

df300_9 = df300
df300_9 = df300_9[df300_9.Q>8]
df300_9 = df300_9[df300_9.Q<9]
df300_9 = df300_9[df300_9.ds>0.]

df300_12 = df300
df300_12 = df300_12[df300_12.Q>11]
df300_12 = df300_12[df300_12.Q<12]
df300_12 = df300_12[df300_12.ds>0.]

#------------------------------------------------------------------
df400 = pd.read_csv("../plot_data/E288_400.dat", delim_whitespace=True)
df400_6 = df400
df400_6 = df400_6[df400_6.Q>5]
df400_6 = df400_6[df400_6.Q<6]
df400_6 = df400_6[df400_6.ds>0.]

df400_7 = df400
df400_7 = df400_7[df400_7.Q>6]
df400_7 = df400_7[df400_7.Q<7]
df400_7 = df400_7[df400_7.ds>0.]

df400_8 = df400
df400_8 = df400_8[df400_8.Q>7]
df400_8 = df400_8[df400_8.Q<8]
df400_8 = df400_8[df400_8.ds>0.]

df400_9 = df400
df400_9 = df400_9[df400_9.Q>8]
df400_9 = df400_9[df400_9.Q<9]
df400_9 = df400_9[df400_9.ds>0.]

df400_12 = df400
df400_12 = df400_12[df400_12.Q>11]
df400_12 = df400_12[df400_12.Q<12]
df400_12 = df400_12[df400_12.ds>0.]

df400_13 = df400
df400_13 = df400_13[df400_13.Q>12]
df400_13 = df400_13[df400_13.Q<13]
df400_13 = df400_13[df400_13.ds>0.]

df400_14 = df400
df400_14 = df400_14[df400_14.Q>13]
df400_14 = df400_14[df400_14.Q<14]
df400_14 = df400_14[df400_14.ds>0.]

#--------------------------------------------------------------------
#         Load Experimental Data
#--------------------------------------------------------------------
df200exp = pd.read_csv("../expdata/E288_200.csv", delim_whitespace=True)
df200exp_5 = df200exp
df200exp_5 = df200exp_5[df200exp_5.Q>4]
df200exp_5 = df200exp_5[df200exp_5.Q<5]
df200exp_5 = df200exp_5[df200exp_5.value>0.]

df200exp_6 = df200exp
df200exp_6 = df200exp_6[df200exp_6.Q>5]
df200exp_6 = df200exp_6[df200exp_6.Q<6]
df200exp_6 = df200exp_6[df200exp_6.value>0.]

df200exp_7 = df200exp
df200exp_7 = df200exp_7[df200exp_7.Q>6]
df200exp_7 = df200exp_7[df200exp_7.Q<7]
df200exp_7 = df200exp_7[df200exp_7.value>0.]

df200exp_8 = df200exp
df200exp_8 = df200exp_8[df200exp_8.Q>7]
df200exp_8 = df200exp_8[df200exp_8.Q<8]
df200exp_8 = df200exp_8[df200exp_8.value>0.]

df200exp_9 = df200exp
df200exp_9 = df200exp_9[df200exp_9.Q>8]
df200exp_9 = df200exp_9[df200exp_9.Q<9]
df200exp_9 = df200exp_9[df200exp_9.value>0.]

#------------------------------------------------------------------
df300exp = pd.read_csv("../expdata/E288_300.csv", delim_whitespace=True)
df300exp_5 = df300exp
df300exp_5 = df300exp_5[df300exp_5.Q>4]
df300exp_5 = df300exp_5[df300exp_5.Q<5]
df300exp_5 = df300exp_5[df300exp_5.value>0.]

df300exp_6 = df300exp
df300exp_6 = df300exp_6[df300exp_6.Q>5]
df300exp_6 = df300exp_6[df300exp_6.Q<6]
df300exp_6 = df300exp_6[df300exp_6.value>0.]

df300exp_7 = df300exp
df300exp_7 = df300exp_7[df300exp_7.Q>6]
df300exp_7 = df300exp_7[df300exp_7.Q<7]
df300exp_7 = df300exp_7[df300exp_7.value>0.]

df300exp_8 = df300exp
df300exp_8 = df300exp_8[df300exp_8.Q>7]
df300exp_8 = df300exp_8[df300exp_8.Q<8]
df300exp_8 = df300exp_8[df300exp_8.value>0.]

df300exp_9 = df300exp
df300exp_9 = df300exp_9[df300exp_9.Q>8]
df300exp_9 = df300exp_9[df300exp_9.Q<9]
df300exp_9 = df300exp_9[df300exp_9.value>0.]

df300exp_12 = df300exp
df300exp_12 = df300exp_12[df300exp_12.Q>11]
df300exp_12 = df300exp_12[df300exp_12.Q<12]
df300exp_12 = df300exp_12[df300exp_12.value>0.]


#------------------------------------------------------------------
df400exp = pd.read_csv("../expdata/E288_400.csv", delim_whitespace=True)
df400exp_6 = df400exp
df400exp_6 = df400exp_6[df400exp_6.Q>5]
df400exp_6 = df400exp_6[df400exp_6.Q<6]
df400exp_6 = df400exp_6[df400exp_6.value>0.]

df400exp_7 = df400exp
df400exp_7 = df400exp_7[df400exp_7.Q>6]
df400exp_7 = df400exp_7[df400exp_7.Q<7]
df400exp_7 = df400exp_7[df400exp_7.value>0.]

df400exp_8 = df400exp
df400exp_8 = df400exp_8[df400exp_8.Q>7]
df400exp_8 = df400exp_8[df400exp_8.Q<8]
df400exp_8 = df400exp_8[df400exp_8.value>0.]

df400exp_9 = df400exp
df400exp_9 = df400exp_9[df400exp_9.Q>8]
df400exp_9 = df400exp_9[df400exp_9.Q<9]
df400exp_9 = df400exp_9[df400exp_9.value>0.]

df400exp_12 = df400exp
df400exp_12 = df400exp_12[df400exp_12.Q>11]
df400exp_12 = df400exp_12[df400exp_12.Q<12]
df400exp_12 = df400exp_12[df400exp_12.value>0.]

df400exp_13 = df400exp
df400exp_13 = df400exp_13[df400exp_13.Q>12]
df400exp_13 = df400exp_13[df400exp_13.Q<13]
df400exp_13 = df400exp_13[df400exp_13.value>0.]

df400exp_14 = df400exp
df400exp_14 = df400exp_14[df400exp_14.Q>13]
df400exp_14 = df400exp_14[df400exp_14.Q<14]
df400exp_14 = df400exp_14[df400exp_14.value>0.]

#--------------------------------------------------------------------
#         Generate Plot
#--------------------------------------------------------------------

nrows = 1
ncols = 3
_labelsize = 2.5
fig = py.figure(figsize = (5.*ncols+_labelsize,3.5*nrows+_labelsize))
ax = py.subplot(131)

if normalize == True:
    norm = 1./df200_off_5.ds.iloc[0]*df200exp_5.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_off_5.pt/df200_off_5.Q),np.max(df200_off_5.pt/df200_off_5.Q),100)
fxnew = interp1d(df200_off_5.pt/df200_off_5.Q,2*np.pi*df200_off_5.pt*np.abs(df200_off_5.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'r')

if normalize == True:
    norm = 1./df200_off_6.ds.iloc[0]*df200exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_off_6.pt/df200_off_6.Q),np.max(df200_off_6.pt/df200_off_6.Q),100)
fxnew = interp1d(df200_off_6.pt/df200_off_6.Q,2*np.pi*df200_off_6.pt*np.abs(df200_off_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'm')

if normalize == True:
    norm = 1./df200_off_7.ds.iloc[0]*df200exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_off_7.pt/df200_off_7.Q),np.max(df200_off_7.pt/df200_off_7.Q),100)
fxnew = interp1d(df200_off_7.pt/df200_off_7.Q,2*np.pi*df200_off_7.pt*np.abs(df200_off_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'y')

if normalize == True:
    norm = 1./df200_off_8.ds.iloc[0]*df200exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_off_8.pt/df200_off_8.Q),np.max(df200_off_8.pt/df200_off_8.Q),100)
fxnew = interp1d(df200_off_8.pt/df200_off_8.Q,2*np.pi*df200_off_8.pt*np.abs(df200_off_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'g')

if normalize == True:
    norm = 1./df200_off_9.ds.iloc[0]*df200exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_off_9.pt/df200_off_9.Q),np.max(df200_off_9.pt/df200_off_9.Q),100)
fxnew = interp1d(df200_off_9.pt/df200_off_9.Q,2*np.pi*df200_off_9.pt*np.abs(df200_off_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'b')

if normalize == True:
    norm = 1./df200_5.ds.iloc[0]*df200exp_5.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_5.pt/df200_5.Q),np.max(df200_5.pt/df200_5.Q),100)
fxnew = interp1d(df200_5.pt/df200_5.Q,2*np.pi*df200_5.pt*np.abs(df200_5.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'r')

if normalize == True:
    norm = 1./df200_6.ds.iloc[0]*df200exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_6.pt/df200_6.Q),np.max(df200_6.pt/df200_6.Q),100)
fxnew = interp1d(df200_6.pt/df200_6.Q,2*np.pi*df200_6.pt*np.abs(df200_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'm')

if normalize == True:
    norm = 1./df200_7.ds.iloc[0]*df200exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_7.pt/df200_7.Q),np.max(df200_7.pt/df200_7.Q),100)
fxnew = interp1d(df200_7.pt/df200_7.Q,2*np.pi*df200_7.pt*np.abs(df200_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'y')

if normalize == True:
    norm = 1./df200_8.ds.iloc[0]*df200exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_8.pt/df200_8.Q),np.max(df200_8.pt/df200_8.Q),100)
fxnew = interp1d(df200_8.pt/df200_8.Q,2*np.pi*df200_8.pt*np.abs(df200_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'g')

if normalize == True:
    norm = 1./df200_9.ds.iloc[0]*df200exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df200_9.pt/df200_9.Q),np.max(df200_9.pt/df200_9.Q),100)
fxnew = interp1d(df200_9.pt/df200_9.Q,2*np.pi*df200_9.pt*np.abs(df200_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'b')

ax.errorbar(df200exp_5.pt/df200exp_5.Q,2*np.pi*df200exp_5.pt*df200exp_5.value/1000.,df200exp_5.error/1000.,fmt = 'ro',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df200exp_6.pt/df200exp_6.Q,2*np.pi*df200exp_6.pt*df200exp_6.value/1000.,df200exp_6.error/1000.,fmt = 'mo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df200exp_7.pt/df200exp_7.Q,2*np.pi*df200exp_7.pt*df200exp_7.value/1000.,df200exp_7.error/1000.,fmt = 'yo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df200exp_8.pt/df200exp_8.Q,2*np.pi*df200exp_8.pt*df200exp_8.value/1000.,df200exp_8.error/1000.,fmt = 'go',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df200exp_9.pt/df200exp_9.Q,2*np.pi*df200exp_9.pt*df200exp_9.value/1000.,df200exp_9.error/1000.,fmt = 'bo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.set_title(r'\rm E288(200 GeV)',fontsize = '40',y = 1.025)
ax.semilogy()
ax.set_xlim(0,0.2)
ax.set_ylim(5e-4,3e1)
ax.tick_params('both', length=15, width=1, which='major',direction='in',labelsize = 40)
ax.tick_params('both', length=7 , width=1, which='minor',direction='in')
ax.xaxis.set_major_locator(py.MaxNLocator(nbins = 4,prune = 'upper'))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax = py.subplot(132)

if normalize == True:
    norm = 1./df300_off_5.ds.iloc[0]*df300exp_5.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_5.pt/df300_off_5.Q),np.max(df300_off_5.pt/df300_off_5.Q),100)
fxnew = interp1d(df300_off_5.pt/df300_off_5.Q,2*np.pi*df300_off_5.pt*np.abs(df300_off_5.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'r')

if normalize == True:
    norm = 1./df300_off_6.ds.iloc[0]*df300exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_6.pt/df300_off_6.Q),np.max(df300_off_6.pt/df300_off_6.Q),100)
fxnew = interp1d(df300_off_6.pt/df300_off_6.Q,2*np.pi*df300_off_6.pt*np.abs(df300_off_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'm')

if normalize == True:
    norm = 1./df300_off_7.ds.iloc[0]*df300exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_7.pt/df300_off_7.Q),np.max(df300_off_7.pt/df300_off_7.Q),100)
fxnew = interp1d(df300_off_7.pt/df300_off_7.Q,2*np.pi*df300_off_7.pt*np.abs(df300_off_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'y')

if normalize == True:
    norm = 1./df300_off_8.ds.iloc[0]*df300exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_8.pt/df300_off_8.Q),np.max(df300_off_8.pt/df300_off_8.Q),100)
fxnew = interp1d(df300_off_8.pt/df300_off_8.Q,2*np.pi*df300_off_8.pt*np.abs(df300_off_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'g')

if normalize == True:
    norm = 1./df300_off_9.ds.iloc[0]*df300exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_9.pt/df300_off_9.Q),np.max(df300_off_9.pt/df300_off_9.Q),100)
fxnew = interp1d(df300_off_9.pt/df300_off_9.Q,2*np.pi*df300_off_9.pt*np.abs(df300_off_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'b')

if normalize == True:
    norm = 1./df300_off_12.ds.iloc[0]*df300exp_12.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_off_12.pt/df300_off_12.Q),np.max(df300_off_12.pt/df300_off_12.Q),100)
fxnew = interp1d(df300_off_12.pt/df300_off_12.Q,2*np.pi*df300_off_12.pt*np.abs(df300_off_12.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'gray')


if normalize == True:
    norm = 1./df300_5.ds.iloc[0]*df300exp_5.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_5.pt/df300_5.Q),np.max(df300_5.pt/df300_5.Q),100)
fxnew = interp1d(df300_5.pt/df300_5.Q,2*np.pi*df300_5.pt*np.abs(df300_5.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'r')

if normalize == True:
    norm = 1./df300_6.ds.iloc[0]*df300exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_6.pt/df300_6.Q),np.max(df300_6.pt/df300_6.Q),100)
fxnew = interp1d(df300_6.pt/df300_6.Q,2*np.pi*df300_6.pt*np.abs(df300_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'm')

if normalize == True:
    norm = 1./df300_7.ds.iloc[0]*df300exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_7.pt/df300_7.Q),np.max(df300_7.pt/df300_7.Q),100)
fxnew = interp1d(df300_7.pt/df300_7.Q,2*np.pi*df300_7.pt*np.abs(df300_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'y')

if normalize == True:
    norm = 1./df300_8.ds.iloc[0]*df300exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_8.pt/df300_8.Q),np.max(df300_8.pt/df300_8.Q),100)
fxnew = interp1d(df300_8.pt/df300_8.Q,2*np.pi*df300_8.pt*np.abs(df300_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'g')

if normalize == True:
    norm = 1./df300_9.ds.iloc[0]*df300exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_9.pt/df300_9.Q),np.max(df300_9.pt/df300_9.Q),100)
fxnew = interp1d(df300_9.pt/df300_9.Q,2*np.pi*df300_9.pt*np.abs(df300_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'b')

if normalize == True:
    norm = 1./df300_12.ds.iloc[0]*df300exp_12.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df300_12.pt/df300_12.Q),np.max(df300_12.pt/df300_12.Q),100)
fxnew = interp1d(df300_12.pt/df300_12.Q,2*np.pi*df300_12.pt*np.abs(df300_12.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'gray')

ax.errorbar(df300exp_5.pt/df300exp_5.Q,2*np.pi*df300exp_5.pt*df300exp_5.value/1000.,df300exp_5.error/1000.,fmt = 'ro',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df300exp_6.pt/df300exp_6.Q,2*np.pi*df300exp_6.pt*df300exp_6.value/1000.,df300exp_6.error/1000.,fmt = 'mo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df300exp_7.pt/df300exp_7.Q,2*np.pi*df300exp_7.pt*df300exp_7.value/1000.,df300exp_7.error/1000.,fmt = 'yo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df300exp_8.pt/df300exp_8.Q,2*np.pi*df300exp_8.pt*df300exp_8.value/1000.,df300exp_8.error/1000.,fmt = 'go',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df300exp_9.pt/df300exp_9.Q,2*np.pi*df300exp_9.pt*df300exp_9.value/1000.,df300exp_9.error/1000.,fmt = 'bo',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df300exp_12.pt/df300exp_12.Q,2*np.pi*df300exp_12.pt*df300exp_12.value/1000.,df300exp_12.error/1000.,fmt = 'o',color = 'gray',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.set_title(r'\rm E288(300 GeV)',fontsize = '40',y = 1.025)
ax.semilogy()
ax.set_xlim(0,0.2)
ax.set_ylim(5e-4,3e1)
ax.set_yticks([])
#py.show()
ax.tick_params('both', length=15, width=1, which='major',direction='in',labelsize = 40)
ax.tick_params('both', length=7 , width=1, which='minor',direction='in')
ax.xaxis.set_major_locator(py.MaxNLocator(nbins = 4,prune = 'upper'))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax = py.subplot(133)

if normalize == True:
    norm = 1./df400_off_6.ds.iloc[0]*df400exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_6.pt/df400_off_6.Q),np.max(df400_off_6.pt/df400_off_6.Q),100)
fxnew = interp1d(df400_off_6.pt/df400_off_6.Q,2*np.pi*df400_off_6.pt*np.abs(df400_off_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'm',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_7.ds.iloc[0]*df400exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_7.pt/df400_off_7.Q),np.max(df400_off_7.pt/df400_off_7.Q),100)
fxnew = interp1d(df400_off_7.pt/df400_off_7.Q,2*np.pi*df400_off_7.pt*np.abs(df400_off_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'y',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_8.ds.iloc[0]*df400exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_8.pt/df400_off_8.Q),np.max(df400_off_8.pt/df400_off_8.Q),100)
fxnew = interp1d(df400_off_8.pt/df400_off_8.Q,2*np.pi*df400_off_8.pt*np.abs(df400_off_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'g',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_9.ds.iloc[0]*df400exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_9.pt/df400_off_9.Q),np.max(df400_off_9.pt/df400_off_9.Q),100)
fxnew = interp1d(df400_off_9.pt/df400_off_9.Q,2*np.pi*df400_off_9.pt*np.abs(df400_off_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'b',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_12.ds.iloc[0]*df400exp_12.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_12.pt/df400_off_12.Q),np.max(df400_off_12.pt/df400_off_12.Q),100)
fxnew = interp1d(df400_off_12.pt/df400_off_12.Q,2*np.pi*df400_off_12.pt*np.abs(df400_off_12.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'gray',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_13.ds.iloc[0]*df400exp_13.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_13.pt/df400_off_13.Q),np.max(df400_off_13.pt/df400_off_13.Q),100)
fxnew = interp1d(df400_off_13.pt/df400_off_13.Q,2*np.pi*df400_off_13.pt*np.abs(df400_off_13.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'c',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_off_14.ds.iloc[0]*df400exp_14.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_off_14.pt/df400_off_14.Q),np.max(df400_off_14.pt/df400_off_14.Q),100)
fxnew = interp1d(df400_off_14.pt/df400_off_14.Q,2*np.pi*df400_off_14.pt*np.abs(df400_off_14.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, linestyle = 'dashed', color = 'k',label = '_nolegend_')


if normalize == True:
    norm = 1./df400_6.ds.iloc[0]*df400exp_6.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_6.pt/df400_6.Q),np.max(df400_6.pt/df400_6.Q),100)
fxnew = interp1d(df400_6.pt/df400_6.Q,2*np.pi*df400_6.pt*np.abs(df400_6.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'm',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_7.ds.iloc[0]*df400exp_7.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_7.pt/df400_7.Q),np.max(df400_7.pt/df400_7.Q),100)
fxnew = interp1d(df400_7.pt/df400_7.Q,2*np.pi*df400_7.pt*np.abs(df400_7.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'y',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_8.ds.iloc[0]*df400exp_8.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_8.pt/df400_8.Q),np.max(df400_8.pt/df400_8.Q),100)
fxnew = interp1d(df400_8.pt/df400_8.Q,2*np.pi*df400_8.pt*np.abs(df400_8.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'g',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_9.ds.iloc[0]*df400exp_9.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_9.pt/df400_9.Q),np.max(df400_9.pt/df400_9.Q),100)
fxnew = interp1d(df400_9.pt/df400_9.Q,2*np.pi*df400_9.pt*np.abs(df400_9.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'b',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_12.ds.iloc[0]*df400exp_12.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_12.pt/df400_12.Q),np.max(df400_12.pt/df400_12.Q),100)
fxnew = interp1d(df400_12.pt/df400_12.Q,2*np.pi*df400_12.pt*np.abs(df400_12.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'gray',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_13.ds.iloc[0]*df400exp_13.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_13.pt/df400_13.Q),np.max(df400_13.pt/df400_13.Q),100)
fxnew = interp1d(df400_13.pt/df400_13.Q,2*np.pi*df400_13.pt*np.abs(df400_13.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'c',label = '_nolegend_')

if normalize == True:
    norm = 1./df400_14.ds.iloc[0]*df400exp_14.value.iloc[0]/1000.
else:
    norm = 1
print(norm)
xnew  = np.linspace(np.min(df400_14.pt/df400_14.Q),np.max(df400_14.pt/df400_14.Q),100)
fxnew = interp1d(df400_14.pt/df400_14.Q,2*np.pi*df400_14.pt*np.abs(df400_14.ds)*norm,kind = 3)
ax.plot(xnew,fxnew(xnew),linewidth = 3, zorder = 1, color = 'k',label = '_nolegend_')

ax.errorbar(df300exp_5.pt/df300exp_5.Q,100*2*np.pi*df300exp_5.pt*df300exp_5.value/1000.,df300exp_5.error/1000.,fmt = 'ro' ,label = r'\rm $4<Q<5$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_6.pt/df400exp_6.Q,2*np.pi*df400exp_6.pt*df400exp_6.value/1000.,df400exp_6.error/1000.,fmt = 'mo'     ,label = r'\rm $5<Q<6$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_7.pt/df400exp_7.Q,2*np.pi*df400exp_7.pt*df400exp_7.value/1000.,df400exp_7.error/1000.,fmt = 'yo'     ,label = r'\rm $6<Q<7$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_8.pt/df400exp_8.Q,2*np.pi*df400exp_8.pt*df400exp_8.value/1000.,df400exp_8.error/1000.,fmt = 'go'     ,label = r'\rm $7<Q<8$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_9.pt/df400exp_9.Q,2*np.pi*df400exp_9.pt*df400exp_9.value/1000.,df400exp_9.error/1000.,fmt = 'bo'     ,label = r'\rm $8<Q<9$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_12.pt/df400exp_12.Q,2*np.pi*df400exp_12.pt*df400exp_12.value/1000.,df400exp_12.error/1000.,fmt = 'o', color = 'gray',label = r'\rm $11<Q<12$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_13.pt/df400exp_13.Q,2*np.pi*df400exp_13.pt*df400exp_13.value/1000.,df400exp_13.error/1000.,fmt = 'co',label = r'\rm $12<Q<13$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.errorbar(df400exp_14.pt/df400exp_14.Q,2*np.pi*df400exp_14.pt*df400exp_14.value/1000.,df400exp_14.error/1000.,fmt = 'ko',label = r'\rm $13<Q<14$',markersize = 10, capsize = 5,elinewidth = 2,capthick = 2, zorder = 2)
ax.set_title(r'\rm E288(400 GeV)',fontsize = '40',y = 1.025)
ax.semilogy()
ax.set_xlim(0,0.2)
ax.set_ylim(5e-4,3e1)
ax.set_yticks([])
ax.tick_params('both', length=15, width=1, which='major',direction='in',labelsize = 40)
ax.tick_params('both', length=7 , width=1, which='minor',direction='in')
ax.xaxis.set_major_locator(py.MaxNLocator(nbins = 4,prune = 'upper'))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.legend(loc='center left',frameon = False,fontsize = 40, bbox_to_anchor=(1., 0.5))

fig.text(0.5 ,-0.05, r'\rm $q_{\perp}/Q$', ha='center',fontsize = 40)
fig.text(0.02,0.5  , r'\rm $d\sigma/d q_{\perp}$', va='center', rotation='vertical',fontsize = 40)
#py.tight_layout()
py.subplots_adjust(wspace=0, hspace=0)
fig.savefig('Drell-Yan_E288_normal.pdf', bbox_inches = "tight")
