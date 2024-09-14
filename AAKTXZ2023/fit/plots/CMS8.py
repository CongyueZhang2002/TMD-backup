# Import packages
import pandas as pd
import pylab as py
import numpy as np
import matplotlib

# Style settings
from  matplotlib import rc
from matplotlib import font_manager
from matplotlib.pyplot import gca
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
a = gca()
sizeOfFont = 20
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
fontProperties = {'family':'Helvetica','weight' : 'normal', 'size' : sizeOfFont}
ticks_font = matplotlib.font_manager.FontProperties(family='Helvetica', style='normal',size=sizeOfFont, weight='normal', stretch='normal')
rc('text',usetex=True)
rc('font',**fontProperties)

# Import Pandas dataframe
df = pd.read_csv("../plot_data/CMS8_DY_no.dat", delim_whitespace=True)

# Plot
nrows = 1
ncols = 1
fig = py.figure(figsize=(ncols*4+1.5,nrows*4+1.5))

ax = py.subplot(111)

ax.errorbar(df.pt, df.CX/1000.,df.ERR/1000.,fmt = "r.")
ax.plot(df.pt, df.FUU/1000.)
ax.set_xlim(0,25)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $d\sigma/d p_\perp$ (nb/GeV)',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)
ax.semilogy()
ax.text(0.60,0.9, r'\rm p+Pb $\rightarrow \gamma$ '  ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.8, r'\rm $\sqrt{s} = 8.16$ TeV' ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.7, r'\rm $-2.87<y<1.93$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.6, r'\rm $15<m_{\mu\mu}<60$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.set_ylim(1e-1,50)

py.tight_layout(pad=-0.3)
py.savefig('CMS8_DY_no.pdf',bbox_inches="tight")

# Import Pandas dataframe
df = pd.read_csv("../plot_data/CMS8_DY_fid.dat", delim_whitespace=True)

# Plot
nrows = 1
ncols = 1
fig = py.figure(figsize=(ncols*4+1.5,nrows*4+1.5))

ax = py.subplot(111)

ax.errorbar(df.pt, df.CX/1000.,df.ERR/1000.,fmt = "r.")
ax.plot(df.pt, df.FUU/1000.)
ax.set_xlim(0,25)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $d\sigma/d p_\perp$ (nb/GeV)',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)
ax.semilogy()
ax.text(0.60,0.9, r'\rm p+Pb $\rightarrow \gamma$, fiducial'  ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.8, r'\rm $\sqrt{s} = 8.16$ TeV' ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.7, r'\rm $-2.87<y<1.93$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.6, r'\rm $15<m_{\mu\mu}<60$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.set_ylim(5e-2,50)

py.tight_layout(pad=-0.3)
py.savefig('CMS8_DY_fid.pdf',bbox_inches="tight")

# Import Pandas dataframe
df = pd.read_csv("../plot_data/CMS8_ZZ_fid.dat", delim_whitespace=True)

# Plot
nrows = 1
ncols = 1
fig = py.figure(figsize=(ncols*4+1.5,nrows*4+1.5))

ax = py.subplot(111)

ax.errorbar(df.pt, df.CX/1000.,df.ERR/1000.,fmt = "r.")
ax.plot(df.pt, df.FUU/1000.)
ax.set_xlim(0,25)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $d\sigma/d p_\perp$ (nb/GeV)',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)
ax.semilogy()
ax.text(0.60,0.9, r'\rm p+Pb $\rightarrow \gamma/Z$, fiducial'  ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.8, r'\rm $\sqrt{s} = 8.16$ TeV' ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.7, r'\rm $-2.87<y<1.93$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.6, r'\rm $60<m_{\mu\mu}<120$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)

ax.set_xlim(0,45)

py.tight_layout(pad=-0.3)
py.savefig('CMS8_ZZ_fid.pdf',bbox_inches="tight")
