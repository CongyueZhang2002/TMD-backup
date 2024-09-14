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
df = pd.read_csv("../plot_data/ATLAS5_Y_rat.dat", delim_whitespace=True)

# Plot
nrows = 1
ncols = 1
fig = py.figure(figsize=(ncols*4+1.5,nrows*4+1.5))

ax = py.subplot(111)

ax.errorbar(df.pt, df.CX/1000.,df.ERR/1000.,fmt = "r.")
ax.plot(df.pt, df.FUU/1000.)
ax.set_xlim(0,25)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $1/p_\perp\,d\sigma/d p_\perp$ (nb/GeV)',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)
ax.semilogy()
ax.text(0.60,0.9, r'\rm p+Pb $\rightarrow$ Z'  ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.8, r'\rm $\sqrt{s} = 5.02$ TeV' ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
ax.text(0.60,0.7, r'\rm $-3<y<2$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)

py.tight_layout(pad=-0.3)
py.savefig('ATLAS5_Y1.pdf',bbox_inches="tight")


# Import Pandas dataframe
df2 = pd.read_csv("../plot_data/ATLAS5_Y2.dat", delim_whitespace=True)
df3 = pd.read_csv("../plot_data/ATLAS5_Y3.dat", delim_whitespace=True)

# Plot
nrows = 1
ncols = 1
fig = py.figure(figsize=(ncols*4+1.5,nrows*4+1.5))

ax = py.subplot(111)

ax.errorbar(df2.pt, df2.CX/10000.,df2.ERR/10000.,fmt = "b.",label = r"\rm $-2<y<0\, (\times 10)$")
ax.plot(df2.pt, df2.FUU/10000.,label='_nolegend_',color = "blue")
ax.errorbar(df3.pt, df3.CX/1000. ,df3.ERR/1000. ,fmt = "r.",label = r"\rm $0<y<2$")
ax.plot(df3.pt, df3.FUU/1000. ,label='_nolegend_',color = "red")
ax.set_xlim(0,30)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $1/p_\perp\,d\sigma/d p_\perp$ (nb/GeV)',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)
ax.semilogy()
#ax.text(0.60,0.9, r'\rm p+Pb $\rightarrow$ Z'  ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#ax.text(0.60,0.8, r'\rm $\sqrt{s} = 5.02$ TeV' ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
#ax.text(0.60,0.7, r'\rm $-3<y<2$'              ,fontsize = 25,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)

ax.legend(frameon = False,fontsize = sizeOfFont,loc = 1)

py.tight_layout(pad=-0.3)
py.savefig('ATLAS5_Y23.pdf',bbox_inches="tight")
