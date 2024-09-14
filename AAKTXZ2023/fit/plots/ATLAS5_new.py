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

ax.errorbar(df.pt, df.CX ,df.ERR ,fmt = "r.")
ax.plot(df.pt, df.FUU ,label='_nolegend_',color = "red")
ax.set_xlim(0,30)
ax.set_ylim(0.6,1.1)
ax.tick_params(axis="both",direction="in", length=10)


ax.set_ylabel(r'\rm $d\sigma/d p_\perp (y>0)/d\sigma/d p_\perp (y<0)$',fontsize=sizeOfFont)
ax.set_xlabel(r'\rm $p_\perp$ (GeV)',fontsize=sizeOfFont)

ax.legend(frameon = False,fontsize = sizeOfFont,loc = 1)

py.tight_layout(pad=-0.3)
py.savefig('ATLAS5_RAT.pdf',bbox_inches="tight")
