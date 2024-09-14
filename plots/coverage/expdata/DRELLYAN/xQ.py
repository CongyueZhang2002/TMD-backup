import pandas as pd
import numpy as np
import os
import pylab as py
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib
sizeOfFont = 20
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
#ticks_font = matplotlib.font_manager.FontProperties(family='Helvetica', style='normal',size=sizeOfFont, weight='normal', stretch='normal')

dftot = pd.DataFrame(columns=["Nu","Z","Q2","pt2","MULT-RATIO","STAT","SYS","fuu","fuua","DIS"])

#CMS 5
SCMS5 = 5020.

yminCMS5 = -2.8
ymaxCMS5 =  2.0
QminCMS5 =  60.
QmaxCMS5 = 120.
QbarCMS5 = 0.5*(QminCMS5+QmaxCMS5)

xminCMS5 = QbarCMS5/SCMS5*np.exp(yminCMS5)
xmaxCMS5 = QbarCMS5/SCMS5*np.exp(ymaxCMS5)

#ATLAS 5
ATLS5 = 5020.

yminATL5 = -3.0
ymaxATL5 =  2.0
QminATL5 =  66.
QmaxATL5 = 116.
QbarATL5 = 0.5*(QminATL5+QmaxATL5)

xminATL5 = QbarATL5/ATLS5*np.exp(yminATL5)
xmaxATL5 = QbarATL5/ATLS5*np.exp(ymaxATL5)

#RHIC
SRHIC = 200.

yminRHIC = -2.2
ymaxRHIC =  2.2
QminRHIC =  4.8
QmaxRHIC =  8.2
QbarRHIC =  0.5*(QminRHIC+QmaxRHIC)

xminRHIC = QbarRHIC/SRHIC*np.exp(yminRHIC)
xmaxRHIC = QbarRHIC/SRHIC*np.exp(ymaxRHIC)

#E772
E772S = 2.*800*0.938
xfmin = 0.05
xfmax = 0.30

QminE772 = 4.
QmaxE772 = 9.
QbarE772 = 0.5*(QminE772+QmaxE772)

taumin = QbarE772*QbarE772/E772S
taumax = QbarE772*QbarE772/E772S

xminE772 = (-xfmin+np.sqrt(xfmin**2.+taumin))/2.
xmaxE772 = (-xfmax+np.sqrt(xfmax**2.+taumax))/2.

#E866
E866S = 2.*800*0.938
xfmin = 0.13
xfmax = 0.93

QminE866 = 4.
QmaxE866 = 8.
QbarE866 = 0.5*(QminE866+QmaxE866)

taumin = QbarE866*QbarE866/E866S
taumax = QbarE866*QbarE866/E866S

xminE866 = (-xfmin+np.sqrt(xfmin**2.+taumin))/2.
xmaxE866 = (-xfmax+np.sqrt(xfmax**2.+taumax))/2.

#HERMES
dftot = pd.DataFrame(columns=["Nu","Z","Q2","pt2","MULT-RATIO","STAT","SYS","fuu","fuua","DIS"])

for file in os.listdir("./HERMES_DIS/"):
    if "z.dat" in file:
        df = pd.read_csv("HERMES_DIS/"+file, delim_whitespace=True)
        df = df[df.Z<0.7]
        df = df[df.pt2<0.3]
        dftot = pd.concat([dftot,df])

nuhermx = max(dftot.Nu.tolist())
nuhermn = min(dftot.Nu.tolist())

Q2hermx = max(dftot.Q2.tolist())
Q2hermn = min(dftot.Q2.tolist())

Qhermx = np.sqrt(Q2hermx)
Qhermn = np.sqrt(Q2hermn)

xminherm = Q2hermn/2./0.938/nuhermx
xmaxherm = Q2hermx/2./0.938/nuhermn

ax = py.subplot(111)

CMSxxs = [xminCMS5    ,xmaxCMS5    ]
CMSQns = [QminCMS5**2.,QminCMS5**2.]
CMSQxs = [QmaxCMS5**2.,QmaxCMS5**2.]

ax.fill_between(CMSxxs,CMSQns,CMSQxs,alpha = 0.5,label = r"\rm CMS")

ATLxxs = [xminATL5    ,xmaxATL5    ]
ATLQns = [QminATL5**2.,QminATL5**2.]
ATLQxs = [QmaxATL5**2.,QmaxATL5**2.]

ax.fill_between(ATLxxs,ATLQns,ATLQxs,alpha = 0.5,label = r"\rm ATLAS",hatch = 'x')

RHICxxs = [xminRHIC    ,xmaxRHIC    ]
RHICQns = [QminRHIC**2.,QminRHIC**2.]
RHICQxs = [QmaxRHIC**2.,QmaxRHIC**2.]

ax.fill_between(RHICxxs,RHICQns,RHICQxs,alpha = 0.5,label = r"\rm RHIC")

E772xxs = [xminE772    ,xmaxE772    ]
E772Qns = [QminE772**2.,QminE772**2.]
E772Qxs = [QmaxE772**2.,QmaxE772**2.]

ax.fill_between(E772xxs,E772Qns,E772Qxs,alpha = 0.5,label = r"\rm E772")

E866xxs = [xminE866    ,xmaxE866    ]
E866Qns = [QminE866**2.,QminE866**2.]
E866Qxs = [QmaxE866**2.,QmaxE866**2.]

ax.fill_between(E866xxs,E866Qns,E866Qxs,alpha = 0.5,label = r"\rm E866",hatch = '+')

HERMxxs = [xminherm  ,xmaxherm  ]
HERMQns = [Qhermn**2.,Qhermn**2.]
HERMQxs = [Qhermx**2.,Qhermx**2.]

ax.fill_between(HERMxxs,HERMQns,HERMQxs,alpha = 0.5,label = r"\rm HERMES")

xvals = np.linspace(5e-5,1,int(1e5))
Q2vals = [0.95*140*140*x for x in xvals]

#ax.plot(xvals,Q2vals)
ax.fill_between(xvals,Q2vals,[1.]*len(xvals),alpha = 0.125,label = r"\rm EIC $\sqrt{S} = 140$")

ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(loc = 2,frameon = False,fontsize = 15)

ax.semilogy()
ax.semilogx()

ax.set_xlabel(r"\rm $x$",fontsize = 15)
ax.set_ylabel(r"\rm $Q^2$ (GeV$^2$)",fontsize = 15)

py.tight_layout()

py.savefig("xs.pdf")

