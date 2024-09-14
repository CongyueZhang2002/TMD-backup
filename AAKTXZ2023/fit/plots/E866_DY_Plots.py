#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex = True)


FNAL_E866_FeBe_x2bin_df = pd.read_csv('../plot_data/E866_800_FeBe.dat', sep ='\s+')
FNAL_E866_WBe_x2bin_df = pd.read_csv('../plot_data/E866_800_WBe.dat', sep ='\s+')

FNAL_E866_FeBe_df_off = pd.read_csv('../plot_data/E866_800_FeBe_off.dat', sep ='\s+')
FNAL_E866_WBe_df_off = pd.read_csv('../plot_data/E866_800_WBe_off.dat', sep ='\s+')

##Check Data
#FNAL_E866_FeBe_x2bin_df.head(40)

FNAL_E866_FeBe_Qbin_df = pd.read_csv('../plot_data/E866_800_Qbin_FeBe.dat', sep ='\s+')
FNAL_E866_WBe_Qbin_df = pd.read_csv('../plot_data/E866_800_Qbin_WBe.dat', sep ='\s+')

Qbar = 6.2
class dydata:
    def __init__(self, file, i1, i2):
        self.xdata = file['pt'][i1:i2]
        self.ydata = file['DY-RATIO'][i1:i2]
        self.error = file['error'][i1:i2]
        self.theory = file['R_dy'][i1:i2]
        self.num    = len(self.xdata)

Qbar = 6.2
class dydata3:
    def __init__(self, xdata, ydata, error, theory):
        self.xdata = xdata[(xdata/Qbar < 0.3)]
        self.ydata = ydata[(xdata/Qbar < 0.3)]
        self.error = error[(xdata/Qbar < 0.3)]
        self.theory = theory[(xdata/Qbar < 0.3)]
        self.num    = len(self.xdata)

class dydata2:
    def __init__(self, file, i1, i2):
        self.xdata = file['pt'][i1:i2]
        self.ydata = file['DY-RATIO'][i1:i2]
        self.error = file['error'][i1:i2]

FeBe_Data_Qbin1 = dydata(FNAL_E866_FeBe_Qbin_df, 0,3)
FeBe_Data_Qbin2 = dydata(FNAL_E866_FeBe_Qbin_df, 3,6)
FeBe_Data_Qbin3 = dydata(FNAL_E866_FeBe_Qbin_df, 6,10)
FeBe_Data_Qbin4 = dydata(FNAL_E866_FeBe_Qbin_df, 10,14)

WBe_Data_Qbin1 = dydata(FNAL_E866_WBe_Qbin_df, 0,3)
WBe_Data_Qbin2 = dydata(FNAL_E866_WBe_Qbin_df, 3,6)
WBe_Data_Qbin3 = dydata(FNAL_E866_WBe_Qbin_df, 6,10)
WBe_Data_Qbin4 = dydata(FNAL_E866_WBe_Qbin_df, 10,14)


FeBe_Data_x2bin1 = dydata(FNAL_E866_FeBe_x2bin_df, 0,4)
FeBe_Data_x2bin2 = dydata(FNAL_E866_FeBe_x2bin_df, 4,8)
FeBe_Data_x2bin3 = dydata(FNAL_E866_FeBe_x2bin_df, 8,12)
FeBe_Data_x2bin4 = dydata(FNAL_E866_FeBe_x2bin_df, 12,16)
FeBe_Data_x2bin5 = dydata(FNAL_E866_FeBe_x2bin_df, 16,20)


WBe_Data_x2bin1 = dydata(FNAL_E866_WBe_x2bin_df, 0,4)
WBe_Data_x2bin2 = dydata(FNAL_E866_WBe_x2bin_df, 4,8)
WBe_Data_x2bin3 = dydata(FNAL_E866_WBe_x2bin_df, 8,12)
WBe_Data_x2bin4 = dydata(FNAL_E866_WBe_x2bin_df, 12,16)
WBe_Data_x2bin5 = dydata(FNAL_E866_WBe_x2bin_df, 16,20)

FeBe_Data_off_x2bin1 = dydata3(FNAL_E866_FeBe_df_off['pt'][0:6],
                          FNAL_E866_FeBe_df_off['DY-RATIO'][0:6],
                          FNAL_E866_FeBe_df_off['error'][0:6] ,
                          FNAL_E866_FeBe_df_off['R_dy'][0:6] )

FeBe_Data_off_x2bin2 = dydata3(FNAL_E866_FeBe_df_off['pt'][7:15],
                          FNAL_E866_FeBe_df_off['DY-RATIO'][7:15],
                          FNAL_E866_FeBe_df_off['error'][7:15] ,
                          FNAL_E866_FeBe_df_off['R_dy'][7:15] )

FeBe_Data_off_x2bin3 = dydata3(FNAL_E866_FeBe_df_off['pt'][16:25],
                          FNAL_E866_FeBe_df_off['DY-RATIO'][16:25],
                          FNAL_E866_FeBe_df_off['error'][16:25] ,
                          FNAL_E866_FeBe_df_off['R_dy'][16:25] )

FeBe_Data_off_x2bin4 = dydata3(FNAL_E866_FeBe_df_off['pt'][26:33],
                          FNAL_E866_FeBe_df_off['DY-RATIO'][26:33],
                          FNAL_E866_FeBe_df_off['error'][26:33],
                          FNAL_E866_FeBe_df_off['R_dy'][26:33] )


FeBe_Data_off_x2bin5 = dydata3(FNAL_E866_FeBe_df_off['pt'][34:39],
                          FNAL_E866_FeBe_df_off['DY-RATIO'][34:39],
                          FNAL_E866_FeBe_df_off['error'][34:39],
                          FNAL_E866_FeBe_df_off['R_dy'][34:39] )


WBe_Data_off_x2bin1 = dydata3(FNAL_E866_WBe_df_off['pt'][0:6],
                          FNAL_E866_WBe_df_off['DY-RATIO'][0:6],
                          FNAL_E866_WBe_df_off['error'][0:6],
                          FNAL_E866_WBe_df_off['R_dy'][0:6] )

WBe_Data_off_x2bin2 = dydata3(FNAL_E866_WBe_df_off['pt'][7:15],
                          FNAL_E866_WBe_df_off['DY-RATIO'][7:15],
                          FNAL_E866_WBe_df_off['error'][7:15] ,
                          FNAL_E866_WBe_df_off['R_dy'][7:15] )

WBe_Data_off_x2bin3 = dydata3(FNAL_E866_WBe_df_off['pt'][16:25],
                          FNAL_E866_WBe_df_off['DY-RATIO'][16:25],
                          FNAL_E866_WBe_df_off['error'][16:25] ,
                          FNAL_E866_WBe_df_off['R_dy'][16:25] )

WBe_Data_off_x2bin4 = dydata3(FNAL_E866_WBe_df_off['pt'][26:33],
                          FNAL_E866_WBe_df_off['DY-RATIO'][26:33],
                          FNAL_E866_WBe_df_off['error'][26:33],
                          FNAL_E866_WBe_df_off['R_dy'][26:33] )

WBe_Data_off_x2bin5 = dydata3(FNAL_E866_WBe_df_off['pt'][34:39],
                          FNAL_E866_WBe_df_off['DY-RATIO'][34:39],
                          FNAL_E866_WBe_df_off['error'][34:39],
                          FNAL_E866_WBe_df_off['R_dy'][34:39] )

def dyplot(axes, data, clr):
    axes.errorbar(data.xdata, data.ydata, data.error, fmt = "o",color = clr, ecolor = clr, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False)
    axes.plot(data.xdata, data.theory, color = clr, label = r'\rm Modified nDS')
    axes.tick_params(axis = 'both', direction ='in',labelsize = 25)


def dyplot2(axes, data, clr):
    axes.plot(data.xdata, data.theory, color = clr, label = r'\rm nDS', linestyle = 'dashed')
    axes.tick_params(axis = 'both', direction ='in',labelsize = 25)

def dyplot3(axes, data, clr):
    axes.errorbar(data.xdata, data.ydata, data.error, fmt = "o",color = clr, ecolor = clr, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False)
    axes.tick_params(axis = 'both', direction ='in',labelsize = 25)

def dyplotlabel(axes, data, clr, lbl):
    axes.errorbar(data.xdata, data.ydata, data.error, fmt = "o",color = clr, ecolor = clr, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False,label = lbl)
    axes.plot(data.xdata, data.theory, color = 'black', linestyle = 'dashed')
    axes.tick_params(axis = 'both', direction ='in',labelsize = 25)

def dyplotlabel2(axes, data, clr, lbl):
    axes.errorbar(data.xdata, data.ydata, data.error, fmt = "o",color = clr, ecolor = clr, uplims=False, lolims=False, capsize = 5, capthick =0.8, barsabove = False,label = lbl)
    axes.tick_params(axis = 'both', direction ='in',labelsize = 25)


fig2, axs2 = plt.subplots(2, 5, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})

#Fe/Be
dyplot(axs2[0,0],FeBe_Data_x2bin1,'blue')
dyplot(axs2[0,1],FeBe_Data_x2bin2,'blue')
dyplot(axs2[0,2],FeBe_Data_x2bin3,'blue')
dyplot(axs2[0,3],FeBe_Data_x2bin4,'blue')
dyplot(axs2[0,4],FeBe_Data_x2bin5,'blue')
dyplotlabel(axs2[0,4],FeBe_Data_x2bin5,'blue', r"\rm E866 (800 GeV)")

#W/Be
dyplot(axs2[1,0],WBe_Data_x2bin1,'blue')
dyplot(axs2[1,1],WBe_Data_x2bin2,'blue')
dyplot(axs2[1,2],WBe_Data_x2bin3,'blue')
dyplot(axs2[1,3],WBe_Data_x2bin4,'blue')
dyplot(axs2[1,4],WBe_Data_x2bin5,'blue')

#Proton Width
#Fe/Be
#dyplot2(axs2[0,0],FeBe_Data_off_x2bin1,'black')
#dyplot2(axs2[0,1],FeBe_Data_off_x2bin2,'black')
#dyplot2(axs2[0,2],FeBe_Data_off_x2bin3,'black')
#dyplot2(axs2[0,3],FeBe_Data_off_x2bin4,'black')
#dyplot2(axs2[0,4],FeBe_Data_off_x2bin5,'black')

#W/Be
#dyplot2(axs2[1,0],WBe_Data_off_x2bin1,'black')
#dyplot2(axs2[1,1],WBe_Data_off_x2bin2,'black')
#dyplot2(axs2[1,2],WBe_Data_off_x2bin3,'black')
#dyplot2(axs2[1,3],WBe_Data_off_x2bin4,'black')
#dyplot2(axs2[1,4],WBe_Data_off_x2bin5,'black')

# Axes Limits
# Y axis
axs2[0,0].set_ylim(0.5,1.5)
axs2[1,0].set_ylim(0.5,1.5)
# X axis
axs2[1,0].set_xlim(0,2.2)
axs2[1,1].set_xlim(0.1,2)
axs2[1,2].set_xlim(0.1,2)
axs2[1,3].set_xlim(0.1,2.2)
axs2[1,4].set_xlim(0.1,2.2)

# Set Labels
axs2[0,0].set_ylabel(r"\rm $\sigma^{Fe}/\sigma^{Be}$", fontsize = 25,rotation = 0,labelpad = 55)
axs2[1,0].set_ylabel(r"\rm $\sigma^{W}/\sigma^{Be}$", fontsize = 25,rotation = 0,labelpad = 55)
axs2[1,2].set_xlabel(r"\rm $p_T $ (GeV)", fontsize = 25, labelpad = 20)

# Set Figure Size
fig2.set_size_inches(21,15)
props = dict(boxstyle='round', facecolor='white', alpha=0)
props2 = dict(boxstyle='round', facecolor='gray', alpha=0.4)

axs2[0,0].text(0.2, 1.1,  r'\rm $0.02 < x_2 < 0.04$', transform=axs2[0,0].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs2[0,1].text(0.2, 1.1,  r'\rm $0.04 < x_2 < 0.06$', transform=axs2[0,1].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs2[0,2].text(0.2, 1.1,  r'\rm $0.06 < x_2 < 0.08$', transform=axs2[0,2].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs2[0,3].text(0.2, 1.1,  r'\rm $0.08 < x_2 < 0.10$', transform=axs2[0,3].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs2[0,4].text(0.2, 1.1,  r'\rm $0.10 < x_2 < 0.12$', transform=axs2[0,4].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

#axs2[0,4].legend(frameon = False,fontsize = 25,loc='upper center', bbox_to_anchor=(1.6, 0.4))

#save figs
fig2.savefig('Drell-Yan_E866_x2bins.pdf', bbox_inches = "tight")

#number of data pts lost
total_data_pts = 80

#REMAINING POINTS
remaining_pts =    ( FeBe_Data_x2bin1.num + FeBe_Data_x2bin2.num + FeBe_Data_x2bin3.num  +
                     FeBe_Data_x2bin4.num + FeBe_Data_x2bin5.num + WBe_Data_x2bin1.num   +
                      WBe_Data_x2bin2.num +  WBe_Data_x2bin3.num + WBe_Data_x2bin4.num   +
                      WBe_Data_x2bin4.num +  WBe_Data_x2bin5.num )

cut_pts = total_data_pts - remaining_pts

print('\ntotal points without cuts is:', total_data_pts)
print('total points remaining is:   ', remaining_pts, '(', round(100*remaining_pts/total_data_pts), '%)')
print('total points cut is:         ', cut_pts, '(', round(100*cut_pts/total_data_pts), '%)')


fig3, axs3 = plt.subplots(2, 4, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})

#Fe/Be
dyplot(axs3[0,0],FeBe_Data_Qbin1,'blue')
dyplot(axs3[0,1],FeBe_Data_Qbin2,'blue')
dyplot(axs3[0,2],FeBe_Data_Qbin3,'blue')
dyplot(axs3[0,3],FeBe_Data_Qbin4,'blue')
dyplotlabel2(axs3[0,3],FeBe_Data_Qbin4,'blue', r"\rm E866 (800 GeV)")


#W/Be
dyplot(axs3[1,0],WBe_Data_Qbin1,'blue')
dyplot(axs3[1,1],WBe_Data_Qbin2,'blue')
dyplot(axs3[1,2],WBe_Data_Qbin3,'blue')
dyplot(axs3[1,3],WBe_Data_Qbin4,'blue')

# Axes Limits
# Y axis
axs3[0,0].set_ylim(0.5,1.5)
axs3[1,0].set_ylim(0.5,1.5)
# X axis
axs3[1,0].set_xlim(0,2.2)
axs3[1,1].set_xlim(0.1,2)
axs3[1,2].set_xlim(0.1,2)
axs3[1,3].set_xlim(0.1,2.2)

# Set Labels
axs3[0,0].set_ylabel(r"\rm $\sigma^{Fe}/\sigma^{Be}$", fontsize = 25,rotation = 0,labelpad = 55)
axs3[1,0].set_ylabel(r"\rm $\sigma^{W}/\sigma^{Be}$", fontsize = 25,rotation = 0,labelpad = 55)
axs3[1,2].set_xlabel(r"\rm $p_T $ (GeV)", fontsize = 25, labelpad = 20)

# Set Figure Size
fig3.set_size_inches(21,15)
props = dict(boxstyle='round', facecolor='white', alpha=0)
props2 = dict(boxstyle='round', facecolor='gray', alpha=0.4)


axs3[0,0].text(0.2, 1.1,  r'\rm $4 < Q < 5$ (GeV)', transform=axs3[0,0].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs3[0,1].text(0.2, 1.1,  r'\rm $5 < Q < 6$ (GeV)', transform=axs3[0,1].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs3[0,2].text(0.2, 1.1,  r'\rm $6 < Q < 7$ (GeV)', transform=axs3[0,2].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

axs3[0,3].text(0.2, 1.1,  r'\rm $7 < Q < 8$ (GeV)', transform=axs3[0,3].transAxes, fontsize=25, color = 'black',
    verticalalignment='top', bbox=props2)

#axs3[0,3].legend(frameon = False,fontsize = 25,loc='upper center', bbox_to_anchor=(1.5, 0.4))


#save figs
fig3.savefig('Drell-Yan_E866_Qbins.pdf', bbox_inches = "tight")
