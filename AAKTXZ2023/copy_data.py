import matplotlib
matplotlib.use('Agg')
import os
import pandas as pd
import numpy as np
import pylab as py
import warnings
from  matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib import font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
#from scipy.interpolate import spline,interp1d
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
warnings.filterwarnings('ignore')

src = './fit'

for j in range(101,301): # 200 --> 100 replicas for now
    if len(str(j)) == 1:
        numb = '00'+str(j)
    elif len(str(j)) == 2:
        numb = '0'+str(j)
    else:
        numb = str(j)
    dest = './replica_'+numb
    for filename in os.listdir(src+'/expdata/DRELLYAN/E288/'):
        path2src = src+'/expdata/DRELLYAN/E288/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename, delim_whitespace=True)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['error'][i])
            gas_rand = np.random.normal(scale=sigma)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['value'][i] = df['value'][i]+gas_rand
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
        #    os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/E288/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/E906/'):
        path2src = src+'/expdata/DRELLYAN/E906/'
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['Err'][i])
            gas_rand = np.random.normal(scale=sigma)
            df['R'][i] = df['R'][i]+gas_rand
        rep0dir = dest+'/expdata'
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        path2srcrep = path2srcrep1 + 'E906/'
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/ATLAS5/'):
        path2src = src+'/expdata/DRELLYAN/ATLAS5/'
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(np.sqrt((df['stat+'][i])**2.+(df['sys+'][i])**2.))
            gas_rand = np.random.normal(scale=sigma)
            df['1/pTdsigma/dpT[nb/GeV^2]'][i] = df['1/pTdsigma/dpT[nb/GeV^2]'][i]+gas_rand
        rep0dir = dest+'/expdata'
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        path2srcrep = path2srcrep1 + '/ATLAS5/'
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/ATLAS3/'):
        path2src = src+'/expdata/DRELLYAN/ATLAS3/'
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(np.sqrt((df['stat+'][i])**2.+(df['sys+'][i])**2.))
            gas_rand = np.random.normal(scale=sigma)
            df['(10**7/Nev)*(1/PT)*DN/DPT[GEV**-2]'][i] = df['(10**7/Nev)*(1/PT)*DN/DPT[GEV**-2]'][i]+gas_rand
        rep0dir = dest+'/expdata'
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        path2srcrep = path2srcrep1 + '/ATLAS3/'
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/CMS5/'):
        path2src = src+'/expdata/DRELLYAN/CMS5/'
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(np.sqrt((df['stat+'][i])**2.+(df['sys+'][i])**2.))
            gas_rand = np.random.normal(scale=sigma)
            df['DSIG/DPT[PB/GeV]'][i] = df['DSIG/DPT[PB/GeV]'][i]+gas_rand
        rep0dir = dest+'/expdata'
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        path2srcrep = path2srcrep1 + '/CMS5/'
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/CMS8/'):
        path2src = src+'/expdata/DRELLYAN/CMS8/'
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
#            if j == 0:
#                sigma = 0.
#            else:
#                sigma = abs(np.sqrt((df['stat+'][i])**2.+(df['sys+'][i])**2.))
            sigma = 0.
            gas_rand = np.random.normal(scale=sigma)
            df['sig'][i] = df['sig'][i]+gas_rand
        rep0dir = dest+'/expdata'
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        path2srcrep = path2srcrep1 + '/CMS8/'
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/E605/'):
        path2src = src+'/expdata/DRELLYAN/E605/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['error'][i])
            gas_rand = np.random.normal(scale=sigma)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['ds'][i] = df['ds'][i]+gas_rand
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
            #os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/E605/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/E605/'):
        path2src = src+'/expdata/DRELLYAN/E605/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['error'][i])
            gas_rand = np.random.normal(scale=sigma)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['ds'][i] = df['ds'][i]+gas_rand
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
            #os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/E605/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/E772/'):
        path2src = src+'/expdata/DRELLYAN/E772/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['error'][i])
            gas_rand = np.random.normal(scale=sigma)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['R'][i] = df['R'][i]+gas_rand
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
            #os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/E772/'
            #if not os.path.exists(path2srcrep):
            #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/E866/'):
        path2src = src+'/expdata/DRELLYAN/E866/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma1 = 0.
                sigma2 = 0.
            else:
                sigma1 = abs(df['err1'][i])
                sigma2 = abs(df['err2'][i])
            gas_rand1 = np.random.normal(scale=sigma1)
            gas_rand2 = np.random.normal(scale=sigma2)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['Fe/Be'][i] = df['Fe/Be'][i]+gas_rand1
            df['W/Be'][i] = df['W/Be'][i]+gas_rand2
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
            #os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/E866/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/RHIC/'):
        path2src = src+'/expdata/DRELLYAN/RHIC/'
        #skips = pd.read_csv(path2src+filename, nrows=0).columns.values[0]
        #nskips = skips.split('#comment lines: ')[1]
        df = pd.read_csv(path2src+filename)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = abs(df['error'][i])
            gas_rand = np.random.normal(scale=sigma)
            df['CX'][i] = df['CX'][i]+gas_rand
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep1 = rep0dir+'/DRELLYAN/'
        #if not os.path.exists(path2srcrep1):
            #os.mkdir(path2srcrep1)
        path2srcrep = path2srcrep1 + '/RHIC/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        df = df.round(3)
        df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/DRELLYAN/RHIC2/'):
        if(filename != "RHIC_pp.dat"):
            path2src = src+'/expdata/DRELLYAN/RHIC2/'
            df = pd.read_csv(path2src+filename)
            for i in range(len(df)):
                if j == 0:
                    sigma = 0.
                else:
                    sigma = abs(df['Error'][i])
                gas_rand = np.random.normal(scale=sigma)
                df['RpAU'][i] = df['RpAU'][i]+gas_rand
            rep0dir = dest+'/expdata'
            path2srcrep1 = rep0dir+'/DRELLYAN/'
            path2srcrep = path2srcrep1 + '/RHIC2/'
            df = df.round(3)
            df.to_csv(path2srcrep+filename, sep='\t', index = False)

    for filename in os.listdir(src+'/expdata/SIDIS/HERMES_DIS/'):
        path2src = src+'/expdata/SIDIS/HERMES_DIS/'
        if filename.startswith('.'):  # Skip hidden files
            continue
        df = pd.read_csv(path2src+filename, delim_whitespace=True)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma = np.sqrt(df['STAT'][i]**2.+df['SYS'][i]**2.)
            gas_rand1 = np.random.normal(scale=sigma)
            #gas_rand = 2*sigma
            #while abs(gas_rand)>sigma:
            #    gas_rand = np.random.normal(scale=sigma)
            #initial = df['Sivers'][i]
            df['MULT-RATIO'][i] = df['MULT-RATIO'][i]+gas_rand1
            #final = df['value'][i]
            #if abs(final-initial) > sigma:
            #    print path2src+filename+numb+', something is wrong'
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep = rep0dir+'/SIDIS/HERMES_DIS/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        #df = df.round(3)
        df.to_csv(path2srcrep+filename, index = False)



    for filename in os.listdir(src+'/expdata/SIDIS/JLAB_DIS/'):
        path2src = src+'/expdata/SIDIS/JLAB_DIS/'
        df = pd.read_csv(path2src+filename, delim_whitespace=True)
        for i in range(len(df)):
            if j == 0:
                sigma_C = 0.
                sigma_Fe = 0.
                sigma_Pb = 0.
            else:
                sigma_C  = abs(df['Cerr'][i])
                sigma_Fe = abs(df['FEerr'][i])
                sigma_Pb = abs(df['PBerr'][i])
            gas_rand_C   = np.random.normal(scale=sigma_C)
            gas_rand_Fe  = np.random.normal(scale=sigma_Fe)
            gas_rand_Pb  = np.random.normal(scale=sigma_Pb)
            df['C'][i]   = df['C'][i]+ gas_rand_C
            df['Fe'][i]  = df['Fe'][i]+ gas_rand_Fe
            df['Pb'][i]  = df['Pb'][i]+ gas_rand_Pb
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep = rep0dir+'/SIDIS/JLAB_DIS/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        #df = df.round(3)
        df.to_csv(path2srcrep+filename, index = False)
    for filename in os.listdir(src+'/expdata/SIDIS/JLAB12_DIS/'):
        path2src = src+'/expdata/SIDIS/JLAB12_DIS/'
        df = pd.read_csv(path2src+filename, delim_whitespace=True)
        for i in range(len(df)):
            if j == 0:
                sigma = 0.
            else:
                sigma  = abs(df['ERR'][i])
            gas_rand   = np.random.normal(scale=sigma)
            df['MULT-RATIO'][i]   = df['MULT-RATIO'][i]+ gas_rand
        rep0dir = dest+'/expdata'
        #if not os.path.exists(rep0dir):
        #    os.mkdir(rep0dir)
        path2srcrep = rep0dir+'/SIDIS/JLAB12_DIS/'
        #if not os.path.exists(path2srcrep):
        #    os.mkdir(path2srcrep)
        #df = df.round(3)
        df.to_csv(path2srcrep+filename, index = False)

    for filename in os.listdir(os.path.join(src, 'expdata', 'JLAB2022')):
        if filename in ['pi+.csv', 'pi-.csv']:
            path2src = os.path.join(src, 'expdata', 'JLAB2022')
            df = pd.read_csv(os.path.join(path2src, filename), delimiter=',')
            for i in range(len(df)):
                if j == 0:
                    sigma_C = sigma_Fe = sigma_Pb = 0.
                else:
                    sigma_C = ((df['Cstat'][i])**2 + (df['Csys'][i])**2)**0.5
                    sigma_Fe = ((df['Festat'][i])**2 + (df['Fesys'][i])**2)**0.5
                    sigma_Pb = ((df['Pbstat'][i])**2 + (df['Pbsys'][i])**2)**0.5

                gas_rand_C  = np.random.normal(scale=sigma_C)
                gas_rand_Fe = np.random.normal(scale=sigma_Fe)
                gas_rand_Pb = np.random.normal(scale=sigma_Pb)

                df['C'][i] += gas_rand_C
                df['Fe'][i] += gas_rand_Fe
                df['Pb'][i] += gas_rand_Pb

            rep0dir = os.path.join(dest, 'expdata')
            path2srcrep = os.path.join(rep0dir, 'JLAB2022')
            df.to_csv(os.path.join(path2srcrep, filename), index = False)  
