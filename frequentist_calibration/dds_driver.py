import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pickle
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/")
sys.path.append("/Users/williamrudisill/Documents/Conceptual_Runoff_Model/fortran/")
from ForcingModel import ForcingModel
from snowmodule17 import snowmodule17 as sm17
#from hymod import hymod as hm
from driver import driver as dr
import dds

mpl.use('Qt5Agg')



# load data ...
#start_date = "2017-04-15"
start_date = "2017-10-01"
end_date = "2019-09-28"
data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, daily_precip, qobs, day_len_hrs, daily_swe = frc(start_date, end_date)


# ---- ET DATA  ----
# read PET data
df = pd.read_csv("../data/Anna_etal_East_Flux_Tower.csv")
df['date'] = pd.to_datetime(df['date'])
df = df.set_index('date')
df = df.loc[start_date:end_date]


# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# Load the snow17 parameters
pkl_file = open('./parameter_files/snow17params.pkl', 'rb')
snow17params = pickle.load(pkl_file)
pkl_file.close()

### Do the elevation computing
nlayers = 3
nelev = int(elev.shape[0])
intrv = int(nelev/nlayers)
#intrvD = int(np.argwhere(elev>3700)[0])

# make an array of the fraction of each bin
areas = np.array([intrv,
                  intrv,
                  intrv])/nelev

# array of the mean elevations of each bin
heights = np.array([np.mean(elev[0:intrv]),
                    np.mean(elev[intrv:intrv*2]),
                    np.mean(elev[intrv*2:]),
                    ])

# process dates
dates = pd.date_range(start_date, end_date, freq='D')
jdays = dates.dayofyear.values * 1.0
ntimes=len(dates)

# create PET forcing
PET = np.zeros_like(jdays)
i=0
for l,t in zip(day_len_hrs, daily_temp):
   PET[i] = dr.pet(l,t)
   i+=1


# FUSE Model Parameters #
SNOWOP = 0                       #         ! OPTION   SNOW option
ETOP = 0                         #         ! OPTION   Percolation option
DRAINOP = 0                      #         ! OPTION   Snow option
SATOP = 0                        #         ! OPTION   Saturated Area Option
BASEOP = 0                       #         ! OPTION   Surface Runoff option
SMOP1 = 0                        #         ! OPTION   Soil Water Option -- Top
SMOP2 = 0                        #         ! OPTION   Soil Water Option -- Bottom
ntimes = ntimes                 #         ! FORCING   Number of model timesteps
PET = PET                        #         ! FORCING   Potential Evapotranspiratio
jday = jdays                     #         ! FORCING   Day of Year
tair = daily_temp                #         ! FORCING   Air Temperature = 1 length N
precip = daily_precip            #         ! FORCING   Precipitation = 1 length N
nlayers = 3                      #         ! PARAMETER   SNOW17

##### Run Snow17 Model ######
PET = PET                        #         ! FORCING   Potential Evapotranspiratio
jday = jdays                     #         ! FORCING   Day of Year
tair = daily_temp                #         ! FORCING   Air Temperature = 1 length N
precip = daily_precip            #         ! FORCING   Precipitation = 1 length N
nlayers = 3                      #         ! PARAMETER   SNOW17

# Snow model parameters
rvs = 1                          #         ! PARAMETER   SNOW17
opg_method = 1                   #         ! PARAMETER   SNOW17
dz  = heights - 3096.768                    #         ! PARAMETER   SNOW17
dt  = 24                          #         ! PARAMETER   SNOW17
opg = snow17params.get('opg')#.002                       #         ! PARAMETER   SNOW17
bias  = snow17params.get('bias') #.01                       #         ! PARAMETER   SNOW17
uadj  = snow17params.get('uadj')#.04                       #         ! PARAMETER   SNOW17
mbase = snow17params.get('mbase')#.3                       #         ! PARAMETER   SNOW17
mfmax = snow17params.get('mfmax')#4.0                      #         ! PARAMETER   SNOW17
mfmin =snow17params.get('mfmin')# 1.0                      #         ! PARAMETER   SNOW17
tipm  = snow17params.get('tipm')#.1                        #         ! PARAMETER   SNOW17
nmf   = snow17params.get('nmf')#.4                         #         ! PARAMETER   SNOW17
plwhc = snow17params.get('plwhc')#04                      #         ! PARAMETER   SNOW17
pxtemp  = snow17params.get('pxtemp')#2.                      #         ! PARAMETER   SNOW17
pxtemp1 = snow17params.get('pxtemp1')#-1.                    #         ! PARAMETER   SNOW17
pxtemp2 = snow17params.get('pxtemp2')#3.                     #         ! PARAMETER   SNOW17


# NOT IMPLEMENTED
sm1Fmax = 0.
alpha=0.
lam = 0.

# INITIAL CONDITIONS
sm1i  = 100.
sm2i  = 700.
sm1max = 200.
sm2max = 1000.
psi=0.0
nr =  2.0
kr =  1.0


# hydrologic model parameters ...
parameters= {'Parameter': ["sm1max", "sm2max", "ku", "c", "ks", "lowercasen", "beta"],
             'Init' :     [320.,     1285.,    80.,  1.0,  2.45, 1.0, 1.5],
             'Min' :      [100.,     500.,     1.,   .01, .001, .01,  .1],
             'Max' :      [800.,     1500.,    100., 3.0,  5.,  3.0, 4.0],
             "OnOff" : [True]*7}

parameters["BestValue"] = parameters["Init"]
parameters["ThisValue"] = parameters["Init"]
pdf = pd.DataFrame(parameters)
pdf = pdf.set_index("Parameter")


Nq = 4

def driver_calibrate(x0):
    qchan  = dr.model_driver(snowop   = SNOWOP,     #         ! OPTION   SNOW option
                        etop     = ETOP,       #         ! OPTION   Percolation option
                        drainop  = DRAINOP,    #         ! OPTION   Snow option
                        satop    = SATOP,      #         ! OPTION   Saturated Area Option
                        baseop   = BASEOP,     #         ! OPTION   Surface Runoff option
                        smop1    = SMOP1,      #         ! OPTION   Soil Water Option -- Top
                        smop2    = SMOP2,      #         ! OPTION   Soil Water Option -- Bottom
                        ntimes   = ntimes,     #         ! FORCING   Number of model timesteps
                        pet        = PET,        #         ! FORCING   Potential Evapotranspiratio
                        jday       = jday,       #         ! FORCING   Day of Year
                        tair       = tair,       #         ! FORCING   Air Temperature, length N
                        precip     = precip,     #         ! FORCING   Precipitation, length N
                        nlayers    = nlayers,    #         ! PARAMETER   SNOW17
                        rvs        = rvs,        #         ! PARAMETER   SNOW17
                        opg_method = opg_method, #         ! PARAMETER   SNOW17
                        dz         = dz,         #         ! PARAMETER   SNOW17
                        dt         = dt,         #         ! PARAMETER   SNOW17
                        opg        = opg,        #         ! PARAMETER   SNOW17
                        bias       = bias,       #         ! PARAMETER   SNOW17
                        uadj       = uadj,       #         ! PARAMETER   SNOW17
                        mbase      = mbase,      #         ! PARAMETER   SNOW17
                        mfmax      = mfmax,      #         ! PARAMETER   SNOW17
                        mfmin      = mfmin,      #         ! PARAMETER   SNOW17
                        tipm       = tipm,       #         ! PARAMETER   SNOW17
                        nmf        = nmf,        #         ! PARAMETER   SNOW17
                        plwhc      = plwhc,      #         ! PARAMETER   SNOW17
                        pxtemp     = pxtemp,     #         ! PARAMETER   SNOW17
                        pxtemp1    = pxtemp1,    #         ! PARAMETER   SNOW17
                        pxtemp2    = pxtemp2,    #           ! PARAMETER   SNOW17
                        sm1i       = sm1i,
                        sm2i       = sm2i,
                        sm1max     = x0[0],     #         ! PARAMETER   ?
                        sm2max     = x0[1],     #         ! PARAMETER   ?
                        ku         = x0[2],      #ku,           #         ! PARAMETER   PERCOLATION
                        c          = x0[3],      #c,            #         ! PARAMETER   PERCOLATION
                        sm1fmax    = sm1Fmax,    #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        psi        = psi,        #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        alpha      = alpha,      #         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        ks         = x0[4],#ks,           #         ! PARAMETER   BASEFLOW
                        lam        = lam,      #        ! PARAMETER   BASEFLOW  --- OPTIONAL
                        lowercasen = x0[5],#lowercasen,   #         ! PARAMETER   BASEFLOW  --- OPTIONAL
                        beta       = x0[6],#beta,         #         ! PARAMETER   SFROFF
                        nr         = nr,
                        kr         = kr)

    #qtot, qchan, qb, qsx, et, qin = output
    objective_function = dds.kge(qchan[1], qobs) # + dds.kge(qchan[4], df.values[:,0]))/2.
    return objective_function, qchan[1], qchan[4]   #dds.kge(et, df.values[:,0]))/2., q, et

iterations = 10000
dc = lambda x: driver_calibrate(x)[0]

# now calibrate. ...
parameters, kge_keeper, best_kge_score = dds.DDS(pdf, dc, iterations, r=.2)


# now plot the best
kge0, q0, et0= driver_calibrate(parameters['Init'])
kge, q, et = driver_calibrate(parameters['BestValue'])

# PLOT AND WRAP UP
plt.plot(q0, label='q0')
plt.plot(q, label='q')
plt.plot(qobs.values, label='qobs')
plt.legend()
plt.show()


print("DDS iterations: %s", iterations)
print("Starting kge: %s", kge0)
print("Ending kge score: %s", kge)


# Save
saveflag = input("save (y/n)")
if saveflag in ['y', 'yes']:
    with open("driver_params_calibrated.pkl", 'wb') as pkl:
        pickle.dump(parameters.BestValue, pkl)
    print("saved file to: driver_params_calibrated")
else:
    print("parameters not saved")



