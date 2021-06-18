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
from hymod import hymod as hm
import dds

mpl.use('Qt5Agg')

# read PET data
df = pd.read_csv("./data/Anna_etal_East_Flux_Tower.csv")
df['date'] = pd.to_datetime(df['date'])
df = df.set_index('date')

# load data ...
#start_date = "2017-04-15"
start_date = "2017-10-15"
end_date = "2019-09-28"
data_base = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/"

# load data ...
frc = ForcingModel()
frc.read()

# read in forcings for the desired time period
daily_temp, daily_precip, qobs, day_len_hrs, daily_swe = frc(start_date, end_date)

# read these in (they are fixed in time)
elev = np.load(data_base + "elev.npy")
dz = np.load(data_base + "dz.npy")

# Load the snow17 parameters
pkl_file = open('fortran/parameter_files/snow17params.pkl', 'rb')
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
   PET[i] = hm.pet(l,t)
   i+=1


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



parameters= {'Parameter': ["Kq", "Ks", "Alp", "Huz", "B"],
             'Init' : [0.371335, 0.132172, 0.240687, 13.022648, 0.468337 ],
             'Min' : [.001, .001, .001, 5., .01],
             'Max' : [.999, .999, .999, 1000., .5],
             "OnOff" : [True]*5}
parameters["BestValue"] = parameters["Init"]
parameters["ThisValue"] = parameters["Init"]
pdf = pd.DataFrame(parameters)
pdf = pdf.set_index("Parameter")


Nq = 4

def hymod_calibrate(x):
    q, et = hm.hymod_driver(ntimes     = ntimes,     #         ! FORCING   Number of model timesteps
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
                            pxtemp2    = pxtemp2,    #         ! PARAMETER   SNOW17
                            nq         = Nq,         #         ! PARAMETER   SNOW17
                            kq         = x[0],         #         ! PARAMETER   SNOW17
                            ks         = x[1],         #         ! PARAMETER   SNOW17
                            alp        = x[2],        #         ! PARAMETER   SNOW17
                            huz        = x[3],        #         ! PARAMETER   SNOW17
                            b          = x[4])          #           ! PARAMETER   SNOW17

    return dds.kge(q, qobs), q, et #+ dds.kge(et, df.values[:,0]))/2. , q, et

_hmc_ = lambda x: hymod_calibrate(x)[0]

# now calibrate. ...
parameters, kge_keeper, best_kge_score = dds.DDS(pdf, _hmc_, 100000)



# now plot the best
kge0, q0, et0 = hymod_calibrate(parameters['Init'])
kge, q, et = hymod_calibrate(parameters['BestValue'])

plt.plot(q, label='q0')
plt.plot(q0, label='q')
plt.plot(qobs.values, label='qobs')
plt.legend()
plt.show()




