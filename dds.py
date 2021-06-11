import numpy as np
import pandas as pd
from Hymod import Hymod01
import matplotlib.pyplot as plt
import sys

libPath = "/Users/williamrudisill/Documents/Conceptual_Runoff_Model/"
sys.path.insert(0,libPath)
from ForcingModel import *




def kge(mod, obs):
    #
    mobs = np.mean(obs)
    sobs = np.std(obs)

    # mean ratio
    b = np.mean(mod) / mobs
    # std
    a = np.std(mod) / sobs
    # corr coeff
    r = np.corrcoef(mod, obs)[0, 1]  # corrcoef returns the correlation matrix...
    # the diagonals are 1, the off-diags are the 'r'
    # value that we want
    kgeval = 1-np.sqrt((r - 1.)**2 + (a - 1.)**2 + (b - 1)**2)
    return kgeval


#def main(Nq, Kq, Ks, Alp, Huz, B, Data_name):
# read in observed rainfall-runoff data for one year
Data = pd.read_csv("forcing_data.csv")

Nq = 7

# Initialize states
InState = {'Xq': np.zeros(Nq)+1.,
           'Xs': 0,
           'XHuz': 0}


parameters= {'Parameter': ["Kq", "Ks", "Alp", "Huz", "B"],
             'Init' : [0.371335, 0.132172, 0.240687, 13.022648, 0.468337 ],
             'Min' : [.001, .001, .001, 5., .01],
             'Max' : [.999, .999, .999, 1000., .5],
             "OnOff" : [True]*5}

# add these keys...
parameters["BestValue"] = parameters["Init"]
parameters["ThisValue"] = parameters["Init"]


max_iteration = 10000
r = .2

df = pd.DataFrame(parameters)
df = df.set_index("Parameter")

best_kge_score = -100.
kge_keeper = np.zeros(max_iteration)

for iteration in range(max_iteration):

    # probability of selecting a parameter to calibrate
    prob = 1 - np.log(iteration + 1) / np.log(max_iteration + 1)
    parameter_names = df.index

    # loop through the parameter list
    for param in parameter_names:
        sel = np.random.choice(2, p=[1 - prob, prob])
        if sel == 1:
            df.at[param, 'OnOff'] = True
        else:
            df.at[param, 'OnOff'] = False

    # pick one at random if all are set to false
    if (df.OnOff == False).all():
       backon = np.random.randint(0,4)

       # NOTE -- the following does not work in pandas...
       # df.iloc[backon].at["OnOff"] = True ## this will not assign anything... so dumb
       df.iat[backon, 3] = True # three is the "onoff" column...

    # Get the parameters that are ON
    for turned_on in df.groupby("OnOff").groups[True]:
        best_value = df.loc[turned_on].at["BestValue"]

        # get parameter info... doesnt change
        xj_min = df.loc[turned_on].at["Min"]
        xj_max = df.loc[turned_on].at["Max"]
        xj_init = df.loc[turned_on].at["Init"]

        # define the std of the perturn dist
        sigj = r * (xj_max - xj_min)

        # normally distributed
        x_update = sigj * np.random.randn(1) + xj_init

        x_new = df.at[turned_on, "BestValue"] + x_update

        # If new factor is l.t min, reflect to middle
        if x_new < xj_min:
            x_new = xj_min + (xj_min - x_new)

        # If new factor is g.t max, reflect to middle
        elif x_new > xj_max:
            x_new = xj_max - (x_new - xj_max)


        # update the dataframe
        df.at[turned_on,"ThisValue"] = x_new

    # Only do this if there are any parameters turned off...
    if (df.OnOff == False).any():
        # Make sure the parameters that are turned OFF are set to the 'BestValue'...
        for turned_off in df.groupby("OnOff").groups[False]:
            df.at[turned_off,"ThisValue"] = df.loc[turned_off, "BestValue"]

    # NOW RUN THE MODEL
    x0 = df.ThisValue.values
    Model = Hymod01(Data, x0, InState)

    # Evaluate the Model
    kge_score = kge(Model['Q'],Data['Qobs'])

    # check performance
    if kge_score > best_kge_score:
        df["BestValue"] = df["ThisValue"]
        best_kge_score = kge_score
        # #print(iteration, 'better performance')
        # best_kge_score = kge_score
    # else:
        # do nothing

    kge_keeper[iteration] = kge_score



