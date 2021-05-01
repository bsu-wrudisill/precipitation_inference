# def part(x):
#   foo = dr.model_driver(snowop  = SNOWOP,     #         ! OPTION   SNOW option
#                         etop    = ETOP,       #         ! OPTION   Percolation option
#                         drainop = DRAINOP,    #         ! OPTION   Snow option
#                         satop  = SATOP,      #         ! OPTION   Saturated Area Option
#                         baseop = BASEOP,     #         ! OPTION   Surface Runoff option
#                         smop1  = SMOP1,      #         ! OPTION   Soil Water Option -- Top
#                         smop2  = SMOP2,      #         ! OPTION   Soil Water Option -- Bottom
#                         ntimes   = ntimes,     #         ! FORCING   Number of model timesteps
#                         pet        = PET,        #         ! FORCING   Potential Evapotranspiratio
#                         jday       = jday,       #         ! FORCING   Day of Year
#                         tair       = tair,       #         ! FORCING   Air Temperature, length N
#                         precip     = precip,     #         ! FORCING   Precipitation, length N
#                         nlayers    = nlayers,    #         ! PARAMETER   SNOW17
#                         rvs        = rvs,        #         ! PARAMETER   SNOW17
#                         opg_method = opg_method, #         ! PARAMETER   SNOW17
#                         dz         = dz,         #         ! PARAMETER   SNOW17
#                         dt         = dt,         #         ! PARAMETER   SNOW17
#                         opg        = opg,        #         ! PARAMETER   SNOW17
#                         bias       = bias,       #         ! PARAMETER   SNOW17
#                         uadj       = uadj,       #         ! PARAMETER   SNOW17
#                         mbase      = mbase,      #         ! PARAMETER   SNOW17
#                         mfmax      = mfmax,      #         ! PARAMETER   SNOW17
#                         mfmin      = mfmin,      #         ! PARAMETER   SNOW17
#                         tipm       = tipm,       #         ! PARAMETER   SNOW17
#                         nmf        = nmf,        #         ! PARAMETER   SNOW17
#                         plwhc      = plwhc,      #         ! PARAMETER   SNOW17
#                         pxtemp     = pxtemp,     #         ! PARAMETER   SNOW17
#                         pxtemp1    = pxtemp1,    #         ! PARAMETER   SNOW17
#                         pxtemp2    = pxtemp2,     #           ! PARAMETER   SNOW17
#                         sm1max     = x[0],#m1max,        #         ! PARAMETER   ?
#                         sm2max     = x[1],#m2max,       #         ! PARAMETER   ?
#                         ku         = x[2],#u,           #         ! PARAMETER   PERCOLATION
#                         c          = x[3],#,            #         ! PARAMETER   PERCOLATION
#                         sm1fmax    = x[4],#m1Fmax,      #         ! PARAMETER   PERCOLATION  --- OPTIONAL
#                         psi        = x[5],#si,          #         ! PARAMETER   PERCOLATION  --- OPTIONAL
#                         alpha      = x[6],#lpha,        #         ! PARAMETER   PERCOLATION  --- OPTIONAL
#                         ks         = x[7],#s,           #         ! PARAMETER   BASEFLOW
#                         lam        = x[8],#am,      #        ! PARAMETER   BASEFLOW  --- OPTIONAL
#                         lowercasen = x[9],#owercasen,   #         ! PARAMETER   BASEFLOW  --- OPTIONAL
#                         beta       = x[10])#beta)         #         ! PARAMETER   SFROFF

#   return foo

# def obj(q, fx, x):
#   mod=fx(x)
#   return np.sqrt(np.mean((mod - q)**2))



# x0 = [300.,
# 300.,
# .01 ,
# 1.  ,
# 240 ,
# .1  ,
# .1  ,
# .4  ,
# .1  ,
# 1.  ,
# 2.0]

# objective_func = lambda x: obj(q.values, part, x)

# # x0 =
# #result = scipy.optimize.minimize(objective_func, x0, method='Nelder-Mead')

