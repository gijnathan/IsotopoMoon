import pandas as pd
import numpy as np
import math
import datetime
import isotope_evolution as ie
import time

##################################################################
#                        INPUTS.                                 #
##################################################################

t_step_Myr = 1.e-1 # time step in Myr
end_time_Myr = 200 #4570 # end time in Myr
nr_step = 1.
nr_tol = 1.e5

### RESERVOIRS ###
# deep mantle
DM_mass_f = "N" # "N" if 0 mass. mantle mass S fraction of min and max of total mantle
DM_mass_f_sil = "N" # "N" if 0 mass. mantle mass silicate fraction of min and max of total mantle
DM_ST = "N"
DM_33d = ""
DM_34d = ""
DM_36d = ""
# mantle
M_mass_f = 0.1257 # mantle mass S fraction of min and max [0.1257]
M_mass_f_sil = M_mass_f # mantle mass silicate fraction of min and max of total mantle
M_ST = "N"
M_33d = 0.
M_34d = 0.
M_36d = 0.
# frost (crust)
F_mass_f = 0.
F_ST = "N"
F_33d = ""
F_34d = ""
F_36d = "" 
# silicate-sulfate (crust)
SS_mass_f = 0.
SS_ST = "N"
SS_33d = ""
SS_34d = ""
SS_36d = "" 
# space
S_ST = "N"
S_33d = ""
S_34d = ""
S_36d = ""

rate_f = {"pd":1., # [1. all] "N" = 0.
         "pi":1.,
         "ei":1.,
         "ed":1.,
         "ac":1.,
         "rc":1.,
         "ec":1.}

# options
oscillate = "N" # does the system oscillate or not?
resurf_cm_yr = 1. # resrufacing rate in cm/yr (today rate = 1.0)
sil_mag_S = 0.001 # mass fraction S in mantle melt
thick_C = 40000. # crustal thickness m - min = 30000, max = 50000
f_S2 = 0.2 # fraction of S2 for SO2 due to homogenous gas equilibria [0.2]
f_pl2mo = 0.8 # fraction of mantle melting that goes to plutons [0.8]
f_sq = 0.25 # fraction of sulfate sequestration from gas [0.25]
f_deep = 0. # fraction of mantle melting rate that is added from the deep mantle if present [0.]
f_remobilised = 1. # fraction of f_crust_returned that is remobilised [1.]
f_pu = 0.5 # fraction of material that could be lost to space that is [0.5]


##################################################################
#                        RUN MODEL                               #
##################################################################

# start time
start_time = time.time()
print("Start time:", datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

# Call the model

try:
    ie.iso_evo(
        t_step_Myr, end_time_Myr, nr_step, nr_tol,
        DM_mass_f, DM_mass_f_sil, DM_ST, DM_33d, DM_34d, DM_36d,
        M_mass_f, M_mass_f_sil, M_ST, M_33d, M_34d, M_36d,
        F_mass_f, F_ST, F_33d, F_34d, F_36d,
        SS_mass_f, SS_ST, SS_33d, SS_34d, SS_36d,
        S_ST, S_33d, S_34d, S_36d,
        rate_f, oscillate, resurf_cm_yr, sil_mag_S, thick_C,
        f_S2, f_pl2mo, f_sq, f_deep, f_remobilised, f_pu
    )

except KeyboardInterrupt:
    print("\n[Interrupted by user]")

finally:
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.2f} seconds")

# ie.iso_evo(t_step_Myr,end_time_Myr,nr_step,nr_tol,DM_mass_f,DM_mass_f_sil,DM_ST,DM_33d,DM_34d,DM_36d,M_mass_f,M_mass_f_sil,M_ST,M_33d,M_34d,M_36d,F_mass_f,F_ST,F_33d,F_34d,F_36d,SS_mass_f,SS_ST,SS_33d,SS_34d,SS_36d,S_ST,S_33d,S_34d,S_36d,rate_f,oscillate,resurf_cm_yr,sil_mag_S,thick_C,f_S2,f_pl2mo,f_sq,f_deep,f_remobilised,f_pu)
#ie.iso_evo_simple(t_step_Myr,end_time_Myr,nr_step,nr_tol,DM_mass_f,DM_mass_f_sil,DM_ST,DM_33d,DM_34d,DM_36d,M_mass_f,M_mass_f_sil,M_ST,M_33d,M_34d,M_36d,F_mass_f,F_ST,F_33d,F_34d,F_36d,SS_mass_f,SS_ST,SS_33d,SS_34d,SS_36d,S_ST,S_33d,S_34d,S_36d,rate_f,oscillate,resurf_cm_yr,sil_mag_S,thick_C,f_S2,f_pl2mo,f_sq,f_deep,f_remobilised,f_pu)


###################################################################
#                      ATMOSPHERE CALCULATION                     #
###################################################################
# # atmosphere
# profile = pd.read_csv("A_alt_T.csv")
# F_33R = 0.00912537740657392
# F_34R = 0.0584588083452021
# F_36R = 0.000239082372396731

# output = ie.calc_atm(profile,F_33R,F_34R,F_36R)
# print(output)

###################################################################
# ############### EXTRA CALCULATIONS #################
### range of fractionation factors at a single step ###
# ie.range_frac_fac_ss()

### range of fractionation factors for various proportions of space to burial ###
# x_hg = 0.2
# ie.range_frac_fac_s2b(x_hg)

## S isotope partitioning between core-mantle for various fractions in each #

## inputs
# d_i = -0.02 # bulk d34S of Io
# step_size = 0.01 # step size for fraction of material in core/mantle

# # run model
# ie.calc_iso_variable(d_i,step_size)