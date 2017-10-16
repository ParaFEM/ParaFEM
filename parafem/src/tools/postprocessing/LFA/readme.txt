LFA_CALC_K

Reads the .ttr2 results produced by xx12 and calculates the effective thermal 
conductivity from a laser flash analysis (LFA) simulation.

    Usage: LFA_calc_k <job_name>
    Expects: <job_name>.py in current working directory
    
Add path below to .bash_profile (or similar)
PYTHONPATH=~/PATH_TO_PARAFEM/parafem/bin

<job_name>.py format:
l_mm = float #mm
r_mm = float #mm
m    = float #kg
Cp   = float #J/kg K
Tmax = float #degC

Leave Tmax=0.0 if ttr2 is complete (i.e. simulation has reached equilibrium),
otherwise set Tmax to calculated equilibrium T value.

Author: Llion Evans
Date:   16 October 2017

Copyright University of Manchester 2017

