# Imports numpy to make arrays of values increasing incrementally
# Math is for the logarithmic steps in the creation of the density
# values, if necessary
import numpy as np
from math import *
from scipy.optimize import curve_fit

"""
A .py file containing used constants and values in the curve
creation/fitting process. The original script/file was written
by Andy Crump, a previous REU intern.

The comment to the right of each constant/variable shows the default
value in brackets, or [ ]. The notation of these comments is scientific
notation.
The comment directly above each constant/variable describes the constant
and possible values, if applicable.
The original name of each variable in Crump's script will be in
angled brackets, or < >.

Any variables with a double underscore, __, before their names are private
variables. They are not meant to be accessed from outside of this file.

Some of the comments are the same as Crump's due to either some of
the variables that were described being difficult/impossible to
decipher the meaning of or the variables being easy enough to understand.
Such comments are indicated by the following sequence of characters:
                                [*]

Unlike Crump's script, this does not include an external
.txt file to be read from. The variable names compared to Crump's
constants file are changed to enhance the understandability of the
theory being used to process the meteor data and generate a model that
fits the data.

This file's variables are imported into the scripts that require them.

Credit goes to Andy Crump, who set this up in MATLab.
"""

# Newton's gravitational constant (m^3/kg/s^2)
grav_const = 6.67408 * 10**(-11)

# Mean radius of the Earth (km)
r_e = 6371

# Mass of the Earth (kg)
m_e = 5.972 * 10**(24)

# =========================================================================

# US Standard Atmosphere Density at +80 km (kg/m^3)
# <rho80>
rho_80 = 1.846 * 10**(-5)           # [1.846 * 10^(-5)]

# -------------------------------------------------------------------------

# Atmospheric Pressure at sea level(N/m^2)
# <Pa>
pa_sea = 10.13 * 10**(4)            # [10.13 * 10^(4)]

# Heights to be used in creating a simple
# atmospheric pressure equation(km)
atmo_press_h = np.array([-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15,
                             20, 25, 30, 40, 50, 60, 70, 80])

# Pressures to be used in creating a simple atmospheric
# pressure equation (N/m^2)
atmo_press = np.array([11.39, 10.13, 8.988, 7.950, 7.012, 6.166,
                           5.405, 4.722, 4.111, 3.565, 3.080, 2.650,
                           1.211, 0.5529, 0.2549, 0.1197, 0.0287,
                           0.007978, 0.002196, 0.00052, 0.00011]) * 10**(4)

# Takes the log of the pressures and...
atmo_press_log = np.log10(atmo_press)

# Creates a best-fit line for the heights and log(pressures). Will be
# exponentiated to get back pressures.
p_a_params = np.polyfit(atmo_press_h, atmo_press_log, 1)

# -------------------------------------------------------------------------

# This block calculates the coefficients to the polynomials that are
# log fitted to the atmosphere's height dependent densities and
# temperatures.

# The following tables hold the atmospheric densities and temperatures
# at corresponding heights from the US Standard Atmosphere model and
# the MSISE-90 model.
# See the following website for the values:
# http://www.braeunig.us/space/atmos.htm

# Heights in km
heights = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
           12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80,
           100, 120, 140, 160, 180, 200, 220, 240, 260,
           280, 300, 320, 340, 360, 380, 400]

# Densities in kg/m^3
# NOTE: The change in formatting partway through is due to the
# two models differing in their formatting
densities = [1.47808, 1.347, 1.225, 1.11164, 1.00649, 0.909122,
             0.819129, 0.736116, 0.659697, 0.589501, 0.525168,
             0.466348, 0.412707, 0.310828, 0.193674, 0.0880349,
             0.0394658, 0.0180119, 0.00821392, 0.00385101,
             0.00188129, 0.000977525, 0.000288321, 0.0000742430,
             0.0000157005, 5.08e-7, 1.8e-8, 3.26e-9, 1.18e-9,
             5.51e-10, 2.91e-10, 1.66e-10, 9.91e-11, 6.16e-11,
             3.94e-11, 2.58e-11, 1.72e-11, 1.16e-11, 7.99e-12,
             5.55e-12, 3.89e-12]

# Temperatures in K
# The temperatures are split due to the fitting method not
# doing well with the values. low_temp goes from -2 to 100 km, inclusive.
# high_temp goes from 100 km to 400 km, inclusive.
low_temp = [301.5, 294.65, 288.15, 281.65, 275.15, 268.65, 262.15,
            255.65, 249.15, 242.65, 236.15, 229.65, 223.15,
            216.65, 216.65, 216.65, 221.65, 226.65, 237.05,
            251.05, 265.05, 270.65, 245.45, 217.45, 196.65,
            184.016]

# Temperatures in K
high_temp = [184.016, 374.9715, 635.5703, 787.5532, 877.6729,
             931.2806, 963.2701, 982.4191, 993.9173, 1000.8427,
             1005.0267, 1007.5620, 1009.1030, 1010.0423,
             1010.6166, 1010.9688]

# These arrays are holding their fitted polynomials' coefficients
log_dens_coeff = np.polyfit(heights, np.log10(densities), 7)
log_low_t_coeff = np.polyfit(heights[:26], np.log10(low_temp), 9)
log_high_t_coeff = np.polyfit(heights[25:], np.log10(high_temp), 9)

# -------------------------------------------------------------------------

# Avagadro's Number (particles/mol)
# <avagadroNum>
avagadro_num = 6.0221413 * 10**(23) # [6.02214113 * 10^(23)]

# Boltzmann Constant (J/K)
# <kB>
kb = 1.3806488 * 10**(-23)          # [1.3806488 * 10^(-23)]

# Stefan-Boltzmann Constant (W/m^2/K^4)
# <sigma>
stef_boltz = 5.670373 * 10**(-8)    # [5.670373 * 10^(-8)]

# =========================================================================

# Run script(s) with fewer parameters (debugging).
# True = few paramaters, False = all parameters
# <bShortRun>
b_short_run = False        # [False]

# Plot results after light curves are calculated.
# True = plot, False = don't plot
# <bPlot>
b_plot = True              # [False]

# Account for surviving photomass
# True = account for surviving photomass, False = don't account
# <bIterateMassScalar>
b_iter_mass_scalar = True  # [True]

# Show a progress percent of done iterations vs total iterations
# True = display, False = don't display
# <bDisplayPercentages>
b_display_percent = False  # [False]

# Start a new file, aka don't append to the previously made one
# True = start a new one, False = append
# <bStartFileOver>
b_start_new_file = False   # [True]

# Take every nth data point. Higher n = higher speed
# <stepSize>
step_size = 1              # [1]

# Maximum standard deviation allowed for magnitude
# <maxSigma>
max_sd = 1.0               # [1.0]

# =========================================================================

# Radius meteor is seen at (m)
# <r>
r = 100 * 10**(3)      # [100 * 10^(3)]

# Range of wavelengths (nm)
# <wavelengthRange>
wavelength_range = 380 # [380]

# Scalar luminous efficiency (2003, Campbell-Brown and Koschny, 2003, p.755)
# <tau>
lum_eff = 2 * 10**(-3) # [2 * 10^(-3)]

# Base intensity (Vega's intensity) (W/m^2/nm)
# <I0>
I_0 = 3.67 * 10**(-11) # [3.67 * 10^(-11)]

# =========================================================================

# Atmospheric scale height (km) [*]
# <Hatm>
h_atm = 7.4       # [7.4]

# Atmospheric temperature (K) [*]
# <Ta>
# atmo_temp = 200    # [200]

# Starting temperature of the meteoroid (K)
# <Tinf>
initial_temp = 270  # [270]

# Shape factor. [*]
# 1 =  spherical, > 1 = not spherical
# <A>
shape_fact = 1.21 # [1.21]

# Drag coefficient [*]
# Ranges from 0 to 2 (Campbell-Brown and Koschny, 2003, p.754)
# <gamma>
drag_coeff = 1.0  # [1.0]

# Starting (interface) height (km) [*]
# <hai>
h_ai = 200        # [200]

# Vapor pressure on the meteoroid, negligible (Campbell-Brown and Koschny,
#                                              2003*, p.16)
# <Pv>
vapor_press = 0.0 # [0.0]

# =========================================================================

# Determine which step to use (NOTE: MAY NOT BE USED IN CRUMP'S SCRIPT(S))
# True = logarithmic, False = normal
# <bRhoMMLogStep>
density_log_step = True # [True]

# Starting possible meteor density (kg/m^2) (?)
# <rhoMMStart>
__density_meteor = 100 # [100]

# Density step size (kg/m^2) (?)
# <rhoMMStep>
__density_step = 0.2   # [0.2]

# Ending possible meteor density (kg/m^2) (?)
# <rhoMMEnd>
__density_end = 10000  # [10000]


# Generates the values between the starting and ending densities
# with the specified step size. Dependent on whether or not
# a logarithmic step is specified(kg/m^2) (?)

# Initializes the density array
# <rhoMMLoop>
density_array = []

# Creates log scaled density values if specified. A lot fewer than a normally
# stepped list of density values
if(density_log_step):
    # Creates the interval on a logarithmic step. An additional "step" is added
    # since arange does NOT include the specified end value
    log_density_array = np.arange(log10(__density_meteor),
                                  log10(__density_end) + __density_step,
                                  __density_step)
    # Lops off the final value if it is greater than the end value. This is a
    # "just in case", as floats are a little weird with arange().
    # The 10^(-10) is there because of how Python handles floats
    if log_density_array[-1] > (log10(__density_end) + 10**(-10)):
        log_density_array = log_density_array[:-1]

    # Calculate the actual density values on the log scale
    for i in range(len(log_density_array)):
        density_array.append(10**(log_density_array[i]))

# Just generate the density values normally if a log step is NOT specified.
else:
    density_array = np.arange(__density_meteor, (__density_end + __density_step),
                              __density_step)
    if density_array_[-1] > (__density_end + 10**(-10)):
        density_array = density_array[:-1]

# =========================================================================

# Lowest heat of ablation (J/kg) (Kikwaya, 2011, p.172)
# <LStart>
__heat_of_ablation = 2 * 10**(6) # [2 * 10^(6)]

# Heat of ablation linear step size (J/kg)
# <LStep>
__heat_abl_step = 1 * 10**(6)    # [1 * 10^(6)]

# Highest heat of ablation (J/kg)
# <LEnd>
__heat_abl_end = 9 * 10**(6)     # [9 * 10^(6)]

# Generates the values between the lowest and highest heats of ablation
# with the specified step size (J/kg)
# <LLoop>
heat_abl_array = np.arange(__heat_of_ablation,
                           (__heat_abl_end + __heat_abl_step),
                           __heat_abl_step)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if heat_abl_array[-1] > (__heat_abl_end + 10**(-10)):
    heat_abl_array = heat_abl_array[:-1]

# =========================================================================

# Lowest "glue" release of grains temperature (K)
# <TglueStart>
__glue_temperature = 900 # [900]

# "Glue" release temperature step size (K)
# <TglueStep>
__glue_step = 100        # [100]

# Highest "glue" release temperature (K)
# <TglueEnd>
__glue_end = 1600        # [1600]

# Generates the values between the lowest and highest glue release
# temperatures with the specified step size (K)
# <TglueLoop>
glue_array = np.arange(__glue_temperature,
                       (__glue_end + __glue_step), __glue_step)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if glue_array[-1] > (__glue_end + 10**(-10)):
    glue_array = glue_array[:-1]

# =========================================================================

# Lowest boiling temperature of meteoroid(K)
# <TbStart>
__boil_temperature = 1400 # [1400]

# Boiling temperature step size (K)
# <TbStep>
__boil_step = 200         # [200]

# Highest boiling temperature (K)
# <TbEnd>
__boil_end = 2300         # [2300]

# Generates the values between the lowest and highest meteoroid boiling
# temperatures with the specified step size (K)
# <TbLoop>
boil_array = np.arange(__boil_temperature,
                       (__boil_end + __boil_step), __boil_step)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if boil_array[-1] > (__boil_end + 10**(-10)):
    boil_array = boil_array[:-1]

# =========================================================================

# Lowest thermal conductivity of meteoroid (J/(m*s*K))
# <lambdaCStart>
__thermal_conduct = 0.5 # [0.5]

# Thermal conductivity step size (J/(m*s*K))
# <lambdaCStep>
__thermal_step = 2.5 # [2.5]

# Highest thermal conductivity of meteoroid (J/(m*s*K))
# <lambdaCEnd>
__thermal_end = 3.0  # [3.0]

# Generates the values between the lowest and highest thermal conductivity
# values with the specified step size (J/(m*s*K))
# <lambdaCLoop>
thermal_array = np.arange(__thermal_conduct, (__thermal_end + __thermal_step),
                          __thermal_step)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if thermal_array[-1] > (__thermal_end + 10**(-10)):
    thermal_array = thermal_array[:-1]

# =========================================================================

# Lowest atomic mass of vapor per atom (Campbell-Brown and Kronchsky, 2003*,
# p.19) (amu)
# Ranges from 12 - 56 amu
# <muStart>
__atomic_m_vapor = 23 # [23]

# Atomic mass of vapor step size (amu)
# <muStep>
__atomic_m_step = 1.0 # [1.0]

# Highest atomic mass of vapor (amu)
# <muEnd>
__atomic_m_end = 23   # [23]

# Generates the values between the lowest and highest atomic mass of vapor
# of the meteoroid with the specified step size (amu)
# <muLoop>
atomic_m_array = np.arange(__atomic_m_vapor,
                           (__atomic_m_end + __atomic_m_step),
                           __atomic_m_step)
# amu to kg
atomic_m_array = atomic_m_array * (10**(-3) / avagadro_num)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if atomic_m_array[-1] > (__atomic_m_end + 10**(-10)):
    atomic_m_array = atomic_m_array[:-1]

# =========================================================================

# Specific heat of meteoroid (J/kg/K) [*]
# Ranges from 600 to 1400 J/kg/K
# <cStart>
__spec_heat_meteor = 600 # [600]

# Specific heat step size (J/kg/K)
# <cStep>
__spec_heat_step = 200   # [200]

# Highest specific heat (J/kg/K)
# <cEnd>
__spec_heat_end = 1400   # [1400]

# Generates the values between the lowest and highest specific heat values
# of the meteoroid with the specified step size
# <cLoop>
spec_heat_array = np.arange(__spec_heat_meteor,
                            (__spec_heat_end + __spec_heat_step),
                            __spec_heat_step)

# Checks the last array value to see if it is greater than the specified
# end value. If so, lop it off.
if spec_heat_array[-1] > (__spec_heat_end + 10**(-10)):
    spec_heat_array = spec_heat_array[:-1]

# =========================================================================

# If a short run time is desired (used for debugging), use
# these parameters instead. Each variable points to an array of
# values over which to iterate over. The best value from each of
# the arrays is chosen to form the model of a specific meteor.
# Change 
if(b_short_run):
    print("DEBUGGING: RUNNING WITH FEWER PARAMETERS")
    print("SET b_short_run IN constants.py TO FALSE IF NOT DESIRED.")
    # Possible meteor densities
    #density_array = [2000]
    density_array = np.arange(2000, 5000 + 500, 500)
    # Possible meteor heats of ablation
    #heat_abl_array = [2 * 10**(6)]
    heat_abl_array = np.arange(1 * 10**6, 9 * 10**6, 2*10**6)
    # Possible temperatures at which the meteor glue "fails"
    #glue_array = [900]
    glue_array = np.arange(1000, 1600, 200)
    # Possible meteor boiling temperatures
    boil_array = np.arange(1500, 2300 + 200, 400)
    # Possible meteor thermal conductivity values
    thermal_array = [.5, 3]
    # Possible average meteor vapor atomic masses
    atomic_m_array = [23 * 10**(-3)/avagadro_num]
    # Possible specific heat of the meteor
    spec_heat_array = [600]
    #spec_heat_array = np.arange(600, 1400 + 200, 200)

# =========================================================================

# Sticking factor (Campbell-Brown and Koschny, 2013, p.6) [*]
# 0.5 = stone, 1 =  metals
# <phi>
sticking_fact = 0.5  # [0.5] 

# Scalar meteoroid emissivity
# <epsilon>
met_emissivity = 0.9 # [0.9]

# =========================================================================

# Heat transfer coefficient of the meteoroid material(s)
# Ranges from 0.5 to 6.0
# (https://www.spaceacademy.net.au/watch/debris/metflite.htm)
# <Lambda>
heat_trans_coeff = 1.0           # [1.0]

# A constant that describes the fraction of kinetic energy lost
# due to mass ablation. Used in calculating the change of
# temperature. Ranges from 0 to 1.0 (Kikwaya, 2011 p.171)
heat_loss_coeff = 0.5            # [0.5]

# Smallest mass the meteor breaks up into (kg)
# <minMass>
min_breakup_mass = 1 * 10**(-10) # [1 * 10^(-7)]

# Multiply the photomass by this to account for surviving mass
# 1.0 = all mass survives
# <massScalar>
mass_scalar = 1.0                # [1.0]

# Considered step-size in log mass (used in breakup_model) [*]
# <dLogM>
d_log_mass = 0.1                 # [1.0]

# =========================================================================

# Plot vertical line at breakups [*]
# 0 = off
# <iBreakupLineHeight>
breakup_line_height = 0.1 # [0.1]

# Plot projections before and after shown magnitudes
# <bPreplot>
preplot = False           # [False]

# Ignore maximum magnitude and 0 velocity to calculate residuals
# <bIgnoreMins>
ignore_mins = False       # [False]

# Maximum magnitude allowed aka minimum brightness
# <maxMag>
max_magni = 20            # [20]

# =========================================================================

# File number to start on (not meteor number)
# <startNum>
start_num = 1 # [1]

# File number to end on
# 0 runs through all meteors in the folder
# <endNum>
end_num = 0   # [0]

# Maximum number of topN data points to consider
# <maxN>
max_n = 40    # [40]

# =========================================================================

# Dictionary holding values (from Crump's script) subject to change.
# This is a substitute to MATLab's structs in holding these sorts
# of things
K = {"h_atm":h_atm, "initial_temp":initial_temp,
     "shape_fact":shape_fact, "drag_coeff":drag_coeff, "h_ai":h_ai,
     "vapor_press":vapor_press, "sticking_fact":sticking_fact,
     "met_emissivity":met_emissivity, "heat_trans_coeff":heat_trans_coeff,
     "min_breakup_mass":min_breakup_mass, "mass_scalar":mass_scalar,
     "d_log_mass":d_log_mass}

