# Imports numpy for mathematical arrays and array-based operations
import numpy as np
# Imports mathematical constants and functions
from math import *
# Imports os for path manipulation
import os

import matplotlib.pyplot as plt
# Imports constants and needed values
from lc_constants import *
# Imports the Particles class to make use of it in the breakup_model and
# in here
from lc_particles import Particles
# Imports breakup_model to simulate the breakup of particles into
# other particles
from lc_breakup_model import breakup_model
# Imports erode_particle to simulate the breakup up of particles of the
# main meteoroid body due to erosion
from lc_erode_particle import erode_particle


def get_mass(constants_change, bool_plot, photomass, time, flux_smooth,
             vel_time, vel, heights, zenith, step_size, meteor, working_dir,
             camera_names, avged_fluxes, avged_heights):
    
    """
    The name is more or less a misnomer. Doesn't actually return the mass
    of the meteor. Rather, it churns out the time vs magnitude model and
    the time vs velocity model. From what I can tell, this is called within
    an enormously nested series of loops.
    This script also puts out the height vs magnitude, height vs velocity,
    and height vs temperature of the meteoroid. In the case of mass, the
    total mass of the meteoroid is calculated while temperature is taken
    from the main body of the meteoroid.
    The magnitude information was originally plotted within
    flux_time_driver.py along with the corresponding times,
    but it was determined that in order to plot the magnitude versus the
    height it would be easiest to do the plotting in this script.

    Some of the used equations have been put into their own
    subfunctions for the sake of better organization and ease of editing

    Credit goes to Andy Crump for his getMass36.m MATLab script.

    Parameter(s):
    constants_change - A dictionary of values that change due to being in an
                       enormously nested series of loops.
    bool_plot        - A boolean detailing whether or not to plot. This is
                       not implemented
    photomass        - The photometric mass of the meteor
    time             - An array of time values for the meteor from
                       flux_time_driver()
    flux_smooth      - An array of flux values of the meteor from
                       flux_time_driver(). Each flux value corresponds to
                       each time value by index
    vel_time         - The complete time values. Similar to time, but
                       corresponds only to the velocities in vel. Only
                       used in a single line of code to create a
                       line of best fit
    vel              - The complete velocity values. Corresponds to
                       vel_time. Only used in a single line of code
                       to find a standard deviation
    heights          - All of the height values of the meteor
    zenith           - The zenith angle of the meteor
    step_size        - The step size used to calculate various values.
                       For a step_size n, only every n values will be used.
                       Primarily used for speed
    meteor           - The current meteor's number
    working_dir      - The current working directory
    camera_names     - The names of the cameras that detected the meteor
    avged_fluxes     - The meteor's magnitude as detected. An array with
                       nested arrays
    avged_heights    - The meteor's heights detected by the cameras. An
                       array with nested arrays representing the different
                       cameras

    Return:
    Returns the flux and velocity standard deviations of the created
    models as floats. 
    """
    
    # Creates a Particles object to represent a system of particle types
    particle_sys = Particles()
    
    # Most of the below values change as the function gets called from
    # within a loop. Constants such as the atmospheric scale height are
    # already inside of the imported constants.py

    erosion_coeff = constants_change["erosion_coeff"]
    
    # The density of a particle. NOT THE COMPOSITE METEOROID DENSITY
    grain_dens = constants_change["grain_dens"]
    # The specific heat of the particle
    spec_heat = constants_change["spec_heat"]
    # The thermal conductivity of the particle
    thermal_cond = constants_change["thermal_cond"]
    # Boiling temperature of the meteor material
    boil_temp = constants_change["boil_temp"]
    # Temperature at which material on the meteor will fuse together
    fusion_temp = constants_change["fusion_temp"]
    # Atomic mass per atom of vapor
    atomic_m = constants_change["atomic_m"]
    # Specific heat of the meteor
    spec_heat = constants_change["spec_heat"]
    # Density of the entire meteor
    bulk_dens = constants_change["bulk_dens"]
    # Heat of ablation aka heat at which meteor material on the surface
    # starts vaporizing
    heat_abl = constants_change["heat_abl"]
    # Temperature at which the "glue" within the meteor dissipates,
    # resulting in the other materials held by the glue to no longer
    # be bound
    glue_temp = constants_change["glue_temp"]
    
    
    # The linear regression of the Murray and Beech values for the
    # mass distribution of the mass bins of the particles of a
    # meteoroid. Definition is as follows:
    # s = -2.2735*Fvalue + 2.907
    s = constants_change["s"]

# ==========================================================================

    # Threshold mass. If a calculated mass is lower than this, it is treated
    # differently in the calculations
    ignore_mass = min_breakup_mass / 10
    # Standard deviations of the measured fluxes compared to the generated
    # model, respectively
    flux_s_dev = inf
    vel_s_dev = inf

    # Pre-atmosphere velocity
    vel_inf = constants_change["vel_inf"]

    # The number of particle types in the particle system
    p_len = 1

    # Stores the calculated magnitudes and velocities. It looks like
    # dynamic_vel is used to keep track of the velocity of the largest
    # particle from the breakup_model() output
    dynamic_mag = []
    dynamic_vel = []
    
    # Creates the very first mass bin. The single particle in the mass
    # represents the "core" of the meteoroid
    particle_sys.new_particle(photomass, 1, vel_inf, initial_temp, bulk_dens)

    # Makes a bestfit line between the heights and times passed into
    # this script
    line = np.polyfit(heights, vel_time, 1)

    # Interface time
    inter_time = line[0] * h_ai + line[1]

    # No idea why this variable is here. It'd be very slightly more
    # efficient if inter_time is used
    pre_time = inter_time
    # Makes a small time interval to calculate small changes in particle
    # type velocity and mass
    dt = 0.05
    # Starting time. Will increase by 1 per loop of the outermost loops
    t = 0

    """
    # No idea why this is here
    presteps = ceil((time[0] - inter_time) / dt)
    """

    # Holds the total magnitude of the entire particle system representing
    # the fragmenting meteor
    pre_mag = []

    # Appears to be a remnant of Crump's plotting stuff
    pre_times = []


# ==========================================================================
#                   HIGH ELEVATION BEFORE METEOR IS DETECTED
# ==========================================================================

    # m, n, v, t, and d refer to the mass, number of particles, velocity,
    # temperature, and grain densities of each mass bin.
    # As of the date of modification, these lists were used to plot
    # those values in order to debug/modify Andy Crump's algorithms
    # These particular lists hold the values when the meteoroid is high up
    # in the atmosphere before it is detected
    holder0_m = []
    holder0_n = []
    holder0_v = []
    holder0_t = []
    holder0_d = []
    
    # Runs the equations for the upper atmosphere to find the meteor's
    # initial temperature when it is seen
    while pre_time < time[0]:

        # If the particle's temperature is higher than the glue failing
        # temperature, it is going to break up into smaller particles.
        # Two models: the disruption or an erosion model
        # If set to False, the disruption model will not run.
        if(particle_sys.get_value("TEMPERATURE", 0) > glue_temp and False):
            (p_len, particle_sys) = breakup_model(particle_sys, p_len, 0,
                                                  constants_change, zenith)

        # Saves the amount of mass to be eroded from the meteoroid
        eroded_mass = 0

        # Holds the total intensity of all of the particles
        sum_intensity = 0
        # For each type of particle...
        for p in range(p_len):

            # Get the particle type mass
            particle_mass = particle_sys.get_value("MASS", p)

            # If the particle type mass < the amount of mass to be ignored,
            # move on to the next particle type
            if(particle_mass < ignore_mass):
                continue

            # Get the particle type temperature, velocity, and density
            particle_temp = particle_sys.get_value("TEMPERATURE", p)
            particle_vel = particle_sys.get_value("VELOCITY", p)
            particle_density = particle_sys.get_value("DENSITY", p)
            
            # Gets the current height of the particle
            curr_height = (pre_time - line[1]) / line[0]

            # Gets the atmospheric density at this current height on Earth
            curr_atmo_dens = 10**np.polyval(log_dens_coeff, curr_height)

            # Gets the atmospheric temperature at the current height.
            # To be used in calculating d_parti_temp
            atmo_temp = 0
            if(curr_height >= 100):
                atmo_temp = 10**np.polyval(log_high_t_coeff, curr_height)
            else:
                atmo_temp = 10**np.polyval(log_low_t_coeff, curr_height)

            # ---------------------------------------------------------
            # Calculate dv

            # Gets the change in velocity over a small time increment dt.
            # This equation describes the momentum transfer between the
            # particle and air particles hitting it
            # (Kikwaya, 2011, p.15)
            dv = vel_change(dt, curr_atmo_dens, particle_vel,
                            particle_mass, particle_density)

            """
            # Changes the particle type velocity
            particle_vel -= dv
            """

            # ---------------------------------------------------------
            # Calculate dm

            # Calculates the pressure of the saturated vapor from the particle
            # on the particle. Used to calculate the changes in mass and
            # temperature below
            satur_vapor_press = sat_vapor_press(curr_height, heat_abl, atomic_m,
                                                boil_temp, particle_temp,
                                                curr_atmo_dens)

            # Prevents the saturated vapor pressure from being negative
            if(satur_vapor_press < 0):
                satur_vapor_press = 0

            # Calculates the change in mass over a small increment dt
            # (Kikwaya, 2011, p.23)
            
            dm = small_mass(dt, curr_atmo_dens, heat_abl, particle_mass,
                            particle_density, particle_vel, satur_vapor_press,
                            particle_temp, atomic_m)

            # Makes sure that the meteoroid doesn't gain mass
            if (dm < 0):
                dm = 0

            # Temporarily holds the eroded mass value. Used in calculations
            eroded_mass_1 = 0

            # If the mass bin holds the main particle, calculate
            # the amount of mass to be eroded
            if (p==0 and particle_temp >= glue_temp):
                eroded_mass = erosion_mass(dt, erosion_coeff, particle_mass,
                                           particle_density, particle_vel,
                                           curr_atmo_dens)
                eroded_mass_1 = eroded_mass
            else:
                eroded_mass_1 = 0

            # Changes the particle type mass and if the resulting mass is
            # negative, sets it to 0 instead.
            # dm is the ablated mass
            particle_mass -= (dm + eroded_mass_1)

            # Prevents particles from having negative masses
            if(particle_mass < 0):
                particle_mass = 0

            # ---------------------------------------------------------
            # Calculate dT

            # Calculates the change in temperature over a small increment dt.
            # (Kikwaya, 2011, p.22)
            d_parti_temp = small_temp(dt, spec_heat, particle_mass,
                                      curr_atmo_dens, particle_vel,
                                      particle_density, particle_temp,
                                      heat_abl, dm - eroded_mass_1, atmo_temp)

            if(particle_temp > atmo_temp and
               (particle_temp + d_parti_temp) < atmo_temp):
                particle_temp = atmo_temp
            elif(particle_temp < atmo_temp and d_parti_temp < 0):
                particle_temp = particle_temp
            else:
                # Changes the particle type temperature
                particle_temp += d_parti_temp
                
            # ---------------------------------------------------------
            # Calculate luminosity

            # Calculates the energy emitted from the particle type
            energy_em = dm * lum_eff * particle_vel**(2) / 2 \
                        + lum_eff * particle_mass * dv * particle_vel
            
            # Calculates the observed intensity of teh particle 
            I_obs = energy_em / (dt * 4 * pi * r**(2) * wavelength_range)
            # Adds the observed intensities of ALL of the particles of that
            # particle type
            sum_intensity += I_obs * particle_sys.get_value("PARTICLE_N", p)

            #//////////////////////////////////////////////////////////
            """
            # Changes the particle type velocity
            particle_vel -= dv

            # Changes the particle type mass and if the resulting mass is
            # negative, sets it to 0 instead.
            # dm is the ablated mass
            particle_mass -= (dm + eroded_mass_1)
            

            if(particle_temp > atmo_temp and
               (particle_temp + d_parti_temp) < atmo_temp):
                particle_temp = atmo_temp
            elif(particle_temp < atmo_temp and d_parti_temp < 0):
                particle_temp = particle_temp
            else:
                # Changes the particle type temperature
                particle_temp += d_parti_temp
            """
            #//////////////////////////////////////////////////////////
            
            # Actually does the aforementioned changes
            particle_sys.change_values(p, particle_mass, "SAME", particle_vel,
                                       particle_temp, "SAME")

            # ---------------------------------------------------------

            # Puts the calculated values in the holder lists
            if (len(holder0_m) <= p or not holder0_m):
                while(len(holder0_m) <= p):
                    holder0_m.append([])
                    holder0_n.append([])
                    holder0_v.append([])
                    holder0_t.append([])
            holder0_m[p].append(particle_mass)
            holder0_n[p].append(particle_sys.get_value("PARTICLE_N", p))
            holder0_v[p].append(particle_vel)
            holder0_t[p].append(particle_temp)

        # Appends the magnitude observed to the list holding all of
        # the magnitude data
        # sum_intensity/I_0 is used in the check since log10 doesn't spit out
        # -infinity like MATLab does
        if((sum_intensity / I_0) > 0):
            if((-2.5 * log10(sum_intensity / I_0)) < max_magni):
                pre_mag.append(-2.5 * log10(sum_intensity / I_0))
            else:
                pre_mag.append(max_magni)
        # If the calculated magnitude is larger than the allowed magnitude
        # set it to the allowed magnitude
        else:
            pre_mag.append(max_magni)

        # Actually does the erosion of the main particle. Takes off
        # the calculated eroded_mass from above if it isn't equal to 0.
        if (eroded_mass != 0):
            particle_sys = erode_particle(particle_sys, constants_change,
                                          eroded_mass)
            p_len = particle_sys.get_bin_num()
        

        # Append the used pre-time to the list (I have no idea what
        # this does)
        pre_times.append(pre_time)
        
        # Increments pre-time by the interval dt
        pre_time += dt

        # Increments t by 1
        t += 1

    # Making the velocity of the first mass bin's particle it's original
    # detected velocity again
    for i in range(p_len):
        particle_sys.change_values(i, "SAME", "SAME", vel_inf, "SAME", "SAME")

# ==========================================================================
#                       ELEVATION WHERE METEOR IS DETECTED
# ==========================================================================
    heights[::-1].sort()
    # m, n, v, t, and d refer to the mass, number of particles, velocity,
    # temperature, and grain densities of each mass bin.
    # As of the date of modification, these lists were used to plot
    # those values in order to debug/modify Andy Crump's algorithms
    # These particular lists hold the values when the meteoroid is in the
    # atmosphere and is being detected by the cameras
    holder_m = []
    holder_n = []
    holder_v = []
    holder_t = []
    holder_d = []

    # Tracks the main meteoroid body's velocity
    meteoroid_vel = [vel_inf]

    # Tracks the main meteoroid body's temperature
    meteoroid_temp = [particle_sys.get_value("TEMPERATURE", 0)]
    
    # Crump calls this large loop the "main mass equations." My best guess
    # is that it does what the above large loop does but for when the
    # meteor is in the Earth's atmosphere
    for t in range(0, (time.size-1), step_size):

        # If the particle's temperature is higher than the glue failing
        # temperature, it is going to break up into smaller particles.
        # As such, run breakup_model()
        if(particle_sys.get_value("TEMPERATURE", 0) > glue_temp and False):
            (p_len, particle_sys) = breakup_model(particle_sys, p_len, 0,
                                                  constants_change, zenith)

        # Saves the amount of mass to be eroded from the meteoroid
        eroded_mass = 0

        # Holds the total intensity of all of the particles
        sum_intensity = 0

        # Holds the "index" of the largest particle by mass
        largest_particle = 0

        # Loop looks for the above-mentioned largest particle
        for p in range(1, p_len):
            if(particle_sys.get_value("MASS", p) > \
               particle_sys.get_value("MASS", p) and \
               particle_sys.get_value("PARTICLE_N", p) >= 1):
                largest_particle = p

        # Initializes variables to hold the total particle type mass
        # and the specific "tiny mass"
        total_mass = 0
        tiny_mass = 0

        # For each particle type...
        for p in range(p_len):
            
            # Get the total mass of the combined particle types
            total_mass += particle_sys.get_value("MASS", p) * \
                          particle_sys.get_value("PARTICLE_N", p)

            # If the total mass of a certain particle type is less than
            # the amount of mass to be ignored, add that to the tiny mass
            # tracker and change the mass of each particle of that particle
            # type to 0
            if(particle_sys.get_value("MASS", p) < ignore_mass):
                tiny_mass += particle_sys.get_value("MASS", p) * \
                             particle_sys.get_value("PARTICLE_N", p)
                particle_sys.change_values(p, 0, "SAME", "SAME", "SAME",
                                           "SAME")
                
        # For each particle type...
        for p in range(p_len):
            # If the mass of one particle of a particle type is at least
            # the mass to be ignored, normalize it according to the
            # total and tiny masses
            if(particle_sys.get_value("MASS", p) >= ignore_mass):
                # Adds a check to see if the tiny mass is equal to the
                # total mass to prevent division by 0.
                if((total_mass - tiny_mass) <= 0):
                    continue
                else:
                    normalized_mass = particle_sys.get_value("MASS", p) * \
                                      total_mass / (total_mass - tiny_mass)
                    particle_sys.change_values(p, normalized_mass, "SAME",
                                           "SAME", "SAME", "SAME")
                

        # For each type of particle...
        for p in range(p_len):
            
            # Get the particle type mass
            particle_mass = particle_sys.get_value("MASS", p)

            # Gets the particle type temperature
            particle_temp = particle_sys.get_value("TEMPERATURE", p)

            # If the particle type mass < the amount of mass to be ignored,
            # move on to the next particle type
            if(particle_mass < ignore_mass):
                continue

            # Get the particle type velocity
            particle_vel = particle_sys.get_value("VELOCITY", p)

            # If the speed being recorded is from the largest particle,
            # save the speed
            if(p == largest_particle):
                dynamic_vel.append(particle_vel)

            # Gets the particle type density
            particle_density = particle_sys.get_value("DENSITY", p)


            # Gets the atmospheric density at a height on Earth
            curr_atmo_dens = 10**np.polyval(log_dens_coeff, heights[t])

            # Gets the atmospheric temperature at the current height.
            # To be used in calculating d_parti_temp
            atmo_temp = 0
            if(heights[t] >= 100):
                atmo_temp = 10**np.polyval(log_high_t_coeff, heights[t])
            else:
                atmo_temp = 10**np.polyval(log_low_t_coeff, heights[t])

            # Calculates a different time interval than the previous
            # one in the previous section. This one is based on the different
            # between two "adjacent" times
            dt = time[t+1] - time[t]

            # ---------------------------------------------------------
            # Calculate dv

            # Gets the change in velocity over a small time increment dt.
            # This equation describes the momentum transfer between the
            # particle and air particles hitting it
            # (Kikwaya, 2011, p.15)
            dv = vel_change(dt, curr_atmo_dens, particle_vel,
                            particle_mass, particle_density)

            # Changes the particle type velocity
            particle_vel -= dv

            # ---------------------------------------------------------
            # Calculate dm

            # Calculates the pressure of the saturated vapor from the particle
            # on the particle. Used to calculate the changes in mass and
            # temperature below
            satur_vapor_press = sat_vapor_press(heights[t], heat_abl, atomic_m,
                                                boil_temp, particle_temp,
                                                curr_atmo_dens)

            # Prevents negative vapor pressures
            if(satur_vapor_press < 0):
                satur_vapor_press = 0

            # Calculates the change in mass over a small increment dt
            dm = small_mass(dt, curr_atmo_dens, heat_abl, particle_mass,
                            particle_density, particle_vel, satur_vapor_press,
                            particle_temp, atomic_m)

            # Changes the mass bin mass and if the resulting mass is
            # negative, sets it to 0 instead
            particle_mass -= (dm + eroded_mass_1)
            if(particle_mass < 0):
                particle_mass = 0

            # Makes sure that the meteoroid doesn't gain mass
            if (dm < 0):
                dm = 0

            # Temporarilty holds the eroded mass value. Used in calculations
            eroded_mass_1 = 0

            # If the mass bin holds the main particle, erode that particle
            if (p==0 and particle_temp >= glue_temp):
                eroded_mass = erosion_mass(dt, erosion_coeff, particle_mass,
                                           particle_density, particle_vel,
                                           curr_atmo_dens)
                eroded_mass_1 = eroded_mass
            else:
                eroded_mass_1 = 0

            # ---------------------------------------------------------
            # Calculate dT

            # Calculates the change in temperature over a small increment dt.
            d_parti_temp = small_temp(dt, spec_heat, particle_mass,
                                      curr_atmo_dens, particle_vel,
                                      particle_density, particle_temp,
                                      heat_abl, dm - eroded_mass_1, atmo_temp)
            
            if(particle_temp > atmo_temp and
               (particle_temp + d_parti_temp) < atmo_temp):
                particle_temp = atmo_temp
            elif(particle_temp < atmo_temp and d_parti_temp < 0):
                particle_temp = particle_temp
            else:
                # Changes the particle type temperature
                particle_temp += d_parti_temp

            # ---------------------------------------------------------
            # Calculate luminosity

            # Calculates the energy emitted from the particle type
            energy_em = dm * lum_eff * particle_vel**(2) / 2 \
                        + lum_eff * particle_mass * dv * particle_vel
            
            # Calculates the observed intensity of teh particle 
            I_obs = energy_em / (dt * 4 * pi * r**(2) * wavelength_range)
            
            # Adds the observed intensities of ALL of the particles of that
            # particle type
            sum_intensity += I_obs * particle_sys.get_value("PARTICLE_N", p)

            #//////////////////////////////////////////////////////////
            """
            # Changes the particle type velocity
            particle_vel -= dv

            # Changes the particle type mass and if the resulting mass is
            # negative, sets it to 0 instead.
            # dm is the ablated mass
            particle_mass -= (dm + eroded_mass_1)

            if(particle_temp > atmo_temp and
               (particle_temp + d_parti_temp) < atmo_temp):
                particle_temp = atmo_temp
            elif(particle_temp < atmo_temp and d_parti_temp < 0):
                particle_temp = particle_temp
            else:
                # Changes the particle type temperature
                particle_temp += d_parti_temp
            """

            #//////////////////////////////////////////////////////////
    
            # ---------------------------------------------------------
            # The values saved in this section are plotted in the end
            # plot

            # Tracks the main meteoroid body's temperature
            if(p==0):
                meteoroid_temp.append(particle_temp)

            # Tracks the main meteoroid body's velocity
            if(p==0):
                meteoroid_vel.append(particle_vel)

            # ---------------------------------------------------------

            # Saves the calculated values to the holder lists
            if (len(holder_m) <= p or not holder_m):
                while(len(holder_m) <= p):
                    holder_m.append([])
                    holder_n.append([])
                    holder_v.append([])
                    holder_t.append([])
            holder_m[p].append(particle_mass)
            holder_n[p].append(particle_sys.get_value("PARTICLE_N", p))
            holder_v[p].append(particle_vel)
            holder_t[p].append(particle_temp)

            # Actually does the aforementioned changes
            particle_sys.change_values(p, particle_mass, "SAME", particle_vel,
                                       particle_temp, "SAME")
        
        # Calculates the total magnitude of the particle system. If it is
        # larger than the maximum allowed magnitude, it is set to that max.
        # It is then appended to the end of the list initialized way back
        # above that holds the calculated magnitudes
        if((sum_intensity / I_0) > 0):
            if((-2.5 * log10(sum_intensity / I_0)) <= max_magni):
                dynamic_mag.append(-2.5 * log10(sum_intensity / I_0))
            else:
                dynamic_mag.append(max_magni)
        else:
            dynamic_mag.append(max_magni)

        # Actually does the erosion of the main particle. Takes off
        # the calculated eroded_mass from above if it isn't equal to 0.
        if (eroded_mass != 0):
            particle_sys = erode_particle(particle_sys, constants_change,
                                          eroded_mass)
            p_len = particle_sys.get_bin_num()

# ==========================================================================
    # This following block finds the standard deviation of the magnitude

    # Gets the total number of calculated magnitudes
    end_mag = len(dynamic_mag)

    # Ignores the velocity = 0 and max_magni restrictions (?)
    if(ignore_mins):
        for i in range(end_mag, -1, -1):
            if(dynamic_mag[i] < max_magni):
                end_mag = i
                break

    # ------------------------------------------------------------
    # This section plots the held values in the holder lists. As of
    # the last date of modification, it is being used as a way
    # to look at the various values held by each mass bin in order
    # to debug/modify existing algorithms.
    # The location of this is in a not-so-great spot
    # Remove/add the pound symbols next to the triple quotes above
    # and below a section to comment out/in that section.
    
    if(bool_plot):

        # 
        save_dir = os.path.join(working_dir, "LightcurvePlots")
        
        # Magnitude vs height figure
        fig = plt.figure(1)
        plt.plot(heights[:len(dynamic_mag)], dynamic_mag,
                 "k-", label="Best-Fit Model")
        for i in range(len(avged_fluxes)):
            plt.plot(avged_heights[i], avged_fluxes[i], "o",
                     mfc="none", label=str(camera_names[i]))
        plt.grid(True)
        ax = plt.gca()
        ax.set_title("Apparent Magnitude vs Height")
        ax.set_xlabel("Height (km)")
        ax.set_ylabel("Apparent Magnitude")
        ax.legend()
        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(ax.get_ylim()[::-1])

        # Saves the magnitude vs height figure
        plt.savefig(os.path.join(save_dir, str(meteor)+"_MAG.png"))


        # Temperature vs height figure
        fig = plt.figure(2)
        plt.plot(heights[:len(meteoroid_temp)], meteoroid_temp,
                 "k-", label="Best-Fit Model")
        plt.grid(True)
        ax = plt.gca()
        ax.set_title("Temperature vs Height")
        ax.set_xlabel("Height (km)")
        ax.set_ylabel("Temperature (K)")
        ax.legend()
        ax.set_xlim(ax.get_xlim()[::-1])

        # Saves the temperature vs height figure
        plt.savefig(os.path.join(save_dir, str(meteor)+"_TEMP.png"))


        # Velocity vs height figure
        fig = plt.figure(3)
        plt.plot(heights[:len(meteoroid_vel)], meteoroid_vel,
                 "k-", label="Best-Fit Model")
        plt.grid(True)
        ax = plt.gca()
        ax.set_title("Velocity vs Height")
        ax.set_xlabel("Height (km)")
        ax.set_ylabel("Velocity (m/s)")
        ax.legend()
        ax.set_xlim(ax.get_xlim()[::-1])

        # Saves the velocity vs height figure
        plt.savefig(os.path.join(save_dir, str(meteor)+"_VEL.png"))

        # -----------------------------------------------------------

        # HIGH ELEVATION MAGNITUDE
        fig = plt.figure(11)
        plt.plot(pre_times, pre_mag, "-k")
        plt.plot(10.6, 2, "r")
        ax = plt.gca()
        ax.set_title("High Up Time vs Magnitude")
        ax.set_ylim(ax.get_ylim()[::-1])
        #"""

        # ELEVATION THE METEOR IS SEEN
        #"""
        fig = plt.figure(12)
        ax = plt.gca()
        ax.set_title("Masses")
        ax.set_yscale('log')
        for i in range(len(holder_m)):
            plt.plot(list(range(0, (time.size-1), step_size))
                     [:len(holder_m[i])],
                     holder_m[i])

        fig = plt.figure(4)
        ax = plt.gca()
        ax.set_title("Velocities")
        for i in range(len(holder_m)):
            plt.plot(list(range(0, (time.size-1), step_size))
                     [:len(holder_m[i])],
                     holder_v[i])

        fig = plt.figure(5)
        ax = plt.gca()
        ax.set_title("Temperatures")
        for i in range(len(holder_m)):
            plt.plot(list(range(0, (time.size-1), step_size))
                     [:len(holder_m[i])],
                     holder_t[i])
            
        fig = plt.figure(6)
        ax = plt.gca()
        ax.set_title("Number of Particles per Bin")
        for i in range(len(holder_m)):
            plt.plot(list(range(0, (time.size-1), step_size))
                     [:len(holder_m[i])],
                     holder_n[i])
        #"""

        # HIGH ELEVATION
        #"""
        fig = plt.figure(7)
        ax = plt.gca()
        ax.set_title("High Up Masses")
        ax.set_yscale('log')
        for i in range(len(holder0_m)):
            plt.plot(pre_times[-len(holder0_m[i]):],
                     holder0_m[i])

        fig = plt.figure(8)
        ax = plt.gca()
        ax.set_title("High Up Velocities")
        for i in range(len(holder0_m)):
            plt.plot(pre_times[-len(holder0_m[i]):],
                     holder0_v[i])

        fig = plt.figure(9)
        ax = plt.gca()
        ax.set_title("High Up Temperatures")
        for i in range(len(holder0_m)):
            plt.plot(pre_times[-len(holder0_m[i]):],
                     holder0_t[i])
            
        fig = plt.figure(10)
        ax = plt.gca()
        ax.set_title("High Up Number of Particles per Bin")
        for i in range(len(holder0_m)):
            plt.plot(pre_times[-len(holder0_m[i]):],
                     holder0_n[i])
        
        # -----------------------------------------------------------
        # Note: Use set_yscale('log') to make the y-axis be plotted
        # on a log scale

        # Shows the plots on screen if need be
        plt.show()

    # ------------------------------------------------------------

    # Calculates the standard deviation of the calculated fluxes
    y_resid = np.subtract(np.array(dynamic_mag[0:end_mag]),
                          flux_smooth[0:end_mag])
    ss_resid = np.sum(np.square(y_resid))
    flux_s_dev = np.sqrt(ss_resid / (len(dynamic_mag) - 1))

    # This following block finds the standard deviation of the velocity
    end_vel = len(dynamic_vel)
    # Ignores the velocity = 0 and max_magni restrictions (?)
    if(ignore_mins):
        for i in range(end_vel, 1, -1):
            if(dynamic_vel[i] > 0):
                end_vel = i
                break

    # Note the division by 1000. It's done to convert km/s
    y_resid = np.subtract(np.array(dynamic_vel[0:end_vel]),
                          vel[0:end_vel]) / 1000
    ss_resid = np.sum(np.square(y_resid))
    vel_s_dev = np.sqrt(ss_resid / (len(dynamic_vel) - 1))

    # Returns the magnitude and velocity standard deviations
    return(flux_s_dev, vel_s_dev)
    
# ==========================================================================

# This entire section is composed of equations used throughout the script.
# Since they are so long, I decided that it'd be more organized if they
# were subfunctions down here and not just equations in the script.

def vel_change(dt, curr_atmo_dens, particle_vel,
               particle_mass, particle_density):
    """
    Gets the change in velocity over a small time increment dt.
    This equation describes the momentum transfer between the
    particle and air particles hitting it
    (Kikwaya, 2011, p.15)

    Parameter(s):
    dt               - The interval over which the change in mass
                       is calculated
    curr_atmo_dens   - The current atmospheric density
    particle_vel     - The velocity of a particle of a particle type
    particle_mass    - The mass of the each particle of a particle type
    particle_density - The density of a particle of a particle type

    Return:
    Returns the change in velocity as a float
    """
    
    dv = dt * drag_coeff * curr_atmo_dens * particle_vel**(2) \
         / particle_mass * shape_fact \
         * (particle_mass / particle_density)**(2/3)
    
    return dv


def sat_vapor_press(curr_height, heat_abl, atomic_m, boil_temp, particle_temp,
                    curr_atmo_dens):
    """
    Calculates the pressure of the saturated vapor from the particle
    on the particle. Used to calculate the changes in mass and
    temperature below

    Parameter(s):
    curr_height - The current height of the meteoroid in km
    heat_abl    - The heat of ablation of each particle of a particle type
    atomic_m    - The atomic mass of the particle/meteor vapor(amu)
    boil_temp   - The boiling temperature of each particle of a particle type

    Return:
    Returns the saturated vapor pressure on the meteor as a float
    """
    
    satur_vapor_press = pa_sea \
                        * np.exp(heat_abl * atomic_m / (kb * boil_temp)) \
                        * np.exp(-heat_abl * atomic_m / (kb * particle_temp))

    return satur_vapor_press


def small_mass(dt, curr_atmo_dens, heat_abl, particle_mass, particle_density,
               particle_vel, satur_vapor_press, particle_temp, atomic_m):
    """
    Calculates the loss in mass due to ablation over a small increment dt
    dm: Kikwaya, 2011, p.23
    dm2: https://www.spaceacademy.net.au/watch/debris/metflite.htm

    Parameter(s):
    dt                - The interval over which the change in mass
                        is calculated
    curr_atmo_dens    - The atmospheric density at the meteoroid's current
                        height
    heat_abl          - The heat of ablation of the meteoroid
    particle_mass     - The mass of the each particle of a particle type
    particle_density  - The density of a particle of a particle type
    particle_vel      - The velocity of the particle/grain
    satur_vapor_press - The saturated vapor pressure of the meteoroid's
                        vapor
    particle_temp     - The temperature of the meteoroid.
    atomic_m          - The atomic mass of the particle/meteor vapor(amu)

    Return:
    Returns the change in mass as a float
    """

    if particle_temp < 0:
        particle_temp = 0

    # Campbell-Brown (2004)
    dm = dt * shape_fact * \
         (particle_mass / particle_density)**(2/3) * sticking_fact \
         * (satur_vapor_press - vapor_press) / sqrt(2 * pi * kb \
                                                  * particle_temp \
                                                  / atomic_m)

    # Classical approach
    dm2 = dt * shape_fact * heat_trans_coeff / (2 * heat_abl) \
          * curr_atmo_dens * particle_vel**(3) \
          * (particle_mass / particle_density)**(2/3)

    return dm

def erosion_mass(dt, erosion_coeff, particle_mass, bulk_density,
                 particle_vel, curr_atmo_dens):
    """
    Calculates the mass lost due to erosion. Used in erode_mass.py to
    distribute that said mass into various mass bins in a Particles
    object. Equation from:
        Borovicka, J., Spurny, P., Koten, P., "Atmospheric
            deceleration and light curves of Draconid meteors
            and implications for the structures of cometary dust"
            (2007)

    Parameter(s):
    dt               - The time step over which the change in mass is
                       calculated 
    particle_mass    - The mass of the meteoroid material to be eroded away
    particle_density - The density of the material
    particle_vel     - The velocity of the material
    curr_atmo_dens   - The current atmospheric density at the height of the
                       main body

    Return:
    Returns the mass lost due to erosion within a time step dt.
    """

    erosion_mass = dt * erosion_coeff * drag_coeff * shape_fact \
                   * (particle_mass / bulk_density)**(2/3) \
                   * curr_atmo_dens * particle_vel**(3)

    if(erosion_mass > particle_mass):
        erosion_mass = particle_mass
    
    return erosion_mass

def small_temp(dt, spec_heat, particle_mass, curr_atmo_dens, particle_vel,
               particle_density, particle_temp, heat_abl, dm, atmo_temp):
    """
    Calculates the change in temperature over a small increment dt.
    (Kikwaya, 2011, p.22)

    Parameter(s):
    spec_heat         - The specific heat of the particle
    particle_mass     - The mass of the each particle of a particle type
    curr_atmo_dens    - The current atmospheric density
    particle_vel      - The velocity of a particle of a particle type
    particle_density  - The density of a particle of a particle type
    particle_temp     - The temperature of a ... type
    heat_abl          - The heat of ablation of a ... type
    dm                - The change in mass over a small interval dt of
                        a ... type

    Return:
    Returns a float of the change in temperature over dt
    """
    if(dm < 0):
        dm = 0

    if(particle_mass <= 0):
        return 0
    
    d_parti_temp = dt * 1 / (spec_heat * particle_mass) \
                   * (heat_loss_coeff * curr_atmo_dens \
                      * particle_vel**(3) / 2 * shape_fact \
                      * (particle_mass / particle_density)**(2/3) \
                      - 4 * stef_boltz * met_emissivity \
                      * (particle_temp**(4) - atmo_temp**(4)) \
                      * shape_fact \
                      * (particle_mass / particle_density)**(2/3) \
                      - 1 * heat_abl * dm / dt)

    return d_parti_temp
