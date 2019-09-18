# Imports numpy for mathematical arrays. math is for several math
# functions and operations like cos() and sqrt()
# particles is used to simulate a system of particles
import numpy as np
from math import *
from lc_constants import *
from lc_particles import Particles

def breakup_model(particle_sys, p_len, p, constants_change, zenith):
    """
    This script models a meteor breaking up into smaller particles. It does
    the "breaking up" of particles and returns another Particles object
    that came from the passed in Particles object.
    For documentation on the Particles class, see particles.py.

    Credit goes to Andy Crump for his breakup30.m script.

    Parameters:
    particle_sys     - A Particles object used to simulate a system of
                       particles with varying masses, particle counts,
                       velocities, surface temperatures, and densities
    p_len            - The number of mass bins in Particles
    p                - The specific mass bin. It's really an index
                       used to access the specific particle type
    constants_change - A dictionary holding various values that change
                       due to being in loop(s). These values change to
                       find a model that most closely resembles a particular
                       meteor
    zenith           - The zenith angle of a particular meteor

    Return:
    Returns the passed in Particles object, but changed to model a
    breakup of particle(s). Also returns the number of different
    "types" of particles in the system.
    """

    # Note to self: s = -2.2735*Fvalue + 2.907,
    # LightcurvefitDriver.m, line 103

    # A bunch of constants from the passed in constants dictionary
    # held by initialized variables

    # The grain density of the particle
    grain_dens = constants_change["grain_dens"]
    # The specific heat of the particle
    spec_heat = constants_change["spec_heat"]
    # The thermal conductivity of the particle
    thermal_cond = constants_change["thermal_cond"]
    # The mass distribution of the particles
    s = constants_change["s"]

    # If the mass of the first particle is less than 0 kg, don't do
    # the breakup algorithm
    if(particle_sys.get_value("MASS", 0) <= 0):
        particle_sys.change_values(0, 0, "SAME", "SAME", "SAME", "SAME")
        return (p_len, particle_sys)

# =======================================================================
    
    # Heating depth equation variable initializations

    # Gets the particle's density (=bulk_dens for the main particle aka 0)
    particle_density = particle_sys.get_value("DENSITY", 0)
    
    # The thickness of the meteoroid shell being heated
    x_0 = sqrt(thermal_cond / (particle_density * spec_heat)) * \
          sqrt(h_atm * 1000/ (particle_sys.get_value("VELOCITY", p) * \
                        cos(radians(zenith))))
    
    # Radius of the meteoroid
    r_total = (3/4 * particle_sys.get_value("MASS", 0) / \
               (constants_change["bulk_dens"] * pi))**(1/3)

    # Prevents r_core from being negative
    if(x_0 > r_total):
        x_0 = r_total

    # Radius of the meteoroid "core"
    r_core = r_total - x_0

    # Total mass aka the entire mass of the single, first particle
    m_total = particle_sys.get_value("MASS", 0)

    # Mass of the core. Since only the original particle aka the
    # main body of the meteroid breaks up, the density is the bulk density
    m_core = (4/3) * pi * particle_density * r_core**(3)

    # The total mass of the fragments in the outer shell
    m_shell = m_total - m_core

    # The surface temperature of the meteor
    t_surf = particle_sys.get_value("TEMPERATURE", 0)

    # The core temperature as a function of the surface temperature
    t_core = 1 / exp(1) * (t_surf - initial_temp) + initial_temp

    # The largest possible mass allowed to break from the shell. Based
    # off of the grain density.
    max_mass = (4/3) * pi * grain_dens * (x_0 / 2)**(3)

    # If the maximum allowed mass is smaller than the minimum
    # breakup mass, don't do the breakup algorithm
    if(max_mass < min_breakup_mass):
        # Dummy return values
        return (p_len, particle_sys)

# =======================================================================

    # Adjust the mass and temperature of the original particle
    # to the core values due to breakup
    particle_sys.change_values(0, m_core, "SAME", "SAME", t_core, "SAME")

    # For each different mass between the minimum breakup mass
    # and the maximum breakup mass threshold, add new types of
    # particles to the Particles object
    mass_range = np.arange(log10(max_mass),
                           (log10(min_breakup_mass) + (-d_log_mass)),
                           (-d_log_mass))

    # If the ending mass in the mass bin range is less than the
    # allowed minimum mass, remove it.
    if(mass_range[-1] < log10(min_breakup_mass)):
        mass_range = mass_range[:-1]

    for log_m in mass_range:
        # Adds a new mass, particle count, velocity, temperature,
        # and density (changes done below in that order)
        # Note the scaling in the particle count calculation according to s
        particle_sys.new_particle(10**log_m,
                                  (1/(s-1) * ((10**(log_m
                                                    - d_log_mass ))**(-s+1)
                                              - (10**(log_m))**(-s+1))),
                                  particle_sys.get_value("VELOCITY", 0),
                                  t_surf,
                                  grain_dens)

    # ------------------------------------
    # This block makes sure the total mass of the new particles equals the
    # shell mass

    # The new total number of mass bins in the system
    new_len = particle_sys.get_bin_num()
    
    # Sums up the new particles' masses to rescale the shell mass
    m_sum = 0
    for j in np.arange(p_len, new_len):
        m_sum += particle_sys.get_value("PARTICLE_N", j) * \
              particle_sys.get_value("MASS", j)

    # Rescales the number of particles in each mass bin
    for j in np.arange(p_len, new_len):
        new_n = particle_sys.get_value("PARTICLE_N", j) * m_shell / m_sum
        particle_sys.change_values(j, "SAME", new_n, "SAME", "SAME", "SAME")

    # ------------------------------------

    # Sets the value to be returned to the counted number of mass bins
    p_len = new_len

    # Returns the number of mass bins and the broken up Particles system 
    return (p_len, particle_sys)
