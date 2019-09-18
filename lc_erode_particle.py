# Imports numpy to use the log10 function
import numpy as np

# Imports relevant constants from constants.py and the Particles class
from lc_constants import *
from lc_particles import Particles

def erode_particle(particle_sys, constants_change, mass):
    """
    Creates mass bins from the particles to be eroded from the meteoroid
    aka the particle in the 0th mass bin.
    This has a very similar algorithm to breakup_model.py.

    Parameter(s):
    particle_sys     - A Particles object. Holds the mass bins of the particles
                       of the meteoroid. See particles.py for more information.
    constants_change - A dictionary containing various constants and the
                       free parameters of the meteoroid.
    mass             - The amount of mass to be eroded away.

    Return:
    Returns the updated Particles object
    """

    # If the mass to be eroded away is less than the ignore mass,
    # don't erode.
    if(mass <= 0):
        # Dummy return variable
        return particle_sys

    p_len = particle_sys.get_bin_num()

    # The mass's temperature. Assumed to be the same temperature as the
    # main mass's temperature
    temperature = particle_sys.get_value("TEMPERATURE", 0)

    # Creates a range of masses to iterate over to determine the new mass
    # bins to be added to the system
    mass_range = np.arange(log10(mass),
                           (log10(min_breakup_mass) + (-d_log_mass)),
                           (-d_log_mass))

    if(len(mass_range)==0):
        mass_range = [min_breakup_mass]
    else:
        # If the smallest mass of the range of values is lower than the lowest
        # allowed mass a bin can have, remove it
        if(mass_range[-1] < log10(min_breakup_mass)):
            mass_range = mass_range[:-1]

    # Gets the grain densitys. It is a free parameter
    grain_dens = constants_change["grain_dens"]

    # Gets the mass distribution of the meteoroid. Also a free parameter
    s = constants_change["s"]

    # For each mass in the generated range, create a new mass bin.
    for log_m in mass_range:
        particle_sys.new_particle(10**log_m,
                                  (1/(s-1) * ((10**(log_m - \
                                                     d_log_mass ))**(-s+1) -
                                                (10**(log_m))**(-s+1))),
                                   particle_sys.get_value("VELOCITY", 0),
                                   temperature,
                                   grain_dens)

    #print(particle_sys.get_value("VELOCITY", 0))

    # Gets the new number of mass bins in the system
    new_len = particle_sys.get_bin_num()

    # ============================================================
    # The following section scales the total masses of the mass bins.
    # The total mass will equal the mass eroded from the meteoroid
    m_sum = 0
    for j in np.arange(p_len, new_len):
        m_sum += particle_sys.get_value("PARTICLE_N", j) * \
              particle_sys.get_value("MASS", j)

    for j in np.arange(p_len, new_len):
        new_n = particle_sys.get_value("PARTICLE_N", j) * mass / m_sum
        """
        if(new_n < 1):
            new_n = 0
        """
        particle_sys.change_values(j, "SAME", new_n, "SAME", "SAME", "SAME")


    # ============================================================
    # Returns the Particles object with the eroded mass inside of the mass
    # bins
    return particle_sys
