# Imports numpy for numpy arrays
import numpy as np

class Particles:
    """
    Simulates a system of particle mass bins to be used in a breakup model.
    This is used in place of the particle structs in Andy Crump's MATLab
    script. Each index of an array/list represents one particle.
    So for example, an index of 0 represents the 0th particle type, an
    index of 1 represents the 1st particle type, and so on.

    Attributes:
    mass         - Contains the mass of each mass bin
    particle_num - Contains the number of particles in each bin
    velocity     - Contains the velocities of each bin
    temperature  - Contains the temperatures of each bin
    density      - Contains the densities of each bin

    Methods:
    new_particle    - Adds a new mass bin to the system
    remove_particle - Removes a mass bin
    get_particle    - Gets a specific bin from the system
    get_value       - Gets a specific parameter of a specific bin
                      indicated by an index
    change_values   - Changes certain parameter(s) of a mass bin.
                      If the string SAME is passed in as a new parameter,
                      the corresponding particle parameter will not change
    get_bin_num     - Gets the number of bins in the system
    reset           - Resets the Particle system to have no particless
    """

    def __init__(self):
        """
        The constructor for the class 
        """
        # Empty lists to hold all inputted values. Temperature in this
        # instance refers to the surface temperature.
        self.mass = []
        self.particle_num = []
        self.velocity = []
        self.temperature = []
        self.density = []


    def new_particle(self, new_m, new_n, new_v,
                     new_t, new_rho):
        """
        Creates a mass bin with the inputted parameters

        Parameter(s):
        new_m   - The new particle type's mass
        new_n   - The new particle type's particle count
        new_v   - The new particle type's velocity
        new_t   - The new particle type's temperature
        new_rho - The new particle type's density
        """

        # Adds the new particle bin's values to the lists
        self.mass.append(new_m)
        self.particle_num.append(new_n)
        self.velocity.append(new_v)
        self.temperature.append(new_t)
        self.density.append(new_rho)

    def remove_nones(self):
        """
        Removes a mass bin from the system. Used in breakup_model and
        get_mass to remove bins that are completely gone. 

        Parameter(s):
        index - An integer. Used remove the specific bin, from
                the 0th bin to the nth bin.

        Return:
        None
        """

        # The filter() function is being used to remove the bins with
        # None as their data values and only None. The .__ne__ at the
        # end of None tells Python to ignore 0, since it is included
        # in what is counted as None. Note that this method only works
        # for Python 3 and above.
        self.mass = list(filter(None.__ne__, self.mass))
        self.particle_num = list(filter(None.__ne__, self.particle_num))
        self.velocity = list(filter(None.__ne__, self.velocity))
        self.temperature = list(filter(None.__ne__, self.temperature))
        self.density = list(filter(None.__ne__, self.density))
        

    def get_particle(self, index):
        """
        Gets the parameters of a certain bin

        Parameter(s):
        index - The index of the particle type. It's determined by
                when the particle type was added in terms of order

        Return:
        A list containing the parameters of the particle type. The
        format is as follows:
        [mass, count, velocity, temperature, density]
        """
        return [self.mass[index], self.particle_num[index],
                self.velocity[index], self.temperature[index],
                self.density[index]]


    def get_value(self, name, index):
        """
        Gets a specific paramater of a specific mass bin

        Parameter(s):
        name  - The name of the parameter as a string
        index - The index of the bin

        Return:
        A float containing the requested parameter
        """
        if(name == "MASS"):
            return self.mass[index]
        if(name == "PARTICLE_N"):
            return self.particle_num[index]
        if(name == "VELOCITY"):
            return self.velocity[index]
        if(name == "TEMPERATURE"):
            return self.temperature[index]
        if(name == "DENSITY"):
            return self.density[index]

    def get_bin_num(self):
        """
        Gets the number of bins within the particle system.

        Parameter(s):
        None

        Return:
        Returns the number of mass bins as an integer
        """

        # Just returns the length of list keeping track of all
        # of the masses, since all of the lists should be the
        # same length
        return len(self.mass)


    def change_values(self, index, new_m, new_n, new_v,
                      new_t, new_rho):
        """
        Changes the parameters of an bin.
        If the string SAME is inputted as one of the new parameters,
        that parameter will stay the same.

        Parameter(s):
        new_m   - The particle type's new mass
        new_n   - The particle type's new particle count
        new_v   - The particle type's new velocity
        new_t   - The particle type's new temperature
        new_rho - The particle type's new density

        Return:
        None
        """
        if new_m != "SAME":
            self.mass[index] = new_m
        if new_n != "SAME":
            self.particle_num[index] = new_n
        if new_v != "SAME":
            self.velocity[index] = new_v
        if new_t != "SAME":
            self.temperature[index] = new_t
        if new_rho != "SAME":
            self.density[index] = new_rho


    def reset(self):
        """
        Resets the particle system to have no bins

        Parameter(s):
        None

        Return:
        None
        """

        # Reinitalizes the instance variables to be empty
        # lists
        self.mass = []
        self.particle_num = []
        self.velocity = []
        self.temperature = []
        self.density = []
