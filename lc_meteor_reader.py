# Imports the os module for path manipulation and
# numpy for mathematical arrays
import os
import numpy as np

def meteor_reader(meteor_dir):
    """
    Accesses the text file of a specific detected meteor from a specific
    camera. Then, the meteor's various times at which it was detected,
    the height of the meteor at those times, and its visual magnitude in
    flux at those times are read from the file and returned. The latitudes
    and longitudes are also returned if the velocity model in the text file
    fails.

    Parameter(s):
    meteor_dir - The FULL path to the text file containing the specific
                 meteor.

    Return:
    A collection of arrays containing the specific meteor's range
    of times of detection, the meteor's height corresponding to those
    times, the velocities at those times and its visual magnitude in
    flux at the recorded times. Also returns the recorded latitudes and
    longitudes of the meteor
    """

    # Initializes an empty list to hold the meteor info in
    meteor = []

    # Opens the text file of the specific meteor and saves
    # each line as its own element in the meteor list
    with open(meteor_dir, "r") as met_read:
        meteor = [line.split() for line in met_read]

    # Closes the meteor text file
    met_read.close()

    # Lops off the header information from the read file
    meteor = meteor[8:]
    # Converts the list into a numpy array
    meteor = np.array(meteor)

    # Returns (Times of Detection, Heights,
    #          Velocities, Visual Magnitudes in Flux, Latitudes,
    #          Longitudes)
    return (meteor[:, 1], meteor[:, 3], meteor[:, 4], meteor[:, 7],
            meteor[:, 8], meteor[:, 9])
