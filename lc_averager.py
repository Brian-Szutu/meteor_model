# Imports numpy for mathematical arrays and math for the floor function
import numpy as np
from math import *

def averager(data_arr, cam_num):
    """
    Averages the passed in data array twice. Once over a certain interval
    and the second time over the each resulting "pair" from the first
    average.
    
    This script was originally built into Andy Crump's getFluxVTime34.m
    MATLab script. I decided that it'd be more organizationally sound
    to split this off from that script, as this serves an entirely
    different purpose compared to just getting flux and time.

    Parameter(s):
    data_arr - The array holding all of the data. It is actaully an array
               holding other arrays, with each inner array holding its own
               set of data.
    cam_num  - How many inner arrays there are aka sets of data. cam_nums
               is short for number of cameras, since in this case the
               script is being used to process data from different cameras
               on the same detected meteor.

    Return:
    The return is an array holding the processed data. The format is the
    same as the passed in data_arr.
    """

    # Initializes an empty list to hold the first averages
    data_run_1 = []

    # For each set of data in the data array...
    for i in range(cam_num):
        # Make a temporary list to store a processed set of data
        temp_array = []
        # For each element up to the halfway point in the set of data...
        for j in range(floor(len(data_arr[i]) / 2)):
            # Averages individual "pairs" (I think) together.
            # Then, it puts single average in the temporary list. Notice
            # the 2*j, which makes it so every odd indexed element is
            # paired with its own even indexed element
            temp_array.append((data_arr[i][2*j] + data_arr[i][2*j+1]) / 2)

        # Append the now processed set into the list used to hold all of
        # the processed sets
        data_run_1.append(temp_array)
        
    # Initializes an empty list to hold the second averages
    data_run_2 = []

    # For each set of data in the processed sets of data...
    for i in range(cam_num):
        # Make a temporary list to store a processed set of data
        temp_array = []
        # For each element in the set...
        for j in range(len(data_run_1[i]) - 1):
            # Average adjacent elements together to get another average.
            # Then, this single average is put into the temporary list
            temp_array.append((data_run_1[i][j+1] + data_run_1[i][j]) / 2)

        # Append the now processed set into the list used to hold all of the
        # processed sets
        data_run_2.append(np.array(temp_array))
    
    # Return a numpy array holding arrays of the processed data
    return data_run_2
