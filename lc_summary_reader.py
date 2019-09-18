# Imports the os module for path manipulation and
# numpy for mathematical arrays
import os
import numpy as np

def summary_reader(summary_dir):
    """
    Accesses the SummaryMeteorLog.txt within a folder and finds
    the meteor numbers, dates of detection, times of detection,
    and the zenith angles of the meteors detected on a given
    day.

    Parameter(s):
    summary_dir - The FULL path to SummaryMeteorLog.txt file.

    Return:
    A collection of numpy arrays containing the meteor numbers,
    dates of detection, times of detection, and the zenith angles
    of the meteors detected on a given day. The last array contains
    the f-values of the light curves. The order returned is the order
    mentioned.
    """

    # Initializes an empty list to hold the meteor summary info in
    summary = []

    # Opens the SummaryMeteorLog.txt and saves each line as its own
    # element in the summary list
    with open(summary_dir, "r") as sum_read:
        summary = [line.split() for line in sum_read]

    # Closes the opened SummaryMeteorLog.txt file
    sum_read.close()

    # Lops off the header information that was read from the file
    summary = summary[3:]
    # Converts the list into a numpy array
    summary = np.array(summary)

    # Returns [Meteor Number, Date of Detection,
    #          Time of Detection, Zenith Angle, f-value]
    return (summary[:, 0], summary[:, 1], summary[:, 2],
            summary[:, 29], summary[:, 33])

