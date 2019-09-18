# Imports numpy for mathematical arrays and various related functions
import numpy as np
import matplotlib.pyplot as plt

# Imports various written scripts to make use of their functions
from lc_averager import averager
from lc_time_offsets import get_time_offsets
from lc_flux_offsets import get_flux_offsets
from lc_make_smooth import make_smooth

def flux_time_driver(nec_meteor_info, cam_num):
    """
    Runs through the detected fluxes by one or more cameras of a specific
    meteor and returns the smoothed times and fluxes along with the raw
    flux data. This script doesn't actually do most of the calculations and
    instead delegates that role to the several of the imported scripts.
    If there is more than one data set due to there being more
    than one camera, the times and fluxes are normalized to the
    data set with the longest light curve.

    Credit goes to Andy Crump for the MATLab scripts containing the
    algorithms.

    Parameter(s):
    nec_meteor_info - A numpy array containing the times of detection,
                      heights at those times, and visual magnitudes in
                      flux at those times of a specific meteor. Each
                      nested numpy array represents a different camera
                      that managed to detect that meteor. Name is short
                      for necessary meteor information.
    cam_num         - The number of cameras that managed to detect the
                      meteor. Used for indexing purposes in the imported
                      scripts.

    Return:
    Returns a numpy array with nested arrays containing the smoothed
    and raw flux data and associated times from ALL of the cameras. It is
    in the following format:
    [[Smoothed Time Data], [Raw Flux Data], [Smoothed Flux Data]]
    Also returns the detected magnitudes and heights in the form
    of arrays, with each nested array representing a different camera.
    """

    # Initializes lists to hold the times of detection, heights,
    # and fluxes of the meteor
    times = []
    heights = []
    fluxes = []

    # For each set, fill the each list with the appropriate information.
    # Everything is converted to a float in the process to be able to
    # do math with the values. In case the float conversion
    # in lightcurve_fit_driver is wonky
    for i in range(cam_num):
        times.append(np.array(nec_meteor_info[i][0], dtype=float))
        heights.append(np.array(nec_meteor_info[i][1], dtype=float))
        fluxes.append(np.array(nec_meteor_info[i][2], dtype=float))

    # Calls averager (name is kind of a misnomer) to cut down the number
    # of points to increase processing speed drastically
    avged_times = averager(times, cam_num)
    avged_heights = averager(heights, cam_num)
    avged_fluxes = averager(fluxes, cam_num)

    # For each data set, not including the first, get the time offsets
    # from the first data set's time
    for i in range(1, cam_num):
        time_offset = get_time_offsets(avged_times[0], avged_heights[0],
                                       avged_times[i], avged_heights[i])
        avged_times[i] = avged_times[i] - time_offset
        
    # Gets the flux offsets from the data set with the longest
    # light curve
    (avged_times, avged_fluxes, flux_offsets) = \
                  get_flux_offsets(avged_times, avged_fluxes,
                                   cam_num, False)

    # For each set of fluxes, account for the offset
    for i in range(cam_num):
        avged_fluxes[i] = np.array(avged_fluxes[i]) - flux_offsets[i]

    # Plots the detected fluxes -> calculated magnitudes of the
    # meteoroid 
    #fig = plt.figure(1)
    #plt.ion()
    """
    for i in range(cam_num):
        plt.plot(avged_times[i], avged_fluxes[i], "o", mfc="none",
                 label=str(camera_names[i]))
    fig.canvas.draw()
    """

    # Initializes empty numpy arrays to store ALL of the times and
    # fluxes without using nested arrays
    time_raw = np.array([])
    flux_raw = np.array([])

    # For each data set, put the appropriate data in the appropriate array
    for i in range(cam_num):
        time_raw = np.concatenate((time_raw, avged_times[i]))
        flux_raw = np.concatenate((flux_raw, avged_fluxes[i]))

    # Gives out the indices that would be given if the times were to be
    # sorted from least to greatest
    sorted_ind = np.argsort(time_raw)
    # Sorts the times from lowest to greatest hours
    time_raw = np.sort(time_raw)
    # Moves the fluxes associated with the same times around so that
    # both have the same indices
    flux_raw = flux_raw[sorted_ind]

    # Smooths the times and fluxes
    (time_smooth, flux_smooth) = make_smooth(time_raw, flux_raw, 5)

    # Plots the smoothed magnitudes on top of the detected
    # magnitudes
    """
    plt.plot(time_smooth, flux_smooth, "c-")
    ax = plt.gca()
    ax.set_title("Time vs Magnitude")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Magnitude")
    ax.set_ylim(ax.get_ylim()[::-1])
    """
    
    # Returns the smoothed times and fluxes and the raw fluxes.
    return (time_smooth, flux_raw, flux_smooth, avged_fluxes, avged_heights)
        
