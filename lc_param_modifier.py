# Imports numpy to make use of argsort and sort
import numpy as np

# These other functions are needed to take special averages, find
# the offsets between the each set of camera times, and
# smooth out the data.
from lc_averager import averager
from lc_time_offsets import get_time_offsets
from lc_make_smooth import make_smooth

def param_modifier(times, velocities, heights, cam_num):
    """
    Similar to flux_time_driver. However, this script uses that process
    to simply make it so the number of elements in the lists to be used
    in the calculations are the same length

    Credit goes to Andy Crump for the algorithm used in flux_time_driver,
    which is also used in here.

    Parameter(s):
    times      - The time strings in seconds. May hold nested lists
                 if there are multiple cameras
    velocities - The velocities of the meteor in km/s. May hold nested lists
    heights    - The heights of the meteor. May hold nested lists
    cam_num    - The number of cameras that detected the meteor

    Return:
    Returns the modified times, velocities, and heights
    """
    
    # Gets the special averages of the passed in data
    avg_times = averager(times, cam_num)
    avg_vel = averager(velocities, cam_num)
    avg_heights = averager(heights, cam_num)

    # Zeros the times from the different cameras to the "first" camera
    # in the sets of times
    for i in range(1, cam_num):
        time_offset = get_time_offsets(avg_times[0], avg_heights[0],
                                       avg_times[i], avg_heights[i])
        avg_times[i] = avg_times[i] - time_offset
        
    # Puts all of the cameras' data into appropriate arrays. There
    # are no nested lists here
    
    # Times array
    all_times = np.array([])
    for i in range(cam_num):
        all_times = np.concatenate((all_times, avg_times[i]), axis=0)

    # Velocities array
    all_vel = np.array([])
    for i in range(cam_num):
        all_vel = np.concatenate((all_vel, avg_vel[i]), axis=0)

    # Heights array
    all_heights = np.array([])
    for i in range(cam_num):
        all_heights = np.concatenate((all_heights, avg_heights[i]), axis=0)

    # Gets the indices of the sorted times list
    sorted_ind = np.argsort(all_times)

    # Sorts the velocities and heights so they will stay paired up with
    # the times
    all_vel = all_vel[sorted_ind]
    all_heights = all_heights[sorted_ind]
    all_times = np.sort(all_times)

    # Smooths the times and velocities, since they are needed in calculations
    (all_times, all_vel) = make_smooth(all_times, all_vel, 15)

    # Returns the times, velocities, and heights
    return(all_times, all_vel, all_heights)
