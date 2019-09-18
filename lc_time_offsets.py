# Imports numpy for the polyfit function. get_r_squared is used
# to get the r-squared values between a line and data points
import numpy as np
from lc_get_r_squared import get_r_squared

def get_time_offsets(time_init, height_init, time_off, height_off):
    """
    Calculates the offset between a time and height set deemed to be the
    starting values and another time and height set. Height is used as the
    y-axis since it has a somewhat linear relationship with time at
    the speeds meteors fall through the atmosphere

    Credit goes to Andy Crump's MATLab script.

    Parameter(s):
    time_init   - Time values considered to be the initial starting point
    height_init - Height values considered to be the initial starting point
    time_off    - Time values to be corrected with the output offset
    height_off  - Height values to be corrected with the output offset

    Return:
    Returns the time offset, as a float, between the heights and
    times considered as the starting point and another set of heights
    and times
    """

    # Produces a two element array as so: [m, b]
    # m is the slope of the line and b is the y intercept (y = mx + b)
    init_line = np.polyfit(time_init, height_init, 1)

    # Produces a best-fit line for the non-initial heights and times
    off_line = np.polyfit(time_off, height_off, 1)

    # If the non-initial time has less than two times, make the first
    # height and time the ones used to calculate the offset
    if(time_off.shape[0] < 2):
        t2 = time_off[0]
        h2 = height_off[0]

    # If there are more than two or more times, take the average between
    # the first two heights and times and use those in the offset
    # calculation
    else:
        t2 = (time_off[0] + time_off[1]) / 2
        h2 = (height_off[0] + height_off[1]) / 2

    # Slope of the initial best fit line
    m = init_line[0]
    # y-intercept of the initial best fit line
    b = init_line[1]

    # Calculates the offset between the initial data and the non-initial
    # data
    offset = t2 - (h2 - b) / m

    # Returns the offset
    return offset
