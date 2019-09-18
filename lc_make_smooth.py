# Imports numpy for mathematical arrays
import numpy as np

def make_smooth(x, y, interval):
    """
    Smooths the inputted y values using the concept of moving averages.
    The larger the interval, the more "smooth" the smoothing will be. The
    smaller the interval, the closer the smoothing will be to the actual
    points.
    The interval term is defined like so:
    For an interval n, the moving average is defined by the n-1
    previous points and the current point.

    Credit goes to Andy Crump for laying out the algorithm
    in his MATLab script.

    Parameter(s):
    x        - The corresponding x values to the y values. A numpy array
    y        - The values to be smoothed. A numpy array
    interval - An integer value that determines how "smooth" the
               smoothing will be, as described in the description above.

    Return:
    The passed in x values and the smoothed y values
    """

    # If the interval is larger than the number of points to be smoothed,
    # set the interval size to be the number of points
    if(interval > x.shape[0]):
        interval = x.shape[0]

    # If the interval size is even, make it odd.
    if(interval%2 == 0):
        interval -= 1

    # Used to determine how many points to include in the moving
    # average implementation in the loop below
    dist = (interval-1)/2

    # Initialize an empty list to hold the averages
    avg = []

    # For each point in y, compute a median for a specific range of
    # points and make that the moving average...
    for i in range(y.shape[0]):
        if(i <= dist):
            # Finds the median of the points from the first to the
            # (i+1+dist)th point. The +1 compared is added since Python
            # ignores the final term in array indexing
            avg.append(np.median(y[0:int(i+dist)]))
        elif(i > x.shape[0]-dist):
            # Finds the median of the points from the (i-dist)th point
            # to the last point
            avg.append(np.median(y[int(i-dist):(x.shape[0])]))
        else:
            # Finds the median of the points from the (i-dist)th point
            # to the (i+1+dist)th point
            avg.append(np.median(y[int(i-dist):int(i+dist)]))

    # Return the passed in x values and the moving averages
    return (x, np.array(avg))
