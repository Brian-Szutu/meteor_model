# Imports numpy to use the polyval function
import numpy as np

def remove_outliers(x, y, index, curve):
    """
    This script removes outliers from a data set. If a y-value
    is over twice the average vertical distance, or average delta,
    away from the curve, it will be removed.

    Credit goes to Andy Crump for the algorithm.

    Parameter(s):
    x     - The x values of the data. May hold more than one set of data.
            Nested numpy arrays
    y     - The y values of the data. May hold more than one set of data
            Nested numpy arrays
    index - Indicates which set of data to remove the outliers from.
    curve - The polyfit curve the data is to be compared to.

    Return:
    A numpy array holding the data set specified by the index with
    the outliers removed
    """

    # Initializes two variables to use to calculate the average delta
    # between the data points and the curve is.
    abs_delta_sum = 0
    deltas_taken = 0

    # For each y value in the specific set, calculate the y value delta
    for i in range(x[index].size):
        abs_delta = abs(y[index][i]) - np.polyval(curve, x[index][i])
        abs_delta_sum += abs_delta
        deltas_taken += 1

    # Find the average delta
    avg_delta = abs_delta_sum / deltas_taken

    # Initialize two lists to hold the data set without the outliers
    x_out = []
    y_out = []

    # For each element in the set...
    for i in range(x[index].size):

        # Recalculate the delta
        abs_delta = abs(y[index][i]) - np.polyval(curve, x[index][i])

        # And if it is within two average deltas of the curve, add the
        # data point to the output
        if(abs_delta <= (2 * avg_delta)):
            x_out.append(x[index][i])
            y_out.append(y[index][i])

    # Return the passed in data sets with the outliers of a specific data
    # set removed
    return (np.array(x_out), np.array(y_out))
