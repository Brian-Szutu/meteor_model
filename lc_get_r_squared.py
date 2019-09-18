# Imports numpy to use the statistics-related
# functions in numpy
import numpy as np

def get_r_squared(x, y, line):
    """
    Calculates the coefficient of determination, or r-squared, between
    the line generated by the model and the data points composed of x and
    y. This value indicates how close the model is to the measured points.

    Thanks to Andy Crump for the formulae and expressions for calculating
    the r-squared value.

    Parameter(s):
    x    - The x values of the measured points
    y    - The y values of the measured points
    line - The model-generated line

    Return:
    Returns the r-squared value as a single value
    """

    # Produces various values from the model to be compared to
    y_fit = np.polyval(line, x)

    # Finds the differences between the measured and theoretical values
    y_resid = np.subtract(y, y_fit)

    # Finds the required values to find r-squared
    ss_resid = np.sum(np.square(y_resid))
    ss_total = (y.shape[0] - 1) * np.var(y)

    # Calculates r-squared
    rsq = 1 - ss_resid / ss_total

    # Returns r-squared
    return rsq
