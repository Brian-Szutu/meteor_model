# Imports numpy for mathematical arrays. get_r_squared is used to
# get the r-squared value of the longest/most
# comprehensive light curve. remove_outliers removes the outliers
# in a data set
import numpy as np
from lc_get_r_squared import get_r_squared
from lc_remove_outliers import remove_outliers

def get_flux_offsets(time_in, flux_in, cam_nums, no_outliers):
    """
    Calculates the offsets between the longest curve the other curves,
    if applicable. Will also remove outliers if requested.

    Credit goes to Andy Crump for the algorithms

    Parameter(s):
    time_in     - A list with nested numpy arrays. Contains all of
                  the detection times from all cameras for a single
                  meteor. Each nested array represents a single
                  camera's data
    flux_in     - A list with nested numpy arrays. Contains all of
                  the detected fluxes from all cameras for a
                  single meteor. Each nested array represents a
                  single camera's data
    cam_nums    - The number of cameras in that detected the
                  meteor.
    no_outliers - A boolean used to determine whether or not to
                  remove outliers in the data sets. True = remove,
                  False = don't touch

    Return:
    Returns three numpy arrays. The first array contains the
    processed time data. The second contains the processed flux data.
    The third contains the flux offset values.
    """

    # Initializes the return arrays to be equal to the passed in ones.
    # This is done so that if no outlier removal is required, the original
    # arrays are passed back out
    time_out = time_in
    flux_out = flux_in

    # The index of the longest light curve. The index refers to the light
    # curve's data set, as the inputted arrays contain multiple
    # data sets
    longest_ind = 0

    # For each data set, find out which one has the most elements
    for i in range(1, cam_nums):
        if(time_out[i].shape[0] > time_out[longest_ind].shape[0]):
            longest_ind = i

    # Create a quadratic regression curve for the longest light curve
    longest_curve = np.polyfit(time_out[longest_ind], flux_out[longest_ind],
                               2)

    # Find the r-squared value between the regression curve and
    # the longest light curve
    rsq = get_r_squared(time_out[longest_ind], flux_out[longest_ind],
                        longest_curve)

# ======================================================================

    # If outliers are to be removed, remove them
    if(no_outliers):
        # Remove the outliers in the longest light curve
        (time_out[longest_ind], flux_out[longest_ind]) = \
                                remove_outliers(time_out, flux_out,
                                                longest_ind, longest_curve)
        # Create a new quadratic regression curve for the no outlier
        # longest light curve
        longest_curve = np.polyfit(time_out[longest_ind],
                                   flux_out[longest_ind], 2)

        # Get a new r-squared value between the no outlier light
        # curve and its regression curve
        rsq = get_r_squared(time_out[longest_ind], flux_out[longest_ind],
                            longest_curve)

        # For each light curve that ISN'T the longest light curve
        for i in range(cam_nums):
            if(i != longest_ind):

                # Remove the outliers in them
                (time_out[longest_ind],
                 flux_out[longest_ind]) = \
                 remove_outliers(time_out[i], flux_out[i],
                                 i, longest_curve)

# ======================================================================

    # Initializes a list to contain the average vertical distance, or
    # average deltas, between the longest quadratic regression curve
    # and each data set.
    avg_deltas = [0] * cam_nums

    # For each data set, excluding the longest light curve (since it'd
    # be compared to itself if this weren't the case)
    for i in range(cam_nums):
        if(i != longest_ind):

            # Calculate the average delta between the set and the
            # regression curve
            delta_sum = 0
            deltas_taken = 0

            for j in range(time_out[i].shape[0]):
                delta = flux_out[i][j] - np.polyval(longest_curve,
                                                    time_out[i][j])
                delta_sum += delta
                deltas_taken += 1

            # Store the average delta of the set in a corresponding spot
            # in the list
            avg_deltas[i] = delta_sum / deltas_taken

    # Finds the average of the average deltas
    avg_of_avg_delt = sum(avg_deltas) / cam_nums

    # Finds the offset of each data set by subtracting their corresponding
    # average delta by the average of the average deltas.
    offsets = []
    for element in avg_deltas:
        offsets.append(element - avg_of_avg_delt)

    # Returns the (maybe) processed time and flux arrays in addition to the
    # offset values
    return (np.array(time_out), np.array(flux_out), np.array(offsets))
