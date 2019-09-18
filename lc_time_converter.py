import numpy as np

def time_converter(time_strings):
    """
    Converts the times inside of a passed in list/numpy array from strings
    into times containing just the hour amount. These converted times are
    floats

    Parameter(s):
    time_strings - A list/array of times as strings in UT format.

    Return:
    Returns an array containing the converted UT times. The converted times
    are in seconds and are floats. The first time value matches the seconds
    place of the earliest time value.
    """

    # Initializes a list to hold the converted times
    time_floats = np.array([])

    # Holds the first time value in the passed in strings of time
    time_first_float = None
    seconds_first_float = None

    # For each time...
    for time in time_strings:
        # Retrieves the hours, minutes, and seconds of a specific time
        hours = time[0:2]
        minutes = time[3:5]
        seconds = time[6:]

        # Converts the hours, minutes, and seconds to floats and turns
        # the hours and minutes into seconds
        t_in_s = (float(hours) * 3600) + (float(minutes) * 60) + float(seconds)

        if(time_first_float == None):
            time_first_float = t_in_s
            seconds_first_float = float(seconds)

        # Appends the converted time into the array to be returned
        time_floats = np.append(time_floats,
                                t_in_s - time_first_float + seconds_first_float)

    # Returns the times as seconds
    return time_floats
