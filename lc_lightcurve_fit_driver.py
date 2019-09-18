# os is used for path manipulation. Numpy is used to create
# and manipulate mathematical arrays.
# math is used for functions like ceil(). statistics is used
# for standard deviation and mean calculations
import os
import numpy as np
from math import *
from statistics import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Imports various constants and values from the written
# constants.py file
from lc_constants import *

# Imports various other functions to be used in things
# like file searching/reading and model fitting
from lc_file_finder import file_finder
from lc_summary_reader import summary_reader
from lc_meteor_reader import meteor_reader
from lc_time_converter import time_converter
from lc_param_modifier import param_modifier
from lc_flux_time_driver import flux_time_driver
from lc_get_mass import get_mass

def main():
    # Sets up a way to use the initialized lists in constants.py
    global density_array, heat_abl_array, glue_array, boil_array, \
           thermal_array, atomic_m_array, spec_heat_array, \
           log_dens_coeff, log_low_t_coeff, log_high_t_coeff
    
    """
    The main function. Goes through every folder in the specified working
    directory. It goes through every single individual meteor in each folder
    and finds various theoretical parameters of each meteor, ranging from
    its molar mass to its density, and saves those parameters to a .txt file
    that will hold every single meteor's parameters.

    flux_time_driver and get_mass are the two largest computational functions
    this script calls.

    Credit goes to Andy Crump for his LightcurveFitDriver.m MATLab script.

    Parameter(s):
    None

    Return:
    None
    """

    # Used to get the working directory of the script
    #working_dir = input("remove later: ")
    working_dir = "C:\\Users\Brian\\Documents\\REU\\Script_Translation\\Layer_1"

    # For each folder in the working directory...
    for folder in os.listdir(working_dir):

        # Save the folder directory for convenience
        folder_dir = os.path.join(working_dir, folder)

        # If the thing selected is NOT a folder, move on to the next
        # object in the working directory
        if os.path.isfile(folder_dir):
            continue

        # Find each camera .txt for each meteor and a summary
        # .txt that should be present in the folder
        (image_txts, summary_txt) = file_finder(folder_dir)

        # If there is no summary .txt, move on to the next folder
        if summary_txt == "":
            print("SummaryMeteorLog.txt is not in " + folder + "!")
            print("Moving on to the next folder...")
            continue

        summary_dir = os.path.join(folder_dir, summary_txt)

        # Gets the relevant data from the summary file
        (meteor_nums, dates, time_strings, zeniths, f_vals) = \
                      summary_reader(os.path.join(folder_dir, summary_dir))

        # Converts the zenith angles to floats
        zeniths = list(np.float_(zeniths))

# ============================================================================

        # This block writes the header of the output file

        # Keeps track of whether or not the header has been written.
        # Used as a check to make sure new information isn't appended
        # into a file with the same but older information
        header_tracker = False

        # If a fresh new output file is desired or there isn't an
        # output file already existing...
        if(b_start_new_file or not\
               os.path.isfile(os.path.join(folder_dir, "FitSummary.txt"))):

            # Create the output file
            header_writer(folder_dir)
            # Set the boolean keeping track of whether or not to create
            # a new output file to False
            b_start_file_over = False
            # Set the boolean tracking the header creation to True
            header_tracker = True

        if(not os.path.isdir(os.path.join(working_dir, "LightcurvePlots"))):
           os.mkdir(os.path.join(working_dir, "LightcurvePlots"))

# ============================================================================

        # Creates a counter variable to keep track of the meteor. This
        # number is used as an index to access the information read from
        # the summary file
        meteor_index = 0

        # For each meteor number in the summary file...
        for meteor in meteor_nums:
            
            # A variable referenced below. No idea what it does
            b_all_zeros = False

            # Create a copy of the constants/value dictionary from
            # constants.py
            constants_change = K

            # Initialize a list to hold the camera .txt files pertaining
            # to the specific meteor, a variable to keep track of how
            # many of those files exist, and whether or not the specific
            # meteor was found in the folder
            meteor_files = []
            meteor_file_counter = 0
            found_meteor = False

            # Keeps track of the cameras that observed the meteor
            camera_names = []

            # Sorts the camara .txts so each meteor has its own "block"
            image_txts.sort()

            # For each camera .txt found in the folder...
            for file in image_txts:
                # If the meteor number is in the file name...
                if meteor in file:
                    # Append that file name to the list above and add
                    # 1 to the counter
                    meteor_files.append(file)
                    camera_names.append(file.split(".")[-2])
                    meteor_file_counter += 1
                    # A meteor has been found so set the variable to true
                    found_meteor = True
                    
                # if the files are all sorted (which hopefully
                # they are due to the sort method used above the loop
                # header),
                # it can be assumed that the desired .txts are all in their
                # own "block"
                elif meteor not in file and found_meteor == True:
                    break
                else:
                    continue

            # Gets the f-value/parameter of the meteor's light curve
            f_value = float(f_vals[meteor_index])
            
            # Gets the mass distribution from the f-parameter using
            # Murray & Beech's parameterization values (2000 & 2003)
            s = -2.2735 * f_value + 2.907
            
            # NOTE TO SELF: 1.5 - 3.5, make guess +- 0.2, step size 0.05

            # --------------------------------------------------------

            # This small section checks to see if there is actually
            # camera data for the meteor. If not, go straight to
            # writing the limited information from SummaryMeteorLog.txt
            # and continue to the next meteor.

            if(header_tracker and
               not any(meteor in file for file in image_txts)):

                # Calls empty_writer to do the writing to file
                empty_writer(folder_dir, meteor_index, meteor, dates,
                             time_strings, f_value, s)

                meteor_index += 1
                continue
            
            elif(not any(meteor in file for file in image_txts)):
                meteor_index += 1
                continue

            # --------------------------------------------------------

            # Saves s to the constants_change dictionary
            constants_change["s"] = s

            # Initializes lists to hold the meteor's times of detection, heights,
            # velocities, and fluxes as detected from a specific camera.
            # The latitudes and longitudes are also included for velocity
            # calculations in case of the velocity model used to create
            # the column in the camera data failing.
            # Each nested list represents each camera's set of data
            times = []
            heights = []
            velocities = []
            fluxes = []
            lats = []
            longs = []

            # For each meteor file, read it and return the relevant data
            # as lists
            for file in meteor_files:
                file_dir = os.path.join(folder_dir, file)
                (met_times, met_heights, met_vels, met_fluxes,
                 met_lats, met_longs) = meteor_reader(file_dir)
                
                # Append the relevant data (in the form of lists) into the
                # appropriate lists initialized above. Converts them to
                # floats and converts the UT time to seconds (also floats)
                # Notice the conversion from km/s to m/s in velocities
                times.append(time_converter(met_times))
                heights.append(list(np.float_(met_heights)))
                velocities.append(list(np.float_(met_vels) * 1000))
                fluxes.append(list(np.float_(met_fluxes)))
                lats.append(list(np.float_(met_lats)))
                longs.append(list(np.float_(met_longs)))

            # Creates an outer list to store the information in a different
            # format. In this case, the initialized list holds nested
            # numpy arrays. Each numpy array represents a camera
            # and holds the times, heights, velocities, and fluxes as
            # well
            nec_meteor_info = []
            for i in range(len(times)):
                # This groups each camera's information into their own
                # numpy arrays and puts each of those arrays into the
                # initialized list
                nec_meteor_info.append(np.array([times[i], heights[i],
                                                 fluxes[i]]))

            # Calls flux_time_driver to get a smoothed out version of
            # the time values and fluxes
            (time_smooth, dump, flux_smooth, avged_fluxes,
             avged_heights) = \
             flux_time_driver(nec_meteor_info,
                              meteor_file_counter)

            # ----------------------------------------------------------

            # Calls param_modifier to use averager() and make_smooth()
            # in a similar manner to flux_time_driver(). This is to match
            # the lengths of the returned lists and the ones returned
            # by flux_time_driver(). mod is short for modified
            (mod_times, mod_vel, mod_heights) = \
                        param_modifier(times, velocities, heights,
                                       meteor_file_counter)

            # ----------------------------------------------------------

            # Gets the scalar to multiply the initial mass by if needed
            # to account for the surviving mass
            mass_scalar = K["mass_scalar"]

            # Holds the total photometric mass of the meteor
            init_photomass = 1e-12

            # Saves the photomass of the meteor to plot the model
            # and observed data below
            photomass = 1e-12

            # For each flux value...
            for i in range(time_smooth.shape[0] - 1):
                # Calculate the observed intensities
                I_obs = I_0 * 10 ** (flux_smooth[i] / (-2.5))
                # Calculate a time step dt
                dt = time_smooth[i+1] - time_smooth[i]
                # Calculate the energy emitted over that time step
                E_emit = I_obs * dt * 4 * pi * r ** (2) * wavelength_range
                # Calculate the mass lost due to the energy emission
                dm = E_emit * 2 / (lum_eff * mod_vel[i] ** 2)
                # Add that lost mass to the photometric mass count
                init_photomass += dm

            # If there are less than 4 measured velocities,
            # set a variable from constants.py to True (I have no
            # idea what this is supposed to do)
            if(len(mod_vel) < 4):
                b_all_zeros = True

            # If there are more than 4 velocities...
            else:
                # Gets the average value of the first four velocities
                # and sets it as the meteor's initial velocity
                constants_change["vel_inf"] = np.mean(mod_vel[0:4])

                # Sets the photometric mass as the initial photometric
                # mass
                photomass = init_photomass

                # If the surviving mass is required (or not I guess),
                # multiply the photometric mass by a scalar between 0 and 1.0.
                if(~b_iter_mass_scalar):
                    photomass *= mass_scalar

                # Counts the number of iterations in the nested
                # for loops below
                iteration = 0

                # Initializes a list with 10 nested lists to hold the
                # parameters slowly closing in on the best values
                vals = [[],[],[],[],[],[],[],[],[],[]]

                # These variables hold the smallest flux and velocity
                # standard deviations between the measured and the
                # model's fluxes and velocities
                min_flux_sd = inf
                min_vel_sd = inf

                #print(photomass) 3.4e4
                new_const = constants_change

                #"""
                new_const["bulk_dens"] = 1200
                new_const["grain_dens"] = 3000
                new_const["heat_abl"] = 2.0e6 # 2 - 9e6
                new_const["glue_temp"] = 980 # 900 - 1600
                new_const["boil_temp"] = 1750 # 1400 - 2300
                new_const["fusion_temp"] = new_const["boil_temp"] - 100
                new_const["thermal_cond"] = 0.5 # 0.1 - 1.0
                new_const["atomic_m"] = 23 * 1e-3 / (6.02e23) # 23, 36, 56
                new_const["spec_heat"] = 1200 # 600 - 1400
                new_const["erosion_coeff"] = 0.15e-6# (0.2 - 0.8)e-6. Upper limit from paper
                new_const["s"] = 1.6 # 1.50 - 3.5
                photomass = 3.4e-4
                """

                """
                (flux_sd, vel_sd) = get_mass(new_const,
                                             True, photomass,
                                             time_smooth,
                                             flux_smooth, mod_times,
                                             mod_vel,
                                             mod_heights,
                                             zeniths[meteor_index],
                                             step_size,
                                             meteor, working_dir,
                                             camera_names,
                                             avged_fluxes,
                                             avged_heights)
                input()
                #"""
                
                
# ============================================================================

                # This section finds the best parameters that fit the
                # meteor's light curve. There's almost certainly a
                # way to use SciPy to optimize the parameters used
                # instead of seven layers of loops
                
                # Start of the nested loops in the search of the
                # best parameters...
                # In these loops, the copied dictionary has values added to
                # it. These values are used in creating different models
                # in get_mass().
                constants_change["erosion_coeff"] = 0.2e-6
                for bulk_dens in density_array:
                    # Meteor density
                    constants_change["bulk_dens"] = bulk_dens
                    # Particle density
                    constants_change["grain_dens"] = bulk_dens

                    for heat_abl in heat_abl_array:
                        # Heat of ablation
                        constants_change["heat_abl"] = heat_abl

                        for glue_temp in glue_array:
                            # Temperature at which the glue fails
                            constants_change["glue_temp"] = glue_temp

                            for boil_temp in boil_array:
                                # Meteor boiling temperature
                                constants_change["boil_temp"] = boil_temp
                                # Temperature at which meteor material fuses
                                # together
                                constants_change["fusion_temp"] = boil_temp - 100

                                for thermal_cond in thermal_array:
                                    # Meteor thermal conductivity
                                    constants_change["thermal_cond"] = thermal_cond

                                    for atomic_m in atomic_m_array:
                                        # Meteor vapor average atomic mass
                                        constants_change["atomic_m"] = atomic_m

                                        for spec_heat in spec_heat_array:
                                            # Meteor specific heat
                                            constants_change["spec_heat"] = spec_heat

                                            # Counts the iteration early
                                            iteration += 1

                                            # Gets the flux and velocity standard
                                            # deviation fo the model vs the
                                            # measured values
                                            (flux_sd, vel_sd) = get_mass(constants_change,
                                                                         False, photomass,
                                                                         time_smooth,
                                                                         flux_smooth, mod_times,
                                                                         mod_vel,
                                                                         mod_heights,
                                                                         zeniths[meteor_index],
                                                                         step_size,
                                                                         meteor, working_dir,
                                                                         camera_names,
                                                                         avged_fluxes,
                                                                         avged_heights)

                                            # Saves the used values in the
                                            # initialized list with 10 nested
                                            # lists above
                                            vals[0].append(bulk_dens)
                                            vals[1].append(heat_abl)
                                            vals[2].append(glue_temp)
                                            vals[3].append(boil_temp)
                                            vals[4].append(thermal_cond)
                                            vals[5].append(atomic_m)
                                            vals[6].append(spec_heat)
                                            vals[7].append(flux_sd)
                                            vals[8].append(vel_sd)
                                            vals[9].append(iteration)

                                            # If the standard deviations
                                            # are the new lowest, save them
                                            if(flux_sd < min_flux_sd):
                                                min_flux_sd = flux_sd

                                            if(vel_sd < min_vel_sd):
                                                min_vel_sd = vel_sd

# ============================================================================
            
            # This section is outside of the nested loops. It does the
            # final calculations with the data from the section above
            # in order to get final values that will be written to
            # the output file

            # Creates lists to store the averages and standard deviations
            # of each iterated parameter (density, boiling temperature,
            # etc...)
            averages = []
            sds = []

            # If not b_all_zeros (still no idea what that does)...
            if(not b_all_zeros):

                # Sorts the saved values from the giant nested
                # loops by least to greatest flux standard deviations
                sorted_indices = np.argsort(vals[7])
                for i in range(len(vals)):
                    vals[i] = np.array(vals[i])[sorted_indices]

                # Used to see how tolerant the values of each parameter
                # to be published should be. Essentially, if top_n is the
                # entire set of used parameters, the value published will be
                # the average of the values produced over the entire
                # parameter space.
                top_n = ceil(vals[0].shape[0])

                # If top_n > the specified max_n in constants.py...
                if(top_n > max_n):
                    # Set top_n to max_n
                    top_n = max_n

                # For each sorted flux standard deviation, compare
                # it to the maximum allowed standard deviation in
                # constants.py. If it is larger than the allowed value,
                # lower top_n by 1 and stop looking through the
                # flux standard deviations
                for i in range(top_n):
                    if(i >= 3 and vals[7][i] > max_sd):
                        top_n -= 1
                        break

                # Use top_n to get the nth values with the lowest
                # standard deviations and average them together.
                # The standard deviations themselves are subject
                # to getting averaged as well
                for i in range(9):
                    averages.append(mean(vals[i][0:(top_n)]))
                    sds.append(np.std(vals[i][0:(top_n)]) * \
                               top_n / (top_n - 1))

# ============================================================================

            constants_change["bulk_dens"] = averages[0]
            constants_change["grain_dens"] = averages[0]
            constants_change["heat_abl"] = averages[1]
            constants_change["glue_temp"] = averages[2]
            constants_change["boil_temp"] = averages[3]
            constants_change["fusion_temp"] = averages[3] - 100
            constants_change["thermal_cond"] = averages[4]
            constants_change["atomic_m"] = averages[5]
            constants_change["spec_heat"] = averages[6]
            
            (_, _) = get_mass(constants_change, True, photomass,
                              time_smooth, flux_smooth, mod_times,
                              mod_vel, mod_heights,
                              zeniths[meteor_index], step_size,
                              meteor, working_dir, camera_names,
                              avged_fluxes, avged_heights)

# ============================================================================


            # This section is responsible for the file writing

            if(header_tracker):
                # Opens the file to append into
                fit_summary = open(os.path.join(folder_dir, \
                                                    "FitSummary.txt"), "a+")

                # Writes the specific meteor's best calculated paramters
                # into the file with similar spacings as the header. The
                # order is in the same order as the headers.
                # The following is the order of the information written
                # into the file, if needed again:
                # ---------------------------------------------------
                # meteor number, date of detection, first time of detection,
                # initial velocity, photometric mass, light curve f-parameter,
                # linear regression related to the f-parameter,
                # meteor's density and its standard deviation[0],
                # heat of ablation and its standard deviation[1],
                # glue failure temperature and its standard deviation[2],
                # meteor's boiling temperature and its standard deviation[3],
                # meteor's thermal conductivity and its standard deviation[4],
                # meteor vapor's average atomic mass and its standard
                # deviation[5],
                # meteor's specific heat and its standard deviation[6],
                # flux and velocity standard deviations, and the number
                # of points to get the averages above
                # ---------------------------------------------------

                fit_summary.write("{: <7} {: <11} {: <12} {: <7.3} ".format(
                    meteor, dates[meteor_index],
                    time_strings[meteor_index],
                    (constants_change["vel_inf"] / 1000)))
                fit_summary.write("{: <12.4E} {: <7} {: <7.4} {: <9.2f} ".format(
                    init_photomass, f_value, constants_change[""],
                    averages[0]))
                fit_summary.write("{: <9.3f} {: <9.2E} {: <9.2E} {: <9} ".format(
                    sds[0], averages[1], sds[1], averages[2]))
                fit_summary.write("{: <7.3f} {: <9} {: <9.3f} {: <9.3} ".format(
                    sds[2], averages[3], sds[3], averages[4]))
                fit_summary.write("{: <7.3} {: <9} {: <9.3} {: <9} ".format(
                    sds[4], (averages[5] * avagadro_num * 10 ** (3)),
                    (sds[5] * avagadro_num * 10 ** (3)),
                    averages[6]))
                fit_summary.write("{: <9.3f} {: <7.4} {: <7.4} {: <8}\n".format(
                    sds[6], averages[7], averages[8], top_n))
                
                fit_summary.close()
            
            # Adds 1 to the meteor index.
            meteor_index += 1

# ============================================================================

# This section just contains the writer subfunctions

def header_writer(folder_dir):
    """
    Writes the header of the output file. The header is just the column
    titles with a row of seperating dashes.

    Parameter(s):
    folder_dir - The full path to the folder being worked in

    Return:
    None
    """
    
    # Creates the file to be appended into
    fit_summary = open(os.path.join(folder_dir, \
                                                "FitSummary.txt"), "w")

    # A list of lists containing the header of the file
    table_header = [["Meteor", "Observed", "Ref Time",
                     "Vinf", "Photometric", "F-skew", "Grain",
                     "Density", "+/-", "Heat Abl", "+/-",
                     "Temp", "+/-", "Boiling", "+/-",
                     "Th Cond", "+/-", "AtomMass", "+/-",
                     "SpecHeat", "+/-", "Flux",
                     "Veloc", "Number"],
                    ["Number", "Date", "(UT)", "(km/s)",
                     "Mass (kg)", "Ratio", "Distri", "(kg/m^3)",
                     "sigma", "(J/kg)", "sigma", "Glue (K)",
                     "sigma", "Temp (K)",  "sigma", "(W/m/K)",
                     "sigma", "(g/mol)", "sigma", "(J/kg/K)",
                     "sigma", "sigma", "sigma", "of Fits", ]]

    # For each row in table_header...
    for row in table_header:
        # Write to file the file with specific spacings
                    # so it looks organized when viewed by a person
        fit_summary.write(("{: <7} {: <11} {: <12} {: <7} " +
                           "{: <12} {: <7} {: <7} {: <9} " +
                           "{: <9} {: <9} {: <9} {: <9} " +
                           "{: <7} {: <9} {: <9} {: <9} " +
                           "{: <7} {: <9} {: <9} {: <9} " +
                           "{: <9} {: <7} {: <7} {: <8} \n").format(*row))

    # Puts in some dividers below the header and above the
    # information to be printed. 
    fit_summary.write("-" * 6 + "  " + "-" * 10 + "  " + \
                      "-" * 11 + "  " + "-" * 6 + "  " + \
                      "-" * 11 + "  " + "-" * 6 + "  " + \
                      "-" * 6 + "  " + "-" * 8 + "  " + \
                      "-" * 8 + "  " + "-" * 8 + "  " + \
                      "-" * 8 + "  " + "-" * 8 + "  " + \
                      "-" * 6 + "  " + "-" * 8 + "  " + \
                      "-" * 8 + "  " + "-" * 8 + "  " + \
                      "-" * 6 + "  " + "-" * 8 + "  " + \
                      "-" * 8 + "  " + "-" * 8 + "  " + \
                      "-" * 8 + "  " + "-" * 6 + "  " + \
                      "-" * 6 + "  " + "-" * 7 + " \n")

    # Closes the file
    fit_summary.close()

def empty_writer(folder_dir, meteor_index, meteor, dates, time_strings, f_value,
                 s):
    """
    If for some reason there isn't any camera data for a meteor, write
    out the information given from SummaryMeteorLog.txt

    Parameter(s):
    folder_dir   - The full path to the folder being worked in
    meteor_index - An integer to access the elements in each read column
                   of SummaryMeteorLog.txt
    meteor       - A list of meteor numbers
    dates        - A list of dates corresponding to said meteors
    time_strings - A list of time values corresponding to said meteors
    f_value      - A list of f-values/parameters corresponding ...
    s   - A value comes from an equation involving the
                   meteor's f-value

    Return:
    None
    """

    # Opens the file to write into
    fit_summary = open(os.path.join(folder_dir,"FitSummary.txt"), "a+")

    # Writes the given data. A dash (-) stands in for the information
    # that cannot be obtained without the camera data
    fit_summary.write("{: <7} {: <11} {: <12} {: <7} ".format(
        meteor, dates[meteor_index], time_strings[meteor_index],
        "-"))
    fit_summary.write("{: <12} {: <7} {: <7.4} {: <9} ".format(
        "-", f_value, s, "-"))
    fit_summary.write("{: <9} {: <9} {: <9} {: <9} ".format(
        "-", "-", "-", "-"))
    fit_summary.write("{: <7} {: <9} {: <9} {: <9.3} ".format(
        "-", "-", "-", "-"))
    fit_summary.write("{: <7} {: <9} {: <9} {: <9} ".format(
        "-", "-", "-", "-"))
    fit_summary.write("{: <9} {: <7} {: <7} {: <8}\n".format(
        "-", "-", "-", "-"))

    # Closes the file
    fit_summary.close()

main()                   
                            

            
