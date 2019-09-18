# Imports the os module for path manipulation
import os

def file_finder(folder_dir):
    """
    Skims a specified FOLDER for image MEAS .txts and a SummaryMeteorLog.txt.
    It then passes back a list containing the list of image .txts and the
    SummaryMeteorLog.txt.

    This function should be called within a loop that goes to all
    available folders in a directory.
    """

    # Initializes an empty string to hold the summary file
    summary_txt = ""

    # Initializes an empty list to hold all of the image .txts
    image_txts = []

    # For each file in the folder...
    for file in os.listdir(folder_dir):
        # If the file is an image .txt,
        if "MEAS" in file and file.endswith(".txt"):
            # Put it in the list to be returned
            image_txts.append(file)
        # If the file is the summary file,
        elif file == "SummaryMeteorLog.txt":
            # Save its name in the variable to be returned
            summary_txt = file
        # Continue to the next file if the current one isn't
        # a needed file
        else:
            continue

    # Return the names of the image .txts and the summary .txt
    return (image_txts, summary_txt)
