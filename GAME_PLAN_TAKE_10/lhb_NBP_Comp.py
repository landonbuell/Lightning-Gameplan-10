"""
Landon Buell
Lightning Research
Game Plan Take 10
30 Oct 2018
"""
   
            #################
            #### IMPORTS ####

import lhb_NBP_Base as Base
import numpy as np
import os

            ###################################
            #### COMPOSITE LEVEL FUNCTIONS ####

def Initialize ():
    """
    Initializes full script
    --------------------------------
    Return dictionary of directories
    """
    int_dir = os.getcwd()                           # establish initial directory
    
            #### Enter & Test Reading Directory Path ####
    while True:                                     # Setup to establish a reading directory
        #read = Base.Input_Directory('READ FROM')    # accept user input
        read = 'C:/Users/Landon/Documents/Lightning Research/data_rawINTF_201-400'
        path = Base.Change_Directory(read)          # attempt to change path       
        if path == True:                            # if successful,
            break                                   # break the loop
   
            #### Enter & Test Reading Directory Path ####
    while True:                                     # Setup to establish a reading directory
        #write = Base.Input_Directory('WRITE TO')    # accept user input
        write = 'C:/Users/Landon/Documents/Lightning Research/INTF_201-400'
        path = Base.Change_Directory(write)         # attempt to change path
        if path == True:                            # if successful,
            break                                   # break the loop
        else:                                       # if failure (path DNE)
            path = Base.Input_Create_Directory()    # prompt user to create path
            if path == True:                        # if yes,
                os.mkdir(write)                     # make the directory
                break

            #### Sub-Directory Paths ####
    paths = Base.Make_Sub_Dirs(write)               # create sub paths and dictionary to store them
    paths['readdir'] = read                         # add read dir to paths dictionary
    paths['intdir'] = int_dir                       # add int dir to paths dictionary
    paths['writedir'] = write                       # add write dir to paths dictionary
    return paths                                    # return dictionary of needed directories

def Analyze_I (file,bnd,dir_dict):
    """
    First layer of analysis for INTF data files. Examies nature of files
    --------------------------------
    file (str) : name of file to 'decode' into array
    bnd (float): upper bound criteria to serve limit low amplitudes
    dir_dict (dict) : dictionary of important directory paths
    --------------------------------
    returns  list of valid INTF files and invalid INTF files
    """
    os.chdir(dir_dict['readdir'])       # change to reading directory path
    if os.path.getsize(file) > 1e8:     # size file larger than 100,000 KB
        return [],False                 # return the empty array and False boolean
    data = Base.ATVT(file)              # Convert file in np array
    data = Base.Modify_Filedata(data)   # adjust data to set parameters
    os.chdir(dir_dict['intdir'])        # change back to initial directory
    event,low,high = \
        Base.Extract_Event(data)        # create an event based on max abs val
    if max(abs(event)) < bnd:           # If max of abs value less than bound
        return [] , False               # return the empty array and False boolean
    event = np.append(event,(low,high)) # add bounds to event array
    return event , True                 # return the evnt otherwise w/ True Booleain

def Analyze_II (data):
    """
    Second Layer of analysis for INTF data. Examine for NBP qualities
    --------------------------------
    data (array) : 1D array of event data. Last 2 idx are bounds of the array
    --------------------------------
    returns a matrix of [xdata,ydata] and True/False if pass/fail respectivly
    """
    x_data = np.arange(data[-2],data[-1])       # make x axis data
    y_data = data[:-2]                          # INTF data is all but last 2 idx

    y_smooth = Base.Curve_Smoother(y_data)      # smooth out the curve  
    markers = Base.Place_Markers(y_data)        # create a list of amplitude marker
    matrix = np.array([x_data,y_data])          # make data into single matrix

    """ Test for Quiet time - 10 us to 20 us"""
    quiet = Base.Test_Quiet_Time(y_data,2000)   # Test quiet time of the y-data set
    print("\tQuiet Time test:",quiet)           # pass/fail quiet time test
    if quiet == False:                          # of condition fails,
        print("\tRise Time test: N/A")          # message
        return matrix,False,markers             # return the matrix,False exit function

    """Test for Rise time - 1.8 us to 3.8 us"""
    rise = Base.Test_Rise_Time(markers,20)         # Test for valid rise time
    print("\tRise Time test:",rise)             # pass/fail quiet time test
    if rise == False:                           # if rise time condition fails:
        return matrix,False,markers             # return the matrix & False exit function
    """If all tests passed"""
    return matrix,True,markers                  # return the 2xN matrix if all tests pass