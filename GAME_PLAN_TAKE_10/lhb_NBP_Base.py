"""
Landon Buell
Lightning Research
Game Plan Take 10
30 Oct 2018
"""

            #################
            #### IMPORTS ####

import mas_raw
import numpy as np
import matplotlib.pyplot as plt
import os 

            ##############################
            #### USER INPUT FUNCTIONS ####

def Input_Directory (text):
    """
    Accept user input to identify a directory path
    -------------------------------
    text (str) : text to complete a input question sentence
    -------------------------------
    returns a user inputted directory path
    """
    while True:
        try:
            path = str(input("\nPlease enter the directory you wish to "+text+" : "))
            if not path:
                print("\n\tERROR! - Please input a directory")
            if path == ' ':
                print("\n\tERROR! - Please input a value")
            else:
                return path
                break
        except:
            print("\n\tERROR! - Inavlid value type")

def Input_Create_Directory ():
    """
    Accept user confirmation to create a firectory path
    -------------------------------
    path (str): directory path to create if desired
    -------------------------------
    returns True/False if yes/no respectivly
    """
    yes = ['y','yes','Y','YES']
    no = ['n','no','N','NO']
    while True:
        try:
            path = str(input("\nCreate this directory path?: "))
            if not path:
                print("\n\tERROR! - Please input a value!")
            if path in yes:
                return True
                break
            if path in no:
                return False
                break
            else:
                print("\n\tERROR! - That input is not valid")
        except:
            print("\n\tERROR! - Inavlid value type")

def Input_Threshold ():
    """
    Accept user input to idenify a minimum threshold for E-Feild
    -------------------------------
    Raises error if input is less than or equal to zero or greater than 160
    Returns an E-Feild boundary otherwise
    """
    print("\nThe Electric Feild Bound should be between 0 & 160")
    print("\t\tValues between 4-10 will produce the best results")
    while True:
        try:
            value = float(input("Input the E-Feild Threshold Limitation: "))
            if not value:
                print("\n\tERROR! - Please input a value")
            if value >= (160):
                print("\n\tERRROR! - That value is too large")
            if value <= 0:
                print("\n\tERRROR! - That value is too small")
            else:
                break
        except:
            print("\n\tERROR! - Inavlid value type")
    return value

            ##############################
            #### PRINT , PLOT & GRAPH ####

def Process_Times(int,fin,num):
    """
    Print Out timing data for program
    -------------------------------
    int (float) : inital script time reading
    fin (float) : final script time reading
    num (int) : number of files processed
    """
    print("\n",'-'*64)
    print("Process time for",num,"files:")
    print("\t",(fin-int),"seconds")
    print("Averaged",((fin-int)/num),"seconds per file")

def Single_Graph(x_data,y_data,title,save=False,show=False):             
    """
    Save a .png file of a plot  of 1 list for external use
    -------------------------------
    x_data (array): list to serve as the x-values for the plot
    y_data (array): list to serve as the y-values for the plot
    title (str): title of the plot an name to save as
    save (bool): Save plot to current directory as .png if True. False by default
    show (bool): Show plot of the current figure if True. False by default
    -------------------------------
    """
            ### Initialize ###
    plt.figure(figsize=(16,12))                     # set figure & size
    plt.ticklabel_format(style='sci',\
                axis='x', scilimits=(0,0))          # Sci notation on x-axis
    plt.title(title,fontsize=36,fontweight='bold')  # add title
    plt.xlabel("Data Index",\
        fontsize=20,fontweight='bold')              # label x-axis
    plt.ylabel("E-Feild Amp [V/m]",\
        fontsize=20,fontweight='bold')              # label y-axis
            ### Recreate Indices for X-axis ###
    plt.plot(x_data,y_data,color='purple')          # plot the dataset
    plt.grid(True)                                  # Create gridlines
    if save == True:                                # if the plot wants to be saved
        plt.savefig(title+'.png')                   # save the event figure
    if show == True:                                # if the plot wants to be shown
        plt.show()                                  # display the plot
    plt.close()                                     # close the figure

def Markers_Graph(x_data,y_data,xmarks,title,save=False,show=False):             
    """
    Save a .png file of a plot w/ markers of 1 list for external use. 
    -------------------------------
    x_data (array): list to serve as the x-values for the plot
    y_data (array): list to serve as the y-values for the plot
    title (str): title of the plot an name to save as
    save (bool): Save plot to current directory as .png if True. False by default
    show (bool): Show plot of the current figure if True. False by default
    -------------------------------
    """
            ### Initialize ###
    plt.figure(figsize=(16,12))                     # set figure & size
    plt.ticklabel_format(style='sci',\
                axis='x', scilimits=(0,0))          # Sci notation on x-axis
    plt.title(title,fontsize=36,fontweight='bold')  # add title
    plt.xlabel("Data Index",\
        fontsize=20,fontweight='bold')              # label x-axis
    plt.ylabel("E-Feild Amp [V/m]",\
        fontsize=20,fontweight='bold')              # label y-axis
    ymarks = np.array([y_data[I] for I in xmarks])  # create y axis markers       
    xmarks = np.array([x_data[I] for I in xmarks])  # create x axis markers
    plt.plot(x_data,y_data,color='purple')          # plot the dataset
    plt.plot(xmarks,ymarks,'ro')          # plot rise marks
    plt.grid(True)                                  # Create gridlines
    if save == True:                                # if the plot wants to be saved
        plt.savefig(title+'.png')                   # save the event figure
    if show == True:                                # if the plot wants to be shown
        plt.show()                                  # display the plot
    plt.close()                                     # close the figure

def Write_Array(x_data,y_data,name):
    """
    Save a .txt file of xdata,ydata
    --------------------------------
    xdata (array) : 1D array of x-axis data, to be written to txt file
    ydata (array) : 1D array of y-axis data, to be written to txt file
    name (str) : name to indentify .txt file by
    --------------------------------
    """
    matrix = np.array([x_data,y_data],dtype=str)# put into one matrix, make strings
    matrix = np.transpose(matrix)               # tranpose matrix into cols
    outfile = open(name+'.txt',mode='w')        # create the text file in write mode
    for I in range (len(matrix)):               # in the matrix 
        outfile.write(matrix[I][0])             # write x value
        outfile.write(',')                      # delimiter
        outfile.write(matrix[I][1])             # write y-value
        outfile.write('\n')                     # newline
    outfile.close()                             # close the file
            ###############################
            #### OPERATIONAL FUNCTIONS ####

def ATVT (file):               
    """
    Converts a Channel A,B,C or D file into a numpy array for further use
    -------------------------------
    file (str) : name of a raw file to be 'decoded'
    -------------------------------
    returns A 1D numpy array of datapoints
    """
    atvt = mas_raw.atvt()               # initiate the mas_raw Library
    atvt.load_file(file)                # decode the file
    data = atvt.raw_i                   # Create a datalist with the information
    data = np.array(data,dtype=float)   # Converts data into an array of floats
    return data                         # Return the list

def Change_Directory (directory):
    """
    Attempt to change working directory path
    -------------------------------
    directory (str) : directory to attempt a change to
    -------------------------------
    returns True/False for success/failure respectivly
    """
    try:                                        # Attempt
        os.chdir(directory)                     # change the path
        print("\n\tSuccessfully changed Directory paths")
        return True                             # return True
    except:                                     # if fails
        print("\n\tERROR! - Could not change directory path!")
        return False                            # return False

def Curve_Smoother (datalist,box_pts=10):
    """
    Smooths an Array through convolution
    -------------------------------
    datalist (array) : original dataset to be smoothed
    box_pts (int) : size of box of convolution
    -------------------------------
    returns a convolved, 'smoothed' array
    """
    box = np.ones(box_pts) / box_pts
    smoothed = np.convolve(datalist,box,mode='same')
    return smoothed

def Extract_Event (data,low=4000,up=6000):
    """
    Isolate maximum amplitude and create smaller dataset to work with
    -------------------------------
    data (array) : list or array of INTF E-feild data
    low (int) : Number of indicies to store before max amp (rec.4000-6000)
    up (int) : Number of indicies to store after max amplitude (rec.6000-14000)
    -------------------------------
    returns an np array of condensed set of data
    """
    max_idx = np.argmax(np.abs(data))   # index of maximum value
    try:                                # try to make low bound
        low_bnd = max_idx - low         # back low idx
    except:                             # if out of range
        low_bnd = data[0]               # low bound  is start of data
    try:                                # try to make upper bound
        up_bnd = max_idx + up           # forward up idx
    except:                             # if out of range
        up_bnd = data[-1]               # upper bound is end of data      
    """Bounds have been established, make an event from the dataset"""
    event = data[low_bnd:up_bnd]        # event from data chunck   
    return event,low_bnd,up_bnd         # return the array & bounds

def Extract_Files (directory,extension):
    """
    Find all files of a certain extension in a given directory path
    -------------------------------
    directory (str) : full directory path to seach given file types for
    extension (str) : indicated the desired file types to search for
    -------------------------------
    Returns a list of files (str) that contains all the files the designated type 
    """
    filelist = []                                   # Empty list to store filenames
    for root, dirs, files in os.walk(directory):    # elements in directory
        for file in files:                          # in those elements           
            if file.endswith(extension):            # If desired file type
                filelist.append(file)               # add to list of files           
    return filelist                                 # Return the list of files

def Make_Sub_Dirs (path):
    """
    Creates a series of directories within a specified path
    -------------------------------
    path (str) : parent directory to create paths within
    -------------------------------
    returns a dictionary of of sub directory paths
    """
            #### Create all Sub dirs ####
    passplots = path+'/NBEs/Plots'
    passarrays = path+'/NBEs/Arrays'
    failplots = path+'/Non-NBEs/Plots'
    failarrays = path+'/Non-NBEs/Arrays'
    paths = [passplots,passarrays,failplots,failarrays]
    
    for path in paths:                  # for each directory path
        try:                            # try to create path:
            os.makedirs(path)           # make the path (And other dirs)
        except:                         # if failure,
            pass                        # do nothing

    dict = {'passplots':passplots,
            'passarrays':passarrays,    
            'failplots':failplots,
            'failarrays':failarrays}    # create a dictionary of directory paths
    return dict                         # return the dictionary

def Place_Markers (data):      
    """
    Isolates markers in a given dataset that meet amplitude criteria
    -------------------------------
    data (array) : 1D array to test rise time of
    -------------------------------
    Returns a list of indexs of 0% to 90% risetimes
    """
    data = np.absolute(data)
    max_amp = max(data)                 # max of datalist
    max_idx = np.argmax(data)           # start at maximum amp index
    marks = []                          # list to store markers
    N = 0.0                             # set int N (0%)
    idx = 0
    while N <= 1:                       # while less than 100%
        amp = N*max_amp                 # desired amplitude is N% of marks
        for I in range (max_idx):       # search all idx up to max amp
            if data[I] <= amp:          # if less than or equal to desired amp
                idx = I                 # that is the desired index
            else:                       # otherwise
                continue                # skip the next itterations
        marks.append(idx)               # add the index to the list of idxs
        N += 0.1                        # incriment by 10%
    marks.sort()                        # sort min to max values   
    return marks                        # return the markers list

def Modify_Filedata(data):
    """
    Format the E-Feild np array for algorithm
    -------------------------------
    data (array) : original dataset , decoded
    -------------------------------
    Returns the modifed dataset for computations
    """
    data = data * (-1)              # flip over x-axis
    average = np.average(data)      # find average value of the dataset
    data -= average                 # Zero the dataset
    """160 V/m = 2^15 Digital units
    1 DU = 0.0048828125 V/m"""
    data *= (0.0048828125)          # convert all data to V/m
    return data                     # Returnt the dataset

            ###########################
            #### TESTING FUNCTIONS ####

def Test_Quiet_Time (data,index=1000,max_amp=0.1):
    """
    Tests first chosen number of idx of an array for quiet time
    -------------------------------
    datalist (array) : 1D array to test quiet time of
    index (int) : idx number to test noise up to (1500 by default)
    max_amp (float) : Amp percentage threshold (rec 0.1 - 0.4) - (set to 0.1 by default)
    -------------------------------
    Returns True / False if pass/fail respectivly
    """
    if max_amp <= 0.0:                  # if max <= 0% (unreasonable)
        max_amp = 0.01                  # set to 1%
    if max_amp >= 1.0:                  # if max >= 100% (unreasonable)
        max_amp = 0.99                  # set to 99%
    data = np.absolute(data)            # absolute value of data
    total_max = max(data)               # max of abs. of data
    sample = data[:index]               # extract set of 'N' idxs
    if max(sample) > max_amp*total_max: # if max of sample is > N% total max
        return False                    # fails quiet time conditionmax
    else:                               # otherwise
        return True                     # passes quiet time condition

def Test_Rise_Time(marks,min=100,max=1000):
    """
    Test a dataset fot NBP - like rise times and pre-positioned markers
    -------------------------------
    marks (array) : list of amplitude value markers
    min (int) : minimum idx difference for NBP rise time (100 by default)
    min (int) : maximum idx difference for NBP rise time (1500 by default)
    -------------------------------
    Returns True / False if pass/fail respectivly
    """
    srt = marks[1]              # index of 10% value
    end = marks[-2]             # index of 90% value
    """For an NBP, the rise time should be 180 idx to 380 us"""
    diff = abs(end - srt)       # find the index difference
    if diff < min:              # if less than minimum idx diff
        return False            # the rise is too fast to be NBP
    if diff > max:              # if more than maximum idx diff
        return False            # the rise is too slow for NBP
    else:                       # otherwise
        return True             # passed test!