"""
Landon Buell
Lightning Research
Game Plan Take 10
30 Oct 2018
"""

            #################
            #### IMPORTS ####

import lhb_NBP_Base as Base
import lhb_NBP_Comp as Comp
import numpy as np
import os

            ######################
            #### USER PAHSE I ####

def PHASE_I ():
    """
    Phase I serves to initalize the program. All starting paramters are
    inputted or defaulted to here. 
    """
    paths = Comp.Initialize()                   # Dictionary of important paths
    #bound = Base.Input_Threshold()              # User input amplitude threshold
    bound = 4.0
    filelist = Base.Extract_Files(paths['readdir'],'.chD')
    print("\nThere are",len(filelist),"files in the requested directory path.")
    print("Expected Computation time:")
    print("Between",int((len(filelist)*0.5)),\
        "and",int((len(filelist)*0.8)),"seconds")

    return paths,bound,filelist

def PHASE_II (filelist,bnd,dir_dict):
    """
    Phase II serves as the filter and testing functions.
    """
    for I in range (len(filelist)):         # for all of the files

        os.chdir(dir_dict['readdir'])       # always change back to the read dir
        name = 'File_'+str(I+1)+'_'+str(filelist[I])
        print('-'*33)                       # spacer
        print(name)                         # title
        
            #### Initial Testings ####
        event , tests = Comp.Analyze_I(filelist[I],bnd,dir_dict)  
                                                # extract an event from the datafile
        if tests == False:                      # if failed either test (last idx is False)
            print("\tInitial tests: Failed")    # indicated failure of int. tests
            print("\tQuiet time tests: N/A")     
            print("\tRise time tests: N/A")     
            continue                            # next itteration of loop
        else:
            print("\tInitial tests: Passed")# indicated success of int. tests

            #### Primary Testings ####    
        matrix,tests,markers = Comp.Analyze_II(event)# test the event for NBP qualities
        xdata = matrix[0]                           # extract x-axis data
        ydata = matrix[1]                           # extract y-axis data
        """
        print("Marked Idxs:",markers)
        print("\tDifference",int(markers[-2]-markers[1]))
        Base.Markers_Graph(xdata,ydata,markers,name,show=True)
        """
        if tests == True:                           # if all tests passes,
            os.chdir(dir_dict['passplots'])         # change to pass plots path
            Base.Single_Graph(xdata,ydata,name,save=True)
            os.chdir(dir_dict['passarrays'])        # change to pass arrays
            Base.Write_Array(xdata,ydata,name)      # write array to .txtfile

        if tests == False:                          # if a test failed
            os.chdir(dir_dict['failplots'])         # change to pass plots path
            Base.Single_Graph(xdata,ydata,name,save=True)
            os.chdir(dir_dict['failarrays'])        # change to pass arrays
            Base.Write_Array(xdata,ydata,name)      # write array to .txtfile
