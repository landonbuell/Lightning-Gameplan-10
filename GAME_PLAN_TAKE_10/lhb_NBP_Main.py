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
import lhb_NBP_Phases as Phases
import time

            
if __name__ == '__main__':

    paths,bound,filelist = \
        Phases.PHASE_I()                    # Execute Phase I

    t_0 = time.process_time()               # current process time
    
    Phases.PHASE_II(filelist,bound,paths)   # Execute Phase II

    t_f = time.process_time ()              # current process time

    Phases.Base.Process_Times(\
        t_0,t_f,len(filelist))              # time data

    print("\nTotal Program Time:",time.process_time())