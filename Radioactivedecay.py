''' Radioactive decay functionality '''

import numpy
import Boosts
import Fourvector
import Particle
import PDG_ID
import random
import MRS
import sys
import Constants
import scipy.optimize
import Availabledecays
import copy
import Stringdecay
import Decaydictionary
#Decaytables = Decaydictionary.Decays()


def radioactivedecay4D(tau):
    ''' Function that finds a time for a decay of a particle to take place using monte carlo methods '''
    # need tau, rnd number, thts it? - maybe take the minimum distance - write two functions.
    rnd = random.random()
    Time = -1.*tau*numpy.log(1-rnd)
    
    return Time
    
    
    
''' Test code '''



