''' Baryon content module '''

import numpy
#import Boosts
import Fourvector
import Particle
import PDG_ID
import random
import MRS
import sys
import Constants
#import Stringdecay
import copy
import Twobodydecay
import Decaydictionary


class Baryoncontent(object):
    def __init__(self,quark1,quark2,quark3):
        self.__quark1 = quark1
        self.__quark2 = quark2
        self.__quark3 = quark3
        
    def Getquark1(self):
        return self.__quark1
        
    def Getquark2(self):
        return self.__quark2
        
    def Getquark3(self):
        return self.__quark3
        
    def Getquarks(self):
        return [self.__quark1,self.__quark2,self.__quark3]
        
        
class Baryontables(object):
    
    def __init__(self,quark1,quark2,quark3):
        
        self.__Baryoncontent = {
        
        2212: Baryoncontent(2,2,1),
        2112: Baryoncontent(2,1,1),
        
        2224: Baryoncontent(2,2,2),
        2214: Baryoncontent(2,2,1),
        2114: Baryoncontent(2,1,1),
        1114: Baryoncontent(1,1,1),
        
        3122: Baryoncontent(3,2,1),
        
        3222: Baryoncontent(3,2,2),
        3212: Baryoncontent(3,2,1),
        3112: Baryoncontent(3,1,1),
        3224: Baryoncontent(3,2,2),
        3214: Baryoncontent(3,2,1),
        3114: Baryoncontent(3,1,1),
        
        3322: Baryoncontent(3,3,2),
        3312: Baryoncontent(3,3,1),
        3324: Baryoncontent(3,3,2),
        3314: Baryoncontent(3,3,1),
        
        3334: Baryoncontent(3,3,3),
        
        }
        
        
    def Getcontent(self,code):
        code = abs(code)
        if self.__Baryoncontent.has_key(code):
            content = self.__Baryoncontent[code].Getquarks()
            
            return content
            
        print "Baryon content does not contain that code"
        
        
    
        
        