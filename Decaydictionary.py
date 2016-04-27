''' Available Decays v2. Making the system more general '''

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
import Baryoncontent


class Decayinput(object):
    def __init__(self, *args):
        self.__arguments = args
        
    def Getdecays(self):
        return self.__arguments
        
    def Printdecays(self):
        print self.__arguments
        
        
class Baryonwfinput(object):
    def __init__(self, *args):
        self.__arguments = args
        
    def Getwfs(self):
        return self.__arguments
        
    def Getdiquarks(self):
        diquarks = []
        for i in range(len(self.__arguments)):
            diquarks.append(self.__arguments[i][0])
        return diquarks
            
        
        
class Mesoncontent(object):
    def __init__(self,quark_no,antiquark_no):
        self.__qID = quark_no
        self.__qbarID = antiquark_no
        
    def GetqID(self):
        return self.__qID
        
    def GetqbarID(self):
        return self.__qbarID
        
class Failsafe2probs(object): # redundant (semi) as has the same format as decayinput. Merely a different name so easi8er to look back on.
    def __init__(self, *args):
        self.__arguments = args
        
    def Getprobs(self):
        return self.__arguments
        
    def Printprobs(self):
        print self.__arguments


class Probabilities(object):
    def __init__(self,mass,weight,wf,scaling): # wavefunction here is the squared value - is in all the code everywhere
        self.__mass      = mass
        self.__weight    = weight
        self.__wavefunc  = wf
        self.__scaling   = scaling
        
    def Printall(self):
        print "mass =!", self.__mass, "weight =", self.__weight, "wf =", self.__wavefunc, "scaling =", self.__scaling
        
    def Getmass(self):
        return self.__mass
        
    def Getweight(self):
        return self.__weight
    
    def Getwf(self):
        return self.__wavefunc
        
    def Getscaling(self):
        return self.__scaling
        

class Decays(object):
    ''' Class containing dictionaries? of the available products producible by quarks at the ends of a string '''
    ''' Negative numbers are antiquarks '''
    ''' DECAYS ARE MASS ORDERED - HIGHEST MASS PRODUCT -> LOWEST MASS PRODUCT '''
    def Getmass(self,code):
        code = abs(code)
        code = int(code)
        #code = float(code)
        ID = PDG_ID.PDG_Type(code)
        ##print ID, "ID in getmass"
        ##print abs(ID), "abs of the ID in getmass"
        return PDG_ID.PDG_Type.mass(ID) # code or ID?
    
    def Getweight(self,code):
        mass = self.Getmass(code)
        weight = numpy.exp(-mass*mass/(Constants.kperpsigma()*Constants.kperpsigma())) # 4 factor 
        return weight
        
    def GetweightFS2(self,code):
        mass = self.Getmass(code)
        weight = numpy.exp(-mass*mass/(Constants.kperpsigma()*Constants.kperpsigma()))
        return weight
        
    def Getweightbaryon(self,code):
        mass = self.Getmass(code)
        weight = numpy.exp(-1* Constants.baryon_exp_factor() * mass*mass/(Constants.kperpsigma()*Constants.kperpsigma()))
        return weight
        
    def Getweightdiquark(self,code):
        mass = self.Getmass(code)
        weight = numpy.exp(-1* Constants.diquark_exp_factor() * mass*mass/(Constants.kperpsigma()*Constants.kperpsigma()))
        return weight
        
    def Finalprob(self,code): # non-mixed states (not eta / eta prime)
        code = abs(code)
        if self.__probinfo.has_key(code):
            weight  = self.__probinfo[code].Getweight()
            wf      = self.__probinfo[code].Getwf()
            scaling = self.__probinfo[code].Getscaling()
            prob = weight*wf*scaling
            return prob
        return "incorrect meson number input"
       
    def Finalprob_mixed(self,code,quark_no): # quark no == 1,2,3 for d,u,s respectively
        code = abs(code)
        if self.__probinfo.has_key(code):
            weight  = self.__probinfo[code].Getweight()
            wf      = self.__probinfo[code].Getwf()[abs(quark_no) - 1]
            scaling = self.__probinfo[code].Getscaling()
            prob = weight*wf*scaling
            return prob
        return "incorrect meson input, or meson is not a mixed state"
        
    def Getbaryonwf(self,baryoncode,diquarkcode):
        baryoncode = abs(baryoncode)
        diquarkcode = abs(diquarkcode)
        if self.__baryonwfs.has_key(baryoncode):
            wfinfo = self.__baryonwfs[baryoncode].Getwfs()
            for i in range(len(wfinfo)):
                if diquarkcode == wfinfo[i][0]:
                    return wfinfo[i][1]
                    
            #print "not found in that baryon wf set, check for errors"
            
    def Getbaryonleftover(self,baryoncode,diquarkcode):
        antitest = baryoncode
        baryoncode = abs(baryoncode)
        diquarkcode = abs(diquarkcode)
        if self.__baryonwfs.has_key(baryoncode):
            wfinfo = self.__baryonwfs[baryoncode].Getwfs()
            for i in range(len(wfinfo)):
                if diquarkcode == wfinfo[i][0]:
                    quarkflavour = wfinfo[i][2]
                    #if antitest > 0:
                        #quarkflavour = quarkflavour *-1
                    return quarkflavour
                        
                    # return the leftover quark flavour
                    
            #print "not found in that baryon wf set, check for errors"
                
            
        
    def Finalprob_baryon(self,baryoncode,diquarkcode):
        baryoncode = abs(baryoncode)
        diquarkcode = abs(diquarkcode)
        if self.__probinfo.has_key(baryoncode): # write that in as a check.
            weight = self.Getweightbaryon(baryoncode)
            wf = self.Getbaryonwf(baryoncode,diquarkcode)
            scaling = self.__probinfo[baryoncode].Getscaling()            
            prob = weight * wf * scaling
            
            return prob
        return "incorrect baryon input, or baryon wf function working incorrectly"
        
    def Finalprob_diquark(self,diquark):
        diquark = abs(diquark)
        if self.__probinfo.has_key(diquark):
            weight  = self.__probinfo[diquark].Getweight()
            wf      = self.__probinfo[diquark].Getwf()
            scaling = self.__probinfo[diquark].Getscaling()
            prob = weight*wf*scaling
            return prob
        return "incorrect meson input, or meson is not a mixed state"
        
        
    def __init__(self):
        self.__decays = { # mass ordered from heaviest on left to lightest on right
            1:   Decayinput(331,313,223,113,-213,3101,3103,1103,2101,2103,221,311,-211,111),
            2:   Decayinput(331,323,223,113,213,3201,3203,2203,2101,2103,221,321,211,111),
            3:   Decayinput(333,331,3303,-313,-323,3201,3203,3101,3103,221,-311,-321),
           -1:   Decayinput(331,-313,223,113,213,-3101,-3103,-1103,-2101,-2103,221,-311,211,111),
           -2:   Decayinput(331,-323,223,113,-213,-3201,-3203,-2203,-2101,-2103,221,-321,-211,111),
           -3:   Decayinput(333,331,-3303,313,323,-3201,-3203,-3101,-3103,221,311,321),
           
         1103:   Decayinput(3114,1114,2114,3112,2112), # dd(1)
         2101:   Decayinput(3122,2112,2212), # ud(0)
         2103:   Decayinput(3214,2214,2114,3212,2112,2212), # ud(1)
         2203:   Decayinput(3224,2224,2214,3222,2212), # uu(1)
         3101:   Decayinput(3312,3112,3212,3122), # sd(0)
         3103:   Decayinput(3314,3114,3214,3312,3112,3212,3122), # sd(1)
         3201:   Decayinput(3322,3212,3222,3122), # su(0)
         3203:   Decayinput(3324,3214,3224,3322,3212,3222,3122), # su(1)
         3303:   Decayinput(3334,3314,3324,3312,3322), # ss(1)
         
         #Anti diquarks. Everything is antimattered as no "self inverse" baryons exist!
        -1103:   Decayinput(-3114,-1114,-2114,-3112,-2112), # dd(1)
        -2101:   Decayinput(-3122,-2112,-2212), # ud(0)
        -2103:   Decayinput(-3214,-2214,-2114,-3212,-2112,-2212), # ud(1)
        -2203:   Decayinput(-3224,-2224,-2214,-3222,-2212), # uu(1)
        -3101:   Decayinput(-3312,-3112,-3212,-3122), # sd(0)
        -3103:   Decayinput(-3314,-3114,-3214,-3312,-3112,-3212,-3122), # sd(1)
        -3201:   Decayinput(-3322,-3212,-3222,-3122), # su(0)
        -3203:   Decayinput(-3324,-3214,-3224,-3322,-3212,-3222,-3122), # su(1)
        -3303:   Decayinput(-3334,-3314,-3324,-3312,-3322) # ss(1)
           
           }
                  
        self.__masses = {
            1:   Decayinput(self.Getmass(331),self.Getmass(313),self.Getmass(223),self.Getmass(113),self.Getmass(-213),self.Getmass(221),self.Getmass(311),self.Getmass(-211),self.Getmass(111)),
            2:   Decayinput(self.Getmass(331),self.Getmass(323),self.Getmass(223),self.Getmass(113),self.Getmass(213),self.Getmass(221),self.Getmass(321),self.Getmass(211),self.Getmass(111)),
            3:   Decayinput(self.Getmass(333),self.Getmass(331),self.Getmass(-313),self.Getmass(-323),self.Getmass(221),self.Getmass(-311),self.Getmass(-321)),
           -1:   Decayinput(self.Getmass(331),self.Getmass(-313),self.Getmass(223),self.Getmass(113),self.Getmass(213),self.Getmass(221),self.Getmass(-311),self.Getmass(211),self.Getmass(111)),
           -2:   Decayinput(self.Getmass(331),self.Getmass(-323),self.Getmass(223),self.Getmass(113),self.Getmass(-213),self.Getmass(221),self.Getmass(-321),self.Getmass(-211),self.Getmass(111)),
           -3:   Decayinput(self.Getmass(333),self.Getmass(331),self.Getmass(313),self.Getmass(323),self.Getmass(221),self.Getmass(311),self.Getmass(321)),
           
           
         1103:   Decayinput(self.Getmass(3114),self.Getmass(1114),self.Getmass(2114),self.Getmass(3112),self.Getmass(2112)), # dd(1)
         2101:   Decayinput(self.Getmass(3122),self.Getmass(2112),self.Getmass(2212)), # ud(0)
         2103:   Decayinput(self.Getmass(3214),self.Getmass(2214),self.Getmass(2114),self.Getmass(3212),self.Getmass(2112),self.Getmass(2212)), # ud(1)
         2203:   Decayinput(self.Getmass(3224),self.Getmass(2224),self.Getmass(2214),self.Getmass(3222),self.Getmass(2212)), # uu(1)
         3101:   Decayinput(self.Getmass(3312),self.Getmass(3112),self.Getmass(3212),self.Getmass(3122)), # sd(0)
         3103:   Decayinput(self.Getmass(3314),self.Getmass(3114),self.Getmass(3214),self.Getmass(3312),self.Getmass(3112),self.Getmass(3212),self.Getmass(3122)), # sd(1)
         3201:   Decayinput(self.Getmass(3322),self.Getmass(3212),self.Getmass(3222),self.Getmass(3122)), # su(0)
         3203:   Decayinput(self.Getmass(3324),self.Getmass(3214),self.Getmass(3224),self.Getmass(3322),self.Getmass(3212),self.Getmass(3222),self.Getmass(3122)), # su(1)
         3303:   Decayinput(self.Getmass(3334),self.Getmass(3314),self.Getmass(3324),self.Getmass(3312),self.Getmass(3322)), # ss(1)
         
         #Anti diquarks. Everything is antimattered as no "self inverse" baryons exist!
        -1103:   Decayinput(self.Getmass(-3114),self.Getmass(-1114),self.Getmass(-2114),self.Getmass(-3112),self.Getmass(-2112)), # dd(1)
        -2101:   Decayinput(self.Getmass(-3122),self.Getmass(-2112),self.Getmass(-2212)), # ud(0)
        -2103:   Decayinput(self.Getmass(-3214),self.Getmass(-2214),self.Getmass(-2114),self.Getmass(-3212),self.Getmass(-2112),self.Getmass(-2212)), # ud(1)
        -2203:   Decayinput(self.Getmass(-3224),self.Getmass(-2224),self.Getmass(-2214),self.Getmass(-3222),self.Getmass(-2212)), # uu(1)
        -3101:   Decayinput(self.Getmass(-3312),self.Getmass(-3112),self.Getmass(-3212),self.Getmass(-3122)), # sd(0)
        -3103:   Decayinput(self.Getmass(-3314),self.Getmass(-3114),self.Getmass(-3214),self.Getmass(-3312),self.Getmass(-3112),self.Getmass(-3212),self.Getmass(-3122)), # sd(1)
        -3201:   Decayinput(self.Getmass(-3322),self.Getmass(-3212),self.Getmass(-3222),self.Getmass(-3122)), # su(0)
        -3203:   Decayinput(self.Getmass(-3324),self.Getmass(-3214),self.Getmass(-3224),self.Getmass(-3322),self.Getmass(-3212),self.Getmass(-3222),self.Getmass(-3122)), # su(1)
        -3303:   Decayinput(self.Getmass(-3334),self.Getmass(-3314),self.Getmass(-3324),self.Getmass(-3312),self.Getmass(-3322)) # ss(1)
           
           }
           
        ''' probinfo: mass, weight, wf^2, scaling '''   
        # Getweight used as. Could use GetweightFS2
        # These are objects producible at the ends of strings. mesons, diquarks and baryons(at diquark ends).
        self.__probinfo = {
            111:   Probabilities(self.Getmass(111),self.Getweight(111),Constants.pi0_wf2(),1),
            113:   Probabilities(self.Getmass(113),self.Getweight(113),Constants.pi0_wf2(),Constants.vector_scaling_factor()),
            211:   Probabilities(self.Getmass(211),self.Getweight(211),1,1),
            213:   Probabilities(self.Getmass(213),self.Getweight(213),1,Constants.vector_scaling_factor()),
            221:   Probabilities(self.Getmass(221),self.Getweight(221),Constants.eta_wf2(),Constants.eta_scaling_factor()),     #eta /eta' wf are tuples of (d,u,s) wfs respectively (in a list)
            223:   Probabilities(self.Getmass(223),self.Getweight(223),Constants.pi0_wf2(),Constants.vector_scaling_factor()), # omega has effectively the same wf as the pi/rho mesons.
            311:   Probabilities(self.Getmass(311),self.Getweight(311),1,Constants.strange_suppression_factor()),
            313:   Probabilities(self.Getmass(313),self.Getweight(313),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            321:   Probabilities(self.Getmass(321),self.Getweight(321),1,Constants.strange_suppression_factor()),
            323:   Probabilities(self.Getmass(323),self.Getweight(323),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            331:   Probabilities(self.Getmass(331),self.Getweight(331),Constants.eta_prime_wf2(),Constants.eta_prime_scaling_factor()),
            333:   Probabilities(self.Getmass(333),self.Getweight(333),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            
            #Diquarks are producible at the ends of the string. These are set such that abs(code taken due to invarince of sign)
            
           1103:   Probabilities(self.Getmass(1103),self.Getweightdiquark(1103),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor()), # dd(1)
           2101:   Probabilities(self.Getmass(2101),self.Getweightdiquark(2101),1,Constants.diquark_scaling_factor()), # ud(0)
           2103:   Probabilities(self.Getmass(2103),self.Getweightdiquark(2103),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor()), # ud(1)
           2203:   Probabilities(self.Getmass(2203),self.Getweightdiquark(2203),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor()), # uu(1)
           3101:   Probabilities(self.Getmass(3101),self.Getweightdiquark(3101),1,Constants.diquark_scaling_factor() * Constants.strange_diquark_scaling_factor()), # sd(0)
           3103:   Probabilities(self.Getmass(3103),self.Getweightdiquark(3103),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor()), # sd(1)
           3201:   Probabilities(self.Getmass(3201),self.Getweightdiquark(3201),1,Constants.diquark_scaling_factor() * Constants.strange_diquark_scaling_factor()), # su(0)
           3203:   Probabilities(self.Getmass(3203),self.Getweightdiquark(3203),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor()), # su(1)
           3303:   Probabilities(self.Getmass(3303),self.Getweightdiquark(3303),1,Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor()), # ss(1)              
            
            # Baryons
            
           2212:   Probabilities(0,0,0,Constants.baryon_scaling_factor()), # dud values put in. These are calculated elsewhere on a percase basis.
           2112:   Probabilities(0,0,0,Constants.baryon_scaling_factor()),
           2224:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           2214:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           2114:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           1114:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           
           3122:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           
           3222:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           3212:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           3112:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           3224:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           3214:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           3114:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           
           3322:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           3312:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor()),
           3324:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           3314:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           
           3334:   Probabilities(0,0,0,Constants.baryon_scaling_factor() * Constants.strange_baryon_scaling_factor() * Constants.vector_baryon_scaling_factor()),
           
           
            
            }
            
        self.__baryonwfs = { # input is a set of lists, i.e. wavefunction overlap with that baryon.
        # Format or arguments are [diquark, wf overlap, leftover quark flavour]
        
           2212:   Baryonwfinput([2203,1./3,1], [2103,1./6,2], [2101,0.5,2]),
           2112:   Baryonwfinput([1103,1./3,2], [2103,1./6,1], [2101,0.5,1]),
           3222:   Baryonwfinput([3203,1./6,2], [3201,0.5,2], [2203,1./3,3]),
           3112:   Baryonwfinput([3103,1./6,1], [3101,0.5,1], [1103,1./3,3]),
           3212:   Baryonwfinput([2103,1./3,3], [3203,1./12,1], [3201,0.25,1], [3103,1./12,2], [3101,0.25,2]),
           3122:   Baryonwfinput([2101,1./3,3], [3203,0.25,1], [3201,1./12,1], [3103,0.25,2], [3101,1./12,2]),
           3322:   Baryonwfinput([3203,1./6,3], [3201,0.5,3], [3303,1./3,2]),
           3312:   Baryonwfinput([3103,1./6,3], [3101,0.5,3], [3303,1./3,1]),
           
           2224:   Baryonwfinput([2203,1,2]),
           2214:   Baryonwfinput([2103,1./3,2], [2203,2./3,1]),
           2114:   Baryonwfinput([2103,1./3,1], [1103,2./3,2]),
           1114:   Baryonwfinput([1103,1,1]),
           3224:   Baryonwfinput([2203,2./3,3], [3203,1./3,2]),
           3214:   Baryonwfinput([2103,4./6,3], [3203,1./6,1], [3103,1./6,2]),
           3114:   Baryonwfinput([1103,2./3,3], [3103,1./3,1]),
           3324:   Baryonwfinput([3203,1./3,3], [3303,2./3,2]),
           3314:   Baryonwfinput([3103,1./3,3], [3303,2./3,1]),
           3334:   Baryonwfinput([3303,1,3]),
           
               
               }
        
            
        self.__finalprobs = {
            1:   Decayinput(self.Finalprob_mixed(331,1),self.Finalprob(313),self.Finalprob_mixed(223,1),self.Finalprob_mixed(113,1),self.Finalprob(-213),self.Finalprob_diquark(3101),self.Finalprob_diquark(3103),self.Finalprob_diquark(1103),self.Finalprob_diquark(2101),self.Finalprob_diquark(2103),self.Finalprob_mixed(221,1),self.Finalprob(311),self.Finalprob(211),self.Finalprob_mixed(111,1)),
            2:   Decayinput(self.Finalprob_mixed(331,2),self.Finalprob(323),self.Finalprob_mixed(223,2),self.Finalprob_mixed(113,2),self.Finalprob(213),self.Finalprob_diquark(3201),self.Finalprob_diquark(3203),self.Finalprob_diquark(2203),self.Finalprob_diquark(2101),self.Finalprob_diquark(2103),self.Finalprob_mixed(221,2),self.Finalprob(321),self.Finalprob(211),self.Finalprob_mixed(111,2)), # use abs for the decays - mass same.
            3:   Decayinput(self.Finalprob(333),self.Finalprob_mixed(331,3),self.Finalprob_diquark(3303),self.Finalprob(-313),self.Finalprob(-323),self.Finalprob_diquark(3201),self.Finalprob_diquark(3203),self.Finalprob_diquark(3101),self.Finalprob_diquark(3103),self.Finalprob_mixed(221,3),self.Finalprob(-311),self.Finalprob(-321)), # abs for mixed states qid, again symmetry argument to be made
           -1:   Decayinput(self.Finalprob_mixed(331,1),self.Finalprob(-313),self.Finalprob_mixed(223,1),self.Finalprob_mixed(113,1),self.Finalprob(213),self.Finalprob_diquark(-3101),self.Finalprob_diquark(-3103),self.Finalprob_diquark(-1103),self.Finalprob_diquark(-2101),self.Finalprob_diquark(-2103),self.Finalprob_mixed(221,1),self.Finalprob(-311),self.Finalprob(211),self.Finalprob_mixed(111,1)),
           -2:   Decayinput(self.Finalprob_mixed(331,2),self.Finalprob(-323),self.Finalprob_mixed(223,2),self.Finalprob_mixed(113,2),self.Finalprob(-213),self.Finalprob_diquark(-3201),self.Finalprob_diquark(-3203),self.Finalprob_diquark(-2203),self.Finalprob_diquark(-2101),self.Finalprob_diquark(-2103),self.Finalprob_mixed(221,2),self.Finalprob(-321),self.Finalprob(-211),self.Finalprob_mixed(111,2)), # use abs for the decays - mass same.
           -3:   Decayinput(self.Finalprob(333),self.Finalprob_mixed(331,3),self.Finalprob_diquark(-3303),self.Finalprob(313),self.Finalprob(323),self.Finalprob_diquark(-3201),self.Finalprob_diquark(-3203),self.Finalprob_diquark(-3101),self.Finalprob_diquark(-3103),self.Finalprob_mixed(221,3),self.Finalprob(311),self.Finalprob(321)), # abs for mixed states qid, again symmetry argument to be made
         
         1103:   Decayinput(self.Finalprob_baryon(3114,1103),self.Finalprob_baryon(1114,1103),self.Finalprob_baryon(2114,1103),self.Finalprob_baryon(3112,1103),self.Finalprob_baryon(2112,1103)),
         2101:   Decayinput(self.Finalprob_baryon(3122,2101),self.Finalprob_baryon(2112,2101),self.Finalprob_baryon(2212,2101)),
         2103:   Decayinput(self.Finalprob_baryon(3214,2103),self.Finalprob_baryon(2214,2103),self.Finalprob_baryon(2114,2103),self.Finalprob_baryon(3212,2103),self.Finalprob_baryon(2112,2103),self.Finalprob_baryon(2212,2103)),
         2203:   Decayinput(self.Finalprob_baryon(3224,2203),self.Finalprob_baryon(2224,2203),self.Finalprob_baryon(2214,2203),self.Finalprob_baryon(3222,2203),self.Finalprob_baryon(2212,2203)),
         3101:   Decayinput(self.Finalprob_baryon(3312,3101),self.Finalprob_baryon(3112,3101),self.Finalprob_baryon(3212,3101),self.Finalprob_baryon(3122,3101)),
         3103:   Decayinput(self.Finalprob_baryon(3314,3103),self.Finalprob_baryon(3114,3103),self.Finalprob_baryon(3214,3103),self.Finalprob_baryon(3312,3103),self.Finalprob_baryon(3112,3103),self.Finalprob_baryon(3212,3103),self.Finalprob_baryon(3122,3103)),
         3201:   Decayinput(self.Finalprob_baryon(3322,3201),self.Finalprob_baryon(3212,3201),self.Finalprob_baryon(3222,3201),self.Finalprob_baryon(3122,3201)),
         3203:   Decayinput(self.Finalprob_baryon(3324,3203),self.Finalprob_baryon(3214,3203),self.Finalprob_baryon(3224,3203),self.Finalprob_baryon(3322,3203),self.Finalprob_baryon(3212,3203),self.Finalprob_baryon(3222,3203),self.Finalprob_baryon(3122,3203)),
         3303:   Decayinput(self.Finalprob_baryon(3334,3303),self.Finalprob_baryon(3314,3303),self.Finalprob_baryon(3324,3303),self.Finalprob_baryon(3312,3303),self.Finalprob_baryon(3322,3303)),
         
        -1103:   Decayinput(self.Finalprob_baryon(-3114,-1103),self.Finalprob_baryon(-1114,-1103),self.Finalprob_baryon(-2114,-1103),self.Finalprob_baryon(-3112,-1103),self.Finalprob_baryon(-2112,-1103)),
        -2101:   Decayinput(self.Finalprob_baryon(-3122,-2101),self.Finalprob_baryon(-2112,-2101),self.Finalprob_baryon(-2212,-2101)),
        -2103:   Decayinput(self.Finalprob_baryon(-3214,-2103),self.Finalprob_baryon(-2214,-2103),self.Finalprob_baryon(-2114,-2103),self.Finalprob_baryon(-3212,-2103),self.Finalprob_baryon(-2112,-2103),self.Finalprob_baryon(-2212,-2103)),
        -2203:   Decayinput(self.Finalprob_baryon(-3224,-2203),self.Finalprob_baryon(-2224,-2203),self.Finalprob_baryon(-2214,-2203),self.Finalprob_baryon(-3222,-2203),self.Finalprob_baryon(-2212,-2203)),
        -3101:   Decayinput(self.Finalprob_baryon(-3312,-3101),self.Finalprob_baryon(-3112,-3101),self.Finalprob_baryon(-3212,-3101),self.Finalprob_baryon(-3122,-3101)),
        -3103:   Decayinput(self.Finalprob_baryon(-3314,-3103),self.Finalprob_baryon(-3114,-3103),self.Finalprob_baryon(-3214,-3103),self.Finalprob_baryon(-3312,-3103),self.Finalprob_baryon(-3112,-3103),self.Finalprob_baryon(-3212,-3103),self.Finalprob_baryon(-3122,-3103)),
        -3201:   Decayinput(self.Finalprob_baryon(-3322,-3201),self.Finalprob_baryon(-3212,-3201),self.Finalprob_baryon(-3222,-3201),self.Finalprob_baryon(-3122,-3201)),
        -3203:   Decayinput(self.Finalprob_baryon(-3324,-3203),self.Finalprob_baryon(-3214,-3203),self.Finalprob_baryon(-3224,-3203),self.Finalprob_baryon(-3322,-3203),self.Finalprob_baryon(-3212,-3203),self.Finalprob_baryon(-3222,-3203),self.Finalprob_baryon(-3122,-3203)),
        -3303:   Decayinput(self.Finalprob_baryon(-3334,-3303),self.Finalprob_baryon(-3314,-3303),self.Finalprob_baryon(-3324,-3303),self.Finalprob_baryon(-3312,-3303),self.Finalprob_baryon(-3322,-3303)),
         
         
         
            }
        
        
        self.__mesoncontent = {
            111:   Mesoncontent(0,0), # 0,0 used for mixed states - need a more sophisticated method - atm can only generate meson with same flavour as the decaying quark of teh string
            113:   Mesoncontent(0,0),
            #130:   Mesoncontent(1,-3),      # K0L supposed to be d sbar for now  MAY NEED TO CORRECT
            211:   Mesoncontent(2,-1),
            213:   Mesoncontent(2,-1),
            221:   Mesoncontent(0,0),
            223:   Mesoncontent(0,0),
            #310:   Mesoncontent(1,-3),     # CURRENTLY OMITTING THE K0L AND K0S FROM THE ARRAYS OF MESON CONTENT. NEUTRAL KAON STILL PRESENT WITH SAME EFFECTIVE CONTENT. ASSIGNMENT COMES LATER.
            311:   Mesoncontent(1,-3),
            313:   Mesoncontent(1,-3),
            321:   Mesoncontent(2,-3),
            323:   Mesoncontent(2,-3),
            331:   Mesoncontent(0,0),
            333:   Mesoncontent(0,0), # ACTUALLY 3,-3 BUT WELL SEE
            411:   Mesoncontent(4,-1),
            413:   Mesoncontent(4,-1),
            421:   Mesoncontent(4,-2),
            423:   Mesoncontent(4,-2),
            431:   Mesoncontent(4,-3),
            433:   Mesoncontent(4,-3),
            #441:   Mesoncontent(4,-4),
            #443:   Mesoncontent(4,-4),
            511:   Mesoncontent(1,-5),
            521:   Mesoncontent(2,-5),
            531:   Mesoncontent(3,-5),
            533:   Mesoncontent(3,-5),
            541:   Mesoncontent(4,-5),
            #551:   Mesoncontent(5,-5),
            #553:   Mesoncontent(5,-5),
            
           #-130:   Mesoncontent(3,-1), 
           -211:   Mesoncontent(1,-2), # ANTIPARTICLES / ANTIMESONS - EACH ANTISYMMETRIC MESON REQUIRES A NEGATIVE PARTNER.
           -213:   Mesoncontent(1,-2), # symmetric mesons have been ignored as to not double count mesons (111, 113, 221, 223, 331, 441, 443, 551, 553) do not have negative counterparts
           #-310:   Mesoncontent(3,-1),
           -311:   Mesoncontent(3,-1),
           -313:   Mesoncontent(3,-1),
           -321:   Mesoncontent(3,-2),
           -323:   Mesoncontent(3,-2),
           -411:   Mesoncontent(1,-4),
           -413:   Mesoncontent(1,-4),
           -421:   Mesoncontent(2,-4),
           -423:   Mesoncontent(2,-4),
           -431:   Mesoncontent(3,-4),
           -433:   Mesoncontent(3,-4),
           -511:   Mesoncontent(5,-1),
           -521:   Mesoncontent(5,-2),
           -531:   Mesoncontent(5,-3),
           -533:   Mesoncontent(5,-3),
           -541:   Mesoncontent(5,-4),
            
            }
            
        self.__Failsafe2probs = { # funtion for finding the w.f. overlap when calculating the probability in Failsafe 2 mechanism. - already squared as a value
        # All values are listed d u s in that order. when referring to them, think of quark no -1 
        # If value not present here, that prob is ==1 as it is a pure state
            
              1:   Failsafe2probs(111,113,221,223,331),
              2:   Failsafe2probs(111,113,221,223,331),
              3:   Failsafe2probs(221,331,333),
              4:   Failsafe2probs(441,443),
              5:   Failsafe2probs(551,553),
            
            # Probabilities (mass, weight, w.f., suppression factors.)
            
            111:   Probabilities(self.Getmass(111),self.GetweightFS2(111),Constants.pi0_wf2(),1), # should be a 1 at the end. trying something.
            113:   Probabilities(self.Getmass(113),self.GetweightFS2(113),Constants.pi0_wf2(),Constants.vector_scaling_factor()),
            211:   Probabilities(self.Getmass(211),self.GetweightFS2(211),1,1),
            213:   Probabilities(self.Getmass(213),self.GetweightFS2(213),1,Constants.vector_scaling_factor()),
            221:   Probabilities(self.Getmass(221),self.GetweightFS2(221),Constants.eta_wf2(),Constants.eta_scaling_factor()),     #eta /eta' wf are tuples of (d,u,s) wfs respectively (in a list)
            223:   Probabilities(self.Getmass(223),self.GetweightFS2(223),Constants.pi0_wf2(),Constants.vector_scaling_factor()),
            311:   Probabilities(self.Getmass(311),self.GetweightFS2(311),1,Constants.strange_suppression_factor()),
            313:   Probabilities(self.Getmass(313),self.GetweightFS2(313),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            321:   Probabilities(self.Getmass(321),self.GetweightFS2(321),1,Constants.strange_suppression_factor()),
            323:   Probabilities(self.Getmass(323),self.GetweightFS2(323),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            331:   Probabilities(self.Getmass(331),self.GetweightFS2(331),Constants.eta_prime_wf2(),Constants.eta_prime_scaling_factor()),
            333:   Probabilities(self.Getmass(333),self.GetweightFS2(333),1,Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            411:   Probabilities(self.Getmass(411),self.GetweightFS2(411),1,Constants.charm_scaling_factor()),
            413:   Probabilities(self.Getmass(413),self.GetweightFS2(413),1,Constants.charm_scaling_factor() * Constants.vector_scaling_factor()),
            421:   Probabilities(self.Getmass(421),self.GetweightFS2(421),1,Constants.charm_scaling_factor()),
            423:   Probabilities(self.Getmass(423),self.GetweightFS2(423),1,Constants.charm_scaling_factor() * Constants.vector_scaling_factor()),
            431:   Probabilities(self.Getmass(431),self.GetweightFS2(431),1,Constants.charm_scaling_factor() * Constants.strange_suppression_factor()),
            433:   Probabilities(self.Getmass(433),self.GetweightFS2(433),1,Constants.charm_scaling_factor() * Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            441:   Probabilities(self.Getmass(441),self.GetweightFS2(441),1,Constants.charm_scaling_factor()),
            443:   Probabilities(self.Getmass(443),self.GetweightFS2(443),1,Constants.charm_scaling_factor() * Constants.vector_scaling_factor()),
            511:   Probabilities(self.Getmass(511),self.GetweightFS2(511),1,Constants.bottom_scaling_factor()),
            521:   Probabilities(self.Getmass(521),self.GetweightFS2(521),1,Constants.bottom_scaling_factor()),
            531:   Probabilities(self.Getmass(531),self.GetweightFS2(531),1,Constants.bottom_scaling_factor() * Constants.strange_suppression_factor()),
            533:   Probabilities(self.Getmass(533),self.GetweightFS2(533),1,Constants.bottom_scaling_factor() * Constants.strange_suppression_factor() * Constants.vector_scaling_factor()),
            541:   Probabilities(self.Getmass(541),self.GetweightFS2(541),1,Constants.bottom_scaling_factor() * Constants.charm_scaling_factor()),
            551:   Probabilities(self.Getmass(551),self.GetweightFS2(551),1,Constants.bottom_scaling_factor()),
            553:   Probabilities(self.Getmass(553),self.GetweightFS2(553),1,Constants.bottom_scaling_factor() * Constants.vector_scaling_factor()),
            
            
            }
            
            
    def Getdecays(self,quark_no): # return the possible hadrons formed at a quark end
        if self.__decays.has_key(quark_no):
            return self.__decays[quark_no].Getdecays()
        return "quark number not dealt with"
        
    def Getmasses(self,quark_no): # returns the masses of the hadrons formed at a quark end
        ''' Returns dictionary entry for masses of the decay products for a quark end of a string '''
        if self.__masses.has_key(quark_no):
            return self.__masses[quark_no].Getdecays()
        return "quark number not dealt with"
           
           
    def Getfinalprobs(self,quark_no):
        ''' Returns dictionary entry for probabilities of the decay products for a quark end of a string '''
        if self.__finalprobs.has_key(quark_no):
            return self.__finalprobs[quark_no].Getdecays()
        return "quark number not dealth with"
        
    def Getmesoncontent(self,quark_no,meson_no):
        if self.__mesoncontent.has_key(meson_no):
            if self.__mesoncontent[meson_no].GetqID() == 0: # have eta/etaprime/pi0 - quark content is the qqbar of the quark flavour (by construction) - DESIGN CHOICE
                quark_number = abs(float(quark_no))
                return [quark_number, -1*quark_number] # the eta / etaprime / pi0 is quark/antiquark of the decaying quark
            quarkflavour = self.__mesoncontent[meson_no].GetqID()
            antiquarkflavour = self.__mesoncontent[meson_no].GetqbarID()
            return [float(quarkflavour), float(antiquarkflavour)]
        return "meson number not found in dictionary"
        
    def Getmindecay(self,quark_no):
        if self.__masses.has_key(quark_no):
            return self.__masses[quark_no].Getdecays()[-1]
        return "quark number not dealt with"
        
        
    def GetFailsafe2hadronlists(self,string):
        ''' Create lists, 2 entries. [hadron1 (quark end), hadron2 (antiquark end)] '''
        quarkIDs = string.Getstringcontent()
        qID = quarkIDs[0]
        qbarID = quarkIDs[1]
        possibledecaypairs = []
        
        xlist = numpy.arange(1,6)
        xbarlist = -1 * xlist # create the list of antiquarks to pair up with the string quark.
        
        ''' Create the hadron list '''
        IDs = PDG_ID.PDG_ID()                    
        hadronlist = IDs.Gethadronlist()
        hadronlist = copy.deepcopy(hadronlist)
        
        ''' create full hadron list - including negative mesons (antiparticles) '''
        neghadronlist = copy.deepcopy(hadronlist)
        for i in range(len(neghadronlist)):
            neghadronlist[i] = neghadronlist[i] *-1
        totalhadronlist = hadronlist + neghadronlist
        
        for i in range(len(xlist)):     # iterate over the number of quark pairings. this is done by looking at what is produced at the quark end and then the antiquark end of the string by constraint
            ##print xlist[i]
            # Find all possible qxbar pairs. append them into a list
            qxbarlist = []
            xqbarlist = [] # list of mesons to be filled
            # Consider the duplicate items, when q = xbar
            if qID == abs(xbarlist[i]): # not to be confused with xqbarlist :(
                qxbarlist = self.__Failsafe2probs[qID].Getprobs() # find the degen list
                
                # Then check for the duplicates on the other side of the string.
                if abs(qbarID) == xlist[i]: # if dup then find that list also
                    xqbarlist = self.__Failsafe2probs[xlist[i]].Getprobs()
                    
                    zippedlist = self.Failsafe2ziplists(qxbarlist,xqbarlist) # zip the lists
                    for k in range(len(zippedlist)):
                        possibledecaypairs.append(zippedlist[k]) # append the lists to the main list
                    
                    continue # skip the rest of this iteration - its done. go to next x value
                    
                # if xbar =/= q then do the normal method of searching over the meson list
                for q in range(len(totalhadronlist)): # go over the hadron list 
                    if self.__mesoncontent.has_key(totalhadronlist[q]):
                        Mesonq = self.__mesoncontent[totalhadronlist[q]].GetqID()
                        Mesonqbar = self.__mesoncontent[totalhadronlist[q]].GetqbarID() # get out the parts
                    
                    
                        if Mesonq == xlist[i] and Mesonqbar == qbarID:      # find the matches for string quark, and supposed antiquark
                            xqbarlist.append(totalhadronlist[q])
                
                
                zippedlist = self.Failsafe2ziplists(qxbarlist,xqbarlist)
                for k in range(len(zippedlist)):
                    possibledecaypairs.append(zippedlist[k])
                    
                continue # skip to the next x value
                
                
            # in the event that qID =/= abs(xbarlist[i])
            # Iterate over the list hadrons for the xvalue to build the qxbar mesons    
            
            for p in range(len(totalhadronlist)):
                
                if self.__mesoncontent.has_key(totalhadronlist[p]): # iterate over the list of hadrons / mesons in the dictionary
                    Mesonq = self.__mesoncontent[totalhadronlist[p]].GetqID()
                    Mesonqbar = self.__mesoncontent[totalhadronlist[p]].GetqbarID() # get out the parts
                    
                    # Find the qxbar and xqbar mesons
                    if Mesonq == qID and Mesonqbar == xbarlist[i]:      # find the matches for string quark, and supposed antiquark
                        qxbarlist.append(totalhadronlist[p])
                    if Mesonq == xlist[i] and Mesonqbar == qbarID:
                        xqbarlist.append(totalhadronlist[p])
            
            if abs(qbarID) == xlist[i]: # if dup then find that list also
                    xqbarlist = self.__Failsafe2probs[xlist[i]].Getprobs()
                    ##print xqbarlist, "xqbarlist"
                    
                    zippedlist = self.Failsafe2ziplists(qxbarlist,xqbarlist) # zip the lists
                    ##print zippedlist, "zippedlist"
                    for k in range(len(zippedlist)):
                        possibledecaypairs.append(zippedlist[k]) # append the lists to the main list
                    
                    continue # skip the rest of this iteration - its done. go to next x value            
                                                
            ##print qxbarlist, "qxbarlist"
            ##print xqbarlist, "xqbarlist"
            zippedlist = self.Failsafe2ziplists(qxbarlist,xqbarlist)
            ##print zippedlist, "zippedlist"
            for k in range(len(zippedlist)):
                possibledecaypairs.append(zippedlist[k])
        
        return possibledecaypairs
            
        
        
        
    def Failsafe2ziplists(self,list1,list2):
        listoflists = []
        for i in range(len(list1)):
            hadronpairlist = []
            for p in range(len(list2)):
                hadronpairlist.append(list1[i])
                hadronpairlist.append(list2[p])
                listoflists.append(hadronpairlist) # append the new pair
                hadronpairlist = [] # reintialise the lsit
        
        return listoflists
        
    def Removeheavymesons(self,string,listofdecaypairs):
        ''' Removes the charm and bottom mesons from the possible decay pairs if there is no charm / bottom quarks at the ends of the string '''
        
        quarkIDs = string.Getstringcontent()
        qID = quarkIDs[0]
        qbarID = quarkIDs[1]
        newlist = []
        
        if qID < 4 and abs(qbarID) < 4: # both string ends are light - no heavy mesons > 400 are allowed
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 400 and abs(listofdecaypairs[i][1]) < 400:
                    newlist.append(listofdecaypairs[i])
            return newlist
        
        if qID < 4 and abs(qbarID) == 4: # q is light, qbar is charm. qbar can form charmed mesons and below, q cannot
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 400 and abs(listofdecaypairs[i][1]) < 500: 
                    newlist.append(listofdecaypairs[i])
            return newlist
            
        if qID < 4 and abs(qbarID) == 5: # q is light, qbar is bottom. qbar can for bottomed mesons and below, q cannot
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 400 and abs(listofdecaypairs[i][1]) < 540:# qbarID can produce any of the allowed mesons - ccbar and bbar pairs cannot be produced in the vacumn.
                    newlist.append(listofdecaypairs[i])
            return newlist
            
        if qID == 4 and abs(qbarID) < 4: # q is charm, qbar is light
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 500 and abs(listofdecaypairs[i][1]) < 400:
                    newlist.append(listofdecaypairs[i])
            return newlist
            
        if qID == 4 and abs(qbarID) == 4: # both mesons charm, anything below 500 producible
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 500 and abs(listofdecaypairs[i][1]) < 500:
                    newlist.append(listofdecaypairs[i])
            return newlist
            
        if qID == 4 and abs(qbarID) == 5: # charm and bottom respectively qbar anyting , q below 500
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 500 and abs(listofdecaypairs[i][1]) < 540:
                    newlist.append(listofdecaypairs[i])
            return newlist
                    
        if qID == 5 and abs(qbarID) < 3: # bottom and light meson
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 540 and abs(listofdecaypairs[i][1]) < 400:
                    newlist.append(listofdecaypairs[i])
            return newlist
            
        if qID == 5 and abs(qbarID) == 4: # bottom and charm meson
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 540 and abs(listofdecaypairs[i][1]) < 500:
                    newlist.append(listofdecaypairs[i])
            return newlist
        
        if qID == 5 and abs(qbarID) == 5: # bottom and bottom meson
            for i in range(len(listofdecaypairs)):
                if abs(listofdecaypairs[i][0]) < 540 and abs(listofdecaypairs[i][1]) < 540:
                    newlist.append(listofdecaypairs[i])
            return newlist
        
        
        return newlist
                    
        
    def Removetoomassivemesons(self,string,flavourpossiblelist):
        ''' Deletes mesonic pairs that are not possible due to insufficient string mass '''
        stringmass = string.Getstringmass()
        finallist = []
        
        for i in range(len(flavourpossiblelist)):
            mass1 = self.Getmass(flavourpossiblelist[i][0])
            mass2 = self.Getmass(flavourpossiblelist[i][1])
            totalmass = mass1 + mass2
            if stringmass > totalmass:
                finallist.append(flavourpossiblelist[i])
                
        return finallist

    def CalculateprobsFS2(self,string,listofhadrons): # takes the list of pairs of MESOns that are possible, and returns a list of their probabilties of being found.
        
        quarkIDs = string.Getstringcontent() # needed for handling the etas etc - have to get the correct w.f. overlap.
        
            
        qID = quarkIDs[0]
        qbarID = quarkIDs[1]
        if qID < 0 and qbarID > 0:
            temp = qID
            qID == qbarID
            qbarID = temp
        probs = []
        mixedstatemesons = [111,113,221,223,331] # Hold the mixed states to iterate over and find the ones that I need to do something additional.
        
        for i in range(len(listofhadrons)): # list of pairs of hadrons
            hadron1code = abs(listofhadrons[i][0]) # Probabilities dont care about antihadrons - the masses and suppression factors / wavefunctions are the same as the real counterpart
            hadron2code = abs(listofhadrons[i][1]) # get out the hadron codes
            
            if self.__Failsafe2probs.has_key(hadron1code) == False:

                print "Hadron", hadron1code, "is not in the failsafe2probs table" # some checks for debugging and safety control.
            
            if self.__Failsafe2probs.has_key(hadron1code) == False:
                print "Hadron", hadron2code, "is not in FS2 table"
                
                       
            hadron1weight = self.__Failsafe2probs[hadron1code].Getweight()
            hadron2weight = self.__Failsafe2probs[hadron2code].Getweight()
            
            hadron1supps = self.__Failsafe2probs[hadron1code].Getscaling()
            hadron2supps = self.__Failsafe2probs[hadron2code].Getscaling()
            
            if hadron1code in mixedstatemesons:
                hadron1wf = self.__Failsafe2probs[hadron1code].Getwf()[int(abs(qID) - 1)]
                
            if hadron1code not in mixedstatemesons:
                hadron1wf = self.__Failsafe2probs[hadron1code].Getwf()
                
            if hadron2code in mixedstatemesons:
                hadron2wf = self.__Failsafe2probs[hadron2code].Getwf()[int(abs(qbarID) - 1)]
                
            if hadron2code not in mixedstatemesons:
                hadron2wf = self.__Failsafe2probs[hadron2code].Getwf()
            
            ##print hadron1supps, "supp 1"
            ##print hadron2supps, "supp 2"
            ##print hadron1weight, "weight 1"
            ##print hadron2weight, "weight 2"
            ##print hadron1wf, "wf 1 "
            ##print hadron2wf, "wf 2"

            finalprob = hadron1supps * hadron2supps * hadron1weight * hadron2weight * hadron1wf * hadron2wf
            
            probs.append(finalprob)
        
        return probs
        
    def ChooseFailsafe2hadron(self,listofpairs,probs): # takes the list of possible pairs and the list of probabilities as arguments.
        
        sumprobs = sum(probs)
        no_probs = len(probs) # how many probs there are
        normalised_probs = probs / sumprobs
        rng = numpy.random.choice(no_probs,1,p=normalised_probs) # chosen index in the array
        
        chosenpair = listofpairs[rng] # chosen pair corresponding to that index
        return chosenpair
        
        
        
    def Failsafe2complete(self,string): # Puts it through the ringer - returning a pair of hadrons. This is the chosen pair for string decay.
        
        if string.Getq().Checkquark() is True and string.Getqbar().Checkquark() is True:
            finalstate = self.DoublequarkFailsafe2(string)
            return finalstate
            
        if string.Getq().Checkdiquark() is True and string.Getqbar().Checkquark() is True:
            finalstate = self.DiquarkquarkFailsafe2(string)
            return finalstate
            
        if string.Getqbar().Checkdiquark() is True and string.Getq().Checkquark() is True:
            finalstate = self.DiquarkquarkFailsafe2(string)
            return finalstate
            
        #print "FAILSAFE2 COMPLETE HAS FAILED SLIGHTLY AND CAUGHT A CASE IT SHOULD NOT HAVE TO DEAL WITH"
        #print string.Getstringcontent(), "STRING CONTENT SEEN BY FAILSAFE2COMPLETE"
        
        ''' TESTING SOME TRUTHS, SEE WHAT HAPPENS '''
        
        #print string.Getq().Checkquark()
        #print string.Getq().Checkdiquark()
        #print string.Getqbar().Checkquark()
        #print string.Getqbar().Checkdiquark()
        
        
            
    def DiquarkquarkFailsafe2(self,string):
        
        ''' Find the diquark and the quark in the string '''
        if string.Getq().Checkdiquark() is True:
            if string.Getqbar().Checkquark() is False:
                #print "BEING FED AN INAPPROPRIATE STRING INTO DIQUARKQUARKFAILSAFE2"
                sys.exit()
            diquark = string.GetqID()
            quark = string.GetqbarID()
            
        if string.Getqbar().Checkdiquark() is True:
            if string.Getq().Checkquark() is False:
                #print "BEING FED AN INAPPROPRIATE STRING INTO DIQUARKQUARKFAILSAFE2"
                sys.exit()
            diquark = string.GetqbarID()
            quark = string.GetqID()
        
        diquarkneg = False # check for negative particles being produced.    
        if diquark < 0:
            diquarkneg = True
            
        quarkneg = False
        if quark < 0:
            quarkneg = True    
        
        baryons = self.Getdecays(diquark) # Get the diquark baryons. Produce each baryon first and then observe what is left.
        FS2pairs = []
        for i in range(len(baryons)):
            currentbaryon = baryons[i]
            remaining_quark = self.Getbaryonleftover(currentbaryon,diquark) # have the leftover. Just need to get the available mesons now.
            if diquarkneg == False: # positive diquark
                remaining_quark = remaining_quark * -1 # Have an antiquark produced. Otherwise unaffected
            if diquarkneg == True:
                currentbaryon = currentbaryon * -1 # Get an antibaryon
            
                

            temp_particle = Particle.particle4D(remaining_quark,string.Getqvector4D()) # any mom will do, just need a particle type.
            quark_particle = Particle.particle4D(quark,string.Getqbarvector4D())
            newstring = MRS.mrs4D(quark_particle,temp_particle)
            newstring = self.Reshuffleflavours(newstring) # Consolidate the flavours of the remaining quarks into a string to find the possible productions.
            mesons = self.GetFailsafe1mesonlist(newstring)
            for p in range(len(mesons)):
                FS2pairs.append([currentbaryon,mesons[p]]) # append each possible pairing of baryon and meson.
        
        #print FS2pairs, "FS2 PAIRS DQQFFS2"    
        newpairswithmass = self.Rdqmp(string,FS2pairs) # remove diquar quark massive pairs.

        pickapair = self.DiquarkquarkFailsafe2probs(diquark,quark,newpairswithmass) # used to be FS2pairs. trying to account for mass now
        if abs(pickapair[0]) == 311:
            pickapair[0] = self.Checkforsinglekaon(pickapair[0])

        if abs(pickapair[1]) == 311:
            pickapair[1] = self.Checkforsinglekaon(pickapair[1])
        #print pickapair, "PICKAPAIR DIQUARKQUARKFAILSAFE2"
        
        return pickapair
            
    def Rdqmp(self,string,pairs): # function for removing pairs from baryon meson pair athat are too massive.
        stringmass = string.Getstringmass()
        finallists = []

        for i in range(len(pairs)):
            mass1 = self.Getmass(pairs[i][0])
            mass2 = self.Getmass(pairs[i][1])
            totalmass = mass1 + mass2
            if stringmass > totalmass:
                finallists.append(pairs[i])

        return finallists
            
                
            
            
        
        
    def DiquarkquarkFailsafe2probs(self,diquark,quark,listoflists): # takes a list of [baryon,meson] lists
        
        finalproblist = []
        #print listoflists, "LISTOFLITS DIQUARKQUARKFAILSAFE2PROBS"
        for i in range(len(listoflists)):
            baryonprob = self.Finalprob_baryon(listoflists[i][0],diquark)
            mixedmesons = [111,113,221,223,331]
            if listoflists[i][1] in mixedmesons:
                mesonprob = self.Finalprob_mixed(listoflists[i][1],int(quark))
            if listoflists[i][1] not in mixedmesons:
                mesonprob = self.Finalprob(listoflists[i][1])
                
            finalprob = baryonprob * mesonprob
            finalproblist.append(finalprob)
        
        sumprobs = sum(finalproblist)
        no_probs = len(finalproblist) # how many probs there are
        normalised_probs = finalproblist / sumprobs
        rng = numpy.random.choice(no_probs,1,p=normalised_probs) # chosen index in the array
        
        chosenpair = listoflists[rng] # chosen pair corresponding to that index
        return chosenpair
             
        
        
    def DoublequarkFailsafe2(self,string):
        Decaytables = Decays()
        initial_lists = Decaytables.GetFailsafe2hadronlists(string) # Create all possible pairings
        
        remove_heavy = Decaytables.Removeheavymesons(string,initial_lists) # Remove the charmed and bottomed mesons that cannot be produced here
        remove_massive = Decaytables.Removetoomassivemesons(string,remove_heavy) # Remove the mesons pairs where the sum of the mass is greater than the string
        
        probabilities = Decaytables.CalculateprobsFS2(string,remove_massive)
        
        hadronpair = Decaytables.ChooseFailsafe2hadron(remove_massive,probabilities)
        #print hadronpair, "hadronpair as selected in soublequarkfailesafe2"
        
        if abs(hadronpair[0]) == 311 and abs(hadronpair[1]) == 311:
            #print "double hadronpair activated"
            
            hadronpair = self.Checkfordoublekaon(hadronpair)
            #print hadronpair, "hadron pair after the first check"
            if hadronpair[0] == 311 or hadronpair[1] == 311:
                print hadronpair, "hadronpair after the checkfordoublekaon stuff"
                #print hadronpair
                #print hadronpair[0]
                #print hadronpair[1]
                print "HERE IT IS"
                print "really is this happending?"
                sys.exit()
                
        if abs(hadronpair[0]) == 311 and abs(hadronpair[1]) != 311:
            #print "single kaon system engaged"
            hadronpair[0] = self.Checkforsinglekaon(hadronpair[0])
        
        if abs(hadronpair[1]) == 311 and abs(hadronpair[0]) != 311:
            #print "single KAON 2 SYSTEM REACHED"
            hadronpair[1] = self.Checkforsinglekaon(hadronpair[1])
            
            
        #print hadronpair
    
        return hadronpair
        
    def Findminandmaxpairmass(self,string): # takes a string, and finds all the possible pairs for a string.
        Decaytables = Decays()

        stringcontent = string.Getstringcontent()
        #print stringcontent, "STRINGCONTENT MINMAX MODULE"
        quark = stringcontent[0]
        quark = PDG_ID.PDG_Type(quark)
        antiquark = stringcontent[1]
        antiquark = PDG_ID.PDG_Type(antiquark)
        
        if quark.isQuark() is True and antiquark.isQuark() is True:
            return self.Doublequarkpairmass(string)
            
        if quark.isDiquark() is True and antiquark.isQuark() is True: # quark, diquark string
            return self.Quarkdiquarkpairmass(string)
            
        if quark.isQuark() is True and antiquark.isDiquark() is True: # diquark, quark string.
            return self.Quarkdiquarkpairmass(string)
            
        if quark.isDiquark() is True and antiquark.isDiquark() is True: # both ends are a diquark somehow
            decays = Decaytables.Getdecays(quark)
            heaviestbaryon = decays[0]
            heaviestbaryonmass = self.Getmass(heaviestbaryon)
            maxpairmass = 2*heaviestbaryonmass
            minpairmass = Constants.doublebaryonmassmodifier() * maxpairmass
            return maxpairmass, minpairmass
            
            
            

        
        
        
    def Doublequarkpairmass(self,string):
        Decaytables = Decays()
        string = self.Reshuffleflavours(string)
        #print string.Getstringcontent(), "STRING CONTENT IN DOUBLE QUARK MODULE"
        initial_lists = Decaytables.GetFailsafe2hadronlists(string) # Create all possible pairings
            
        remove_heavy = Decaytables.Removeheavymesons(string,initial_lists) # Remove the charmed and bottomed mesons that cannot be produced here
        #remove_massive = Decaytables.Removetoomassivemesons(string,remove_heavy) # Remove the mesons pairs where the sum of the mass is greater than the string    
        
        #dont remove the mesons that are too massive, as they are the limiting factor!
        totalmass = []
        for i in range(len(remove_heavy)):
            hadron1 = remove_heavy[i][0]
            hadron2 = remove_heavy[i][1]   
            hadron1mass = self.Getmass(hadron1)
            hadron2mass = self.Getmass(hadron2)
            sumofmasses = hadron1mass + hadron2mass
            totalmass.append(sumofmasses)
        
        minmass = min(totalmass)
        maxmass = max(totalmass)
        
        return minmass, maxmass
        
    def Quarkdiquarkpairmass(self,string):
        Decaytables = Decays()
        stringcontent = string.Getstringcontent()
        quark = stringcontent[0]
        quarkID = PDG_ID.PDG_Type(quark)
        antiquark = stringcontent[1]
        antiquarkID = PDG_ID.PDG_Type(antiquark)
        
        if quarkID.isDiquark() is True: # quark is the diquark, antiquark is a regular quark
            if antiquarkID.isQuark() is False:
                print "Taking illegal arguments, consider revision of quarkdiquarkpairmass function"
                sys.exit()
            decays = Decaytables.Getdecays(float(quark))
            heaviestbaryon = decays[0]
            heaviestbaryonmass = self.Getmass(float(heaviestbaryon))
            vacumnflavour = Decaytables.Getbaryonleftover(heaviestbaryon,quark) # baryonleftover produces abs of particle flavour.
            if quark > 0:
                vacumnflavour = vacumnflavour * -1 # if positive diquark then need the negative vacumn
        
            momen = string.Getqvector4D()
            momen2 = string.Getqbarvector4D()
            pseudoparticle = Particle.particle4D(vacumnflavour,momen)
            otherparticle = Particle.particle4D(antiquark, momen2)
            newstring = MRS.mrs4D(pseudoparticle,otherparticle) # create a string with flavour content to reuse another module
            newstring = self.Reshuffleflavours(newstring) # potentially reshuffle if the orders the wrong way around.
            hadrons = self.GetFailsafe1mesonlist(newstring)
            #print hadrons, "HADRONS"
            #print newstring.Getstringcontent(), "STRING CONTENT GOING INTO MESONHEAVY CODE - quark is diquark"
            heaviestmesoncode = self.Heaviestmesonfinder(hadrons)
            
            #print heaviestmesoncode, "HEAVIESTMESONCODE"
            heaviestmesonmass = self.Getmass(float(heaviestmesoncode))
            totalmax = heaviestmesonmass + heaviestbaryonmass
            totalmin = Constants.minmassmodifier() * totalmax
            return totalmin, totalmax
            
            
           
                
        if antiquarkID.isDiquark() is True: # antiquark is the diquark
            if quarkID.isQuark() is False:
                print "Taking illegal arguments, look at mindiquarkmass function"
                sys.exit()
            decays = Decaytables.Getdecays(float(antiquark))
            heaviestbaryon = decays[0]
            heaviestbaryonmass = self.Getmass(heaviestbaryon)
            vacumnflavour = Decaytables.Getbaryonleftover(heaviestbaryon,antiquark)
            
            if antiquark > 0:
                vacumnflavour = vacumnflavour * -1 # if positive diquark then need the negative vacumn
            
            momen = string.Getqvector4D()
            momen2 = string.Getqbarvector4D()
            pseudoparticle = Particle.particle4D(vacumnflavour,momen)
            otherparticle = Particle.particle4D(quark, momen2)
            newstring = MRS.mrs4D(pseudoparticle,otherparticle)
            newstring = self.Reshuffleflavours(newstring) # potentially reshuffle if the orders the wrong way around.
            hadrons = self.GetFailsafe1mesonlist(newstring)
            #print hadrons, "HADRONS"
            #print newstring.Getstringcontent(), "STRING CONTENT GOING INTO MESONHEAVY CODE - antiquark is diquark"
            heaviestmesoncode = self.Heaviestmesonfinder(hadrons)
            #print heaviestmesoncode, "HEAVIESTMESONCODE"
            heaviestmesonmass = self.Getmass(float(heaviestmesoncode))
            totalmax = heaviestmesonmass + heaviestbaryonmass
            totalmin = Constants.minmassmodifier() * totalmax
            return totalmin, totalmax
          
            
            
    def Heaviestmesonfinder(self,listofmesons):
        if not listofmesons:
            print "EMPTY LIST FED IN TO HEAVIESTMESONFINDER - CONSULT ERROR"
            sys.exit()
            
        if type(listofmesons) is int:
            return listofmesons # already have the heaviest meson. Only one meson possiblwe
            
        masses = []    
        for i in range(len(listofmesons)):
            mass = self.Getmass(listofmesons[i])
            masses.append(mass)
        
        maxmass = max(masses)
        index = masses.index(maxmass)
        code = listofmesons[index]
        
        return code
        
    def Checkforsinglekaon(self,code):
        
        rng = bool(random.getrandbits(1))
        if rng == 0:
            code = 130 
        if rng == 1:
            code = 310  
        return code  
                    
        
    def Checkfordoublekaon(self,hadronpair):
        
        rng = bool(random.getrandbits(1))
        if rng == 0:
            hadronpair[0] = 130
            hadronpair[1] = 310
            
        if rng == 1:
            hadronpair[0] = 310
            hadronpair[1] = 130
            
        return hadronpair
        
    def GetFailsafe1mesonlist(self,string):
        ''' Finds and returns a list of possible meson numbers for a string of q qbar flavour (string collapse) '''
        ''' TO BE USED WITH FAILSAFE 1 IN TWO BODY DECAY - I.E. IT FINDS THE POSSIBLE mesonS TO COLLAPSE INTO '''
        quarkIDs = string.Getstringcontent()
        qID = quarkIDs[0]
        qbarID = quarkIDs[1]
        if qID < 0 and qbarID > 0: # Produced a baryon at some point, got a reversed string order.
            temp = qID
            qID = qbarID
            qbarID = temp # swap the flavours around to find mesons
        possiblemesons = []
        
        
        ''' Mixed string states - these are expressedly hard coded into the system. If adjustments made must go and correct this list! '''
        if abs(qID) == abs(qbarID): # write the exceptions on a case by case for the mixed state mesons. Theyre a mess otherwise.
            if qID == 1: # what states contain a ddbar set of quarks?
                ddbarmesons = [111,113,221,223,331]    
                chosenmeson = self.Findheaviest(string,ddbarmesons)
                return ddbarmesons
                #return chosenmeson    
                
            if qID == 2: # same process but for the uubar quarks
                uubarmesons = [111,113,221,223,331] 
                chosenmeson = self.Findheaviest(string,uubarmesons)
                return uubarmesons
                #return chosenmeson
            
            if qID == 3: # Process for ssbar quarks
                ssbarmesons = [221,331,333]
                chosenmeson = self.Findheaviest(string,ssbarmesons)
                return ssbarmesons
                #return chosenmeson
                
            if qID == 4: # Process for ccbar strings
                ccbarmesons = [441,443]
                chosenmeson = self.Findheaviest(string,ccbarmesons)
                return ccbarmesons
                #return chosenmeson
                
            if qID == 5: # Process for bbbar strings
                bbarmesons = [551,553]
                chosenmeson = self.Findheaviest(string,bbarmesons)
                return bbarmesons
                #return chosenmeson
                
                    
        IDs = PDG_ID.PDG_ID()                    
        mesonlist = IDs.Getmesonlist()
        ''' create negative meson list - i.e. anti mesons '''
        negmesonlist = copy.deepcopy(mesonlist)
        for i in range(len(negmesonlist)):
            negmesonlist[i] = negmesonlist[i] *-1
        totalmesonlist = mesonlist + negmesonlist
        
        ''' Search over the entire list of encoded meson content and find matches. Then append them to a list '''
        for i in range(len(totalmesonlist)):
            if self.__mesoncontent.has_key(totalmesonlist[i]):
                Mesonq = self.__mesoncontent[totalmesonlist[i]].GetqID()
                Mesonqbar = self.__mesoncontent[totalmesonlist[i]].GetqbarID()
                if Mesonq == qID and Mesonqbar == qbarID:
                    possiblemesons.append(totalmesonlist[i])
        
        
        return possiblemesons # used to be possiblemesons   
            
        
                
                            
    def GetFailsafe1hadron(self,string):
        
        if string.Getq().Checkquark() is True and string.Getqbar().Checkquark() is True: # two quarks - meson production
            
            possiblemesons = self.GetFailsafe1mesonlist(string)
            #print possiblemesons, "possiblemesons inside GetFailsafe1hadron"
            chosenmeson = self.Findheaviest(string,possiblemesons) 
            return chosenmeson
            
        if string.Getq().Checkdiquark() is True or string.Getqbar().Checkdiquark() is True: # string has a diquark! need to produce a baryon.
             chosenbaryon = self.Failsafe1baryon(string) 
             return chosenbaryon
        
        
        
        ''' Use the submodule thingy to find the heaviest possible meson from the matching mesons ''' # currently throws error for empty list.
              
        #print "GETFAILSAEF1HADRON CODE HAS FOUND AN EXCEPTION"
        #sys.exit()
            
    def Findheaviest(self,string,hadronlist): # hadronlist here in numerical codes
        stringmass = string.Getstringmass()
        
        if type(hadronlist) is int:
            return hadronlist
        
               
        possiblehadrons = []
        for i in range(len(hadronlist)): # Create the list of mesons with low enough mass
            hadronmass = self.Getmass(hadronlist[i])
            if hadronmass < stringmass:
                possiblehadrons.append(hadronlist[i]) # append the light enough hadrons
                
        if len(possiblehadrons) == 0:
            return 0
            
        
        ''' Find the heaviest of the light hadrons '''
        hadronmasses = []
        for i in range(len(possiblehadrons)):
            hadronmasses.append(self.Getmass(possiblehadrons[i]))
            
        
        
        chosenmass = max(hadronmasses)
        index = hadronmasses.index(chosenmass)
        heaviesthadron = possiblehadrons[index]
        
                
        return heaviesthadron # returns a number


    def Stringdecaysmallesthadron(self,string):
        stringmass = string.Getstringmass()
        string = self.Reshuffleflavours(string)
        hadronlist = self.GetFailsafe1mesonlist(string)
        
        if hadronlist == 0:
            return 0
        if type(hadronlist) is int:
            if self.Getmass(hadronlist) < stringmass:
                return hadronlist
            return 0
            
        if len(hadronlist) == 1:
            if self.Getmass(hadronlist[0]) < stringmass:
                return hadronlist
            return 0
            
        for i in range(len(hadronlist)):
            #print "Hadronlist", hadronlist
            if self.Getmass(hadronlist[i]) < stringmass:
                return 1
        
        '''
        if listofpossible == 0:
            return None
            
        return listofpossible
        #print type(listofpossible), "type of hadrons list"
        #print listofpossible, "listofpossible"
        stringmass = string.Getstringmass()
        
        
        if type(listofpossible) is int:
            if self.Getmass(listofpossible) > stringmass:
                return 1
            return None
            
        '''   
        ''' 
        hadron = self.Findheaviest(string,hadronlist)   
        
        if hadron == 0:
            return None
        '''    
        return 0
            
            
    def Reshuffleflavours(self,string):
        
        temp = string.Getstringcontent()
        
        if string.GetqID() < 0:
            string.SetqID(temp[1])
            string.SetqbarID(temp[0])
            return string
            
        return string
        
    def Failsafe1baryon(self,string):
        ''' Takes a diquark-quark string and finds the heaviest baryon possible to produce from that string configuration '''
        #Baryontables = Baryoncontent.Baryontables()
        Decaytables = Decays()


        
        negative = False
        
        stringcontent = string.Getstringcontent()
        if stringcontent[0] < 0 and stringcontent[1] >0:
            #print stringcontent, "error in the failsafe1baryon moddule, being fed incorrect string identity"
            sys.exit()

        if stringcontent[0] > 0 and stringcontent[1] < 0:
            #print stringcontent, "error in the failsafe1baryon stuff, check the string input"
            sys.exit()

        if stringcontent[0] < 0: # If the content of the string is negative, need to produce an antibaryon.
            negative = True
            
        stringcontent[0] = abs(stringcontent[0])
        stringcontent[1] = abs(stringcontent[1]) # Get the absolute values.
        if string.Getq().Checkdiquark() is True: # The quark entry is the diqurk
            diquark = stringcontent[0]
            quark = stringcontent[1]
            
        if string.Getqbar().Checkdiquark() is True: # antiquark entry is the diquark
            diquark = stringcontent[1]
            quark = stringcontent[0]
        
        
        PDGclass = PDG_ID.PDG_ID()
        Baryonlist = PDGclass.Getbaryonlist() # Get the list of baryons from the code
        Possiblebaryons  = []
        for i in range(len(Baryonlist)): # iterate over the list of baryons.
            Baryoncontent = self.__baryonwfs[Baryonlist[i]].Getwfs() # get the current baryon under consideration. returns [diquark, wf, leftover]
            for n in range(len(Baryoncontent)):
                baryon_diquark = Baryoncontent[n][0]
                baryon_quark = Baryoncontent[n][2]
                if baryon_diquark == diquark and baryon_quark == quark:
                    Possiblebaryons.append(Baryonlist[i])
            
        if not Possiblebaryons: #empty list. No baryons found. Have a look into what error we have.
            return 0 # added return 0. does this work?
            #print "no baryon found. Input here incorrect."
            #sys.exit()
            
        if len(Possiblebaryons) > 1: # Multiple possible combinations. Choose the heaveist baryon as most kinematically preferable
            baryonmasses = []
            #print Possiblebaryons, "POSSIBLE BARYONS"
            for i in range(len(Possiblebaryons)):   
                baryonmass = self.Getmass(Possiblebaryons[i])
                baryonmasses.append(baryonmass)
            
            #print baryonmasses, "BARYONMASSES"
                
            heaviest = max(baryonmasses)
            index = baryonmasses.index(heaviest)
            chosencode = Possiblebaryons[i]
            
        if len(Possiblebaryons) == 1:
            chosencode = Possiblebaryons[0]
              
            
        
        if negative == True:
            chosencode = chosencode * -1
            
                        
        return chosencode            
        
        
        
        
        

Decaytables = Decays()





        
''' Testing the baryon code '''
'''
value = Decaytables.Getbaryonwf(2224,2203)
#print value

'''














