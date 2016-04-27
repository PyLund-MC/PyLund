''' Available decays class - contains the algorithm etc for finding the meson content of a string decay '''

import numpy
import Boosts
import Fourvector
import Particle
import PDG_ID
import random
import MRS
import sys
import Constants
import Stringdecay
import Decaydictionary

#Decaytables = Decaydictionary.Decays() # load in the decay tables locally

class availabledecays4D(object):
    def __init__(self,string):
        
        self.__mrs             = string
        self.__stringmass      = string.Getstringmass()
        self.__totalvector4D   = string.Gettotalvector4D()
        self.__totalmom        = string.Gettotalvector()
        
        self.__q                = string.Getq()
        self.__qID              = str(self.__q.Getcode())
        self.__qmom             = string.Getqvector()
        self.__qvector4D        = string.Getqvector4D()
        self.__qcharge          = string.Getqcharge()
        
        self.__qbar             = string.Getqbar()
        self.__qbarID           = str(self.__qbar.Getcode()) # removed the abs()
        self.__qbarmom          = string.Getqbarvector()
        self.__qbarvector4D     = string.Getqbarvector4D()
        self.__qbarcharge       = string.Getqbarcharge()
        
        self.__totalcharge      = self.__qbarcharge + self.__qcharge
        self.__flavour          = ''.join([self.__qID,self.__qbarID])
        
    def Getend(self):
        rng = random.getrandbits(1) # generate a number which is 0 or 1. Line here in case of method change to speed up code.
        if rng == 0:
            return self.__qID, self.__qbarID
        if rng == 1:
            return self.__qbarID, self.__qID
        
    def Getflavour(self):
        return self.__flavour
        
    def Finddecay(self):
        sdc = Stringdecay.stringdecay4D(self.__mrs) # stringdecayclass
        content = self.Getend()
        chosenquark = content[0] # this is either the quark/antiquark of the string QID, QbarID.
        
        chosenquark = float(chosenquark)
        chosencode = PDG_ID.PDG_Type(chosenquark)
        otherendcode = PDG_ID.PDG_Type(float(content[1]))
        
        if chosencode.isDiquark() is True: # If decaying end is a diquark, doesnt actually matter what the other end is. The diquark does its own thing and leaves a remnant.
            return self.Quarkdecay(chosenquark)
            
        if chosencode.isQuark() is True and otherendcode.isDiquark() is True: # In the event that you have a quark and a diquark, and you choose the quark end.
            # The chosen quark cannot produce a diquark as there is already one at the other end. 
            while True:
                product = self.Quarkdecay(chosenquark)
                productcode = PDG_ID.PDG_Type(float(product[1]))
                #print "AVAIALBLEDECAYS STUCK PLACE 1"
                #print productcode
                #print productcode.isDiquark(), "CODE TRUTH"
                if productcode.isDiquark() is False:
                    return product
                continue
        # In the event that we have two quarks, this part runs. Need to make sure that we conserve spin of diquarks and dont produce something unphysical (doesnt have to make a diquark.
        invaliddiquarks = [1101,2201,3301]
        while True:
            product = self.Quarkdecay(chosenquark)
            productcode = PDG_ID.PDG_Type(float(product[1]))
            #print "AVAILABLEDECAYS STUCK PLACE 2"
            if productcode.isDiquark() is True:
                #print "AVAILDECAYS STUCK 2 LOOP"
                if self.__mrs.Getstringmass() <= Constants.baryonproductionthreshold():
                    "AVAILDECAYS STUCK 2 LOOP BBB"
                    continue
                    
                #print "PASSED CRITERIONS"
                recievingdiquark = sdc.Finddiquarkstringcontent(chosenquark,product[1]) # check that you dont produce an invalid diquark (spin conservation)
                #print "GOT PAST DIQUARK FORMER"
                
                if abs(recievingdiquark) in invaliddiquarks:
                    #print "INVALID DIQUARKS BREAKING"
                    continue
            return product
            
            
        return self.Quarkdecay(chosenquark)
            
            
                
        
        
    def Checkforkaons(self,code):
        if abs(code) == 311:
            rng = random.getrandbits(1) # generate a number which is 0 or 1. Line here in case of method change to speed up code.
            if rng == 0:
                return 130
            if rng == 1:
                return 310
                
        return code
   
        
        
    def Quarkdecay(self,chosenstringend):
        Decaytables = Decaydictionary.Decays()
        stringmass = self.__stringmass    # criterion for kinematically possible decay is set variable in the constants file
        #print chosenstringend, "CHOSENSTRING END"
        decays      = Decaytables.Getdecays(chosenstringend)
        #print decays, "POSSIBLE DECAYS FOR THAT STRING END, CHECKIN THE CONSERVATION"
        decaymasses = Decaytables.Getmasses(chosenstringend)
        probs       = Decaytables.Getfinalprobs(chosenstringend)
        minimumstringmass = stringmass 
        allowedmasses = []
        for i in range(len(decaymasses)):
            if decaymasses[i] < minimumstringmass:
                allowedmasses.append(decaymasses[i])
        no_decays = len(allowedmasses) # find the len of the list (how many decays are there)
        if no_decays == 0:
            return [chosenstringend, decays[-1]]
        allowed_probs = probs[-no_decays:] # slice from the right hand side. find the correct probabilities
        #print allowed_probs, "PROBS HERE"
        sumprobs = sum(allowed_probs)
        allowed_probs       = allowed_probs/sumprobs            # normalise the probs to 1 - theyre still correctly proportioned.
        
        allowed_decays = decays[-no_decays:]
        rng = numpy.random.choice(no_decays,1,p=allowed_probs) # generate a number between 0 - N based on the given probabilities.
        chosenhadronID = allowed_decays[rng] # find the corresponding decay mesonID
        chosenhadronID = float(chosenhadronID)
        return [chosenstringend, chosenhadronID]
        
                
                                
                
'''        
decaymasses = (5,4,3,2,1)
probs = (0.1,0.2,0.3)
rng = numpy.random.choice(3,1,probs)
print rng

threshold = 3
allowed = [i for i in decaymasses if i <= threshold]
number = len(allowed)
print allowed
print decaymasses[-number:]
'''