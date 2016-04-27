''' Full decay class - decays a given string a number of times '''

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
import Twobodydecay
import Hadrondecay
#Decaytables = Decaydictionary.Decays()

class fulldecay4D(object):
    def __init__(self):
        #assert isinstance(mrs,MRS.mrs4D) # throwing an error for now, have a look at in the future
        self.__classname        = "I'm the string decay class!"
        self.__description      = "I contain several functions, used for the boosting of strings, their decay, and the associated kinematics!"
        
    
    def Performmanydecays(self,initialstring, No_runs): # performs a decay of a single input string, N runs through
        listofhadronlists = []
        # Create empty lists to append the hadrons produced by the string decaying 
        for i in range(No_runs):
            decayedhadrons = self.Performdecaygeneral(initialstring)
            
            listofhadronlists.append(decayedhadrons)
            #print "PERFORM MANY DECAYS APPEND SUCCESSFUL"
        return listofhadronlists
                    

    def Performdecaygeneral(self,initialstring): # decays a string into hadronic state. Does not decay hadrons further.
        
        ''' Put the initialisation here. Things that dont get altered during the iterative loop '''
        Decaytables = Decaydictionary.Decays()
        twobodyclass = Twobodydecay.twobodydecay4D()
        decayingstring = initialstring
        no_decays = 0
        string2 = copy.deepcopy(decayingstring) # string 2 is the string that has all the editing done to it.
        strings =[]
        strings.append(string2)
        hadrons = []
        
        
        while True: # First loop uses the initialed values. Values are altered at the end of each loop.
            
            stringmass = string2.Getstringmass()
            #print stringmass, "stringmass"
            minpairmass, maxpairmass = Decaytables.Findminandmaxpairmass(string2)
            #print minpairmass, maxpairmass, "min and max mass"
            
            
            ''' Check for Failsafe 1 mass criterion, If met, perform Failsafe 1 and end the decay '''
            
            if stringmass <= minpairmass:
                #print "FAILSAFE 1 ENGAGED"
                Failsafe1state = twobodyclass.Failsafe1(string2) # produces a hadron and a photon (string collapse) from the string.
                hadrons += Failsafe1state
                if Failsafe1state[0].Getcode() == 311:
                    #print "FOUND THE FUCKER"
                    sys.exit()
                
                
                return hadrons
                
            ''' Check for Failsafe 2 mass criterion, If met perform Failsafe 2 and end the decay '''
            if stringmass <= maxpairmass:
                #print "FAILSAFE 2 ENGAGED"
                Failsafe2state = twobodyclass.Failsafe2(string2) # produces two final hadrons from string (seperation)
                if abs(Failsafe2state[0].Getcode()) == 311:
                    #print "FOUND THE FUCKSHIT"
                    sys.exit()
                if abs(Failsafe2state[1].Getcode()) == 311:
                    #print "DOOT DE FRUIT"
                    sys.exit()
                hadrons += Failsafe2state
                return hadrons
                
            ''' If string sufficiently heavy enough then, proceeds to decay via hadronic emission - "regular" string decay '''
            #print "REACHED PAST THE FAILSAFE CODES"    
            decayedsystem = self.Decaystringonce(string2) # produces [newstring, newhadron]
            string2 = decayedsystem[0]
            producedhadron = decayedsystem[1]
            #print producedhadron, "prodhad"
            if abs(producedhadron.Getcode()) == 311:
                    #print "FOUND THE FUCKER"
                    sys.exit()
            hadrons.append(producedhadron) # add the hadron to the hadron list
            #print "FULL STRINGDECAY MODULE COMPLETED"
            
                        
    def Findsmallestmass(self,decayingquarkno):
        Decaytables = Decaydictionary.Decays()
        return Decaytables.Getmindecay(decayingquarkno)                                                 


    def Sortlists(self,listofstringlists,nodecays): 
        ''' finds string lists of length no decays + 1 '''
        chosendecays = []
        for i in range(len(listofstringlists)):
            if len(listofstringlists[i]) == (nodecays + 1): # nodecays +1 - i.e. the  initial string as well as the decayed strings
                chosendecays.append(listofstringlists[i])
      
        return chosendecays
        
        
    def Findmasses(self,stringlists):
        ''' takes a list of string lists of known (or unknown) length and returns the masses '''
        newstrings = copy.deepcopy(stringlists)
        for i in range(len(newstrings)):
            for p in range(len(newstrings[i])):
                newstrings[i][p] = newstrings[i][p].Getstringmass()
                
        return newstrings
        
        
    def Findaverages(self,stringlists):
        ''' Takes a set of strings, all of the same length (number of decays) and returns the average masses'''
        newstring = copy.deepcopy(stringlists)
        for i in range(len(newstring)):
            for p in range(len(newstring[i])):
                newstring[i][p] = newstring[i][p].Getstringmass()
                
        masses = []        
        listofaverages = []
        if not newstring:
            return []       
        for p in range(len(newstring[i])):
            for i in range(len(newstring)):
                masses.append(newstring[i][p])
            listofaverages.append(numpy.mean(masses))
            masses = []
                    
        return listofaverages
            
            
    def Gethadronfrequencies(self,listofhadronlists,hadroncode):
       ''' Takes the product hadron lists and a hadron code and tells you how many of that hadron were produced '''
       nohadron = 0.
       for i in range(len(listofhadronlists)):
           for p in range(len(listofhadronlists[i])):
               if listofhadronlists[i][p].Getcode() == hadroncode:
                   nohadron += 1
         
       return float(nohadron)
        
    def Gettotalhadronno(self,listofhadronlists):
        totalnohadrons = 0.
        for i in range(len(listofhadronlists)):
            for p in range(len(listofhadronlists[i])):
                totalnohadrons += 1
                
        return float(totalnohadrons)

    def Decaystringonce(self,string):
        
        Availabledecayclass = Availabledecays.availabledecays4D(string) # create available decays class
        decayinfo = Availabledecayclass.Finddecay() # find decayinfo for the decay - [decaying quark, hadronID]
        decayingquark = decayinfo[0]
        hadron_ID = decayinfo[1]
        #print hadron_ID, "DECAYSTRINGONCE HADRON ID"
        
        stringdecayclass2 = Stringdecay.stringdecay4D(string)  #create the stringdecayclass
                    
        decayedsystem2 = stringdecayclass2.Decaystring(string,decayingquark,hadron_ID) # perform the decay - getting [decayedstring, producedmeson]
        
        newstring = decayedsystem2[0]
        producedhadron = decayedsystem2[1]
        
        
        #print producedhadron.Getcode()
        #print "DECAYSTRINGONCE PRINTOUT"
        
        return [newstring,producedhadron]
        
    def Performdecayfull(self,initialstring):
        hadrondecayclass = Hadrondecay.hadrondecay4D()
        
        poststringhadrons = self.Performdecaygeneral(initialstring)
        posthadrondecaylist = hadrondecayclass.Decayparticlelist(poststringhadrons)
        
        return posthadrondecaylist
        
        
    def Performdecayfullmanytimes(self,initialstring,No_runs):
        hadrondecayclass = Hadrondecay.hadrondecay4D()
        
        listofhadronlists = []
        finalhadronlist = []
        # Create empty lists to append the hadrons produced by the string decaying 
        for i in range(No_runs):
            decayedhadrons = self.Performdecaygeneral(initialstring)
            
            listofhadronlists.append(decayedhadrons) # initial undecayed hadrons
            
        for i in range(len(listofhadronlists)):
            finalstate = hadrondecayclass.Decayparticlelist(listofhadronlists[i])
            finalhadronlist.append(finalstate)
            
        
        
        return finalhadronlist
        
        
    def Gethadronpercentages(self,listofhadronlists,produciblehadrons): # takes multiple list arguments, as well as the producible hadrons that you want the percentages of
        
        hadronpercents = []
        totalhadronno = self.Gettotalhadronno(listofhadronlists)
        
        for i in range(len(produciblehadrons)):
            freq = self.Gethadronfrequencies(listofhadronlists,produciblehadrons[i])
            percentage = freq*100/totalhadronno
            hadronpercents.append(percentage)
            
        return hadronpercents
            
    def Gethadronabsoluteabundance(self,listofhadronlists,produciblehadrons):
        hadronabundance = []
        totalnoevents = len(listofhadronlists)
        
        for i in range(len(produciblehadrons)):
            freq = self.Gethadronfrequencies(listofhadronlists,produciblehadrons[i])
            howmany = freq/totalnoevents
            hadronabundance.append(howmany)
            
        return hadronabundance
            
        
    
    def Gethadronnames(self,produciblehadrons):
        
        hadronnames = []
        
        for i in range(len(produciblehadrons)):
            classcode = PDG_ID.PDG_Type(produciblehadrons[i])
            name = classcode.name()
            
            hadronnames.append(name)
            
        return hadronnames
            
        
    