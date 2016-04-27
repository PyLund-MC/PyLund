from __future__ import division
import numpy
import Fourvector
import PDG_ID
import Constants
import copy
# Need to import my other modules as well


class particle4D(object):
    def __init__(self,pdgID,mom=Fourvector.vector4D()): # do i need status / id counting the particles seems useful
        #print "particle class called"
        self.__ID  = pdgID
        self.__PDG_ID = PDG_ID.PDG_Type(pdgID)
        self.__vector4D = mom # this used to be jsut mom
        self.__mom = self.__vector4D.Getvector()
        self.__name = PDG_ID.PDG_Type.name(self.__PDG_ID)
        self.__flavour = str(abs(pdgID))
        self.__invariantmass = numpy.sqrt(Fourvector.vector4D.Dotproduct(self.__vector4D,self.__vector4D)) # invariant masssquared
        self.__particlemass = numpy.sqrt(Fourvector.vector4D.Dotproduct(self.__vector4D,self.__vector4D))
        self.__colour = 0
        self.__anticolour = 0 # Colour is by default zero.
        
        
    def Printself(self):
        print "(", self.__name, "mom =", self.__mom, ")"

    def Setcolour(self,colour):
        self.__colour = colour

    def Getcolour(self):
        return self.__colour

    def Setanticolour(self,anticolour):
        self.__anticolour = anticolour

    def Getanticolour(self):
        return self.__anticolour
        
        
    def Getname(self):
        return self.__name
        
    def GetID(self):
        return self.__ID
        
    def Getinvariantmass(self):
        return self.__invariantmass
        
    def Getparticlemass(self):
        return self.__particlemass
        
    def Getcode(self):
        return PDG_ID.PDG_Type.code(self.__PDG_ID)
        
    def Getmass(self):
        return PDG_ID.PDG_Type.mass(self.__PDG_ID)
        
    def Getwidth(self):
        return PDG_ID.PDG_Type.width(self.__PDG_ID)
        
    def Getlifetime(self):
        width = copy.deepcopy(self.Getwidth())
        return Constants.hbar_gevs() / width
        
    def Getvector4D(self):
        return self.__vector4D
        
    def Getcharge(self):
        return PDG_ID.PDG_Type.charge(self.__PDG_ID)
        
    def Getvector(self):
        return self.__mom
        
    def Getmom(self):
        return self.__mom
        
    def Getwidth(self):
        return PDG_ID.PDG_Type.width(self.__PDG_ID)
            
    def Checkfermion(self):
        return PDG_ID.PDG_Type.isFermion(self.__PDG_ID)
    
    def Checklepton(self):
        return PDG_ID.PDG_Type.isLepton(self.__PDG_ID)     
               
    def Checkhadron(self):
        return PDG_ID.PDG_Type.isHadron(self.__PDG_ID)
           
    def Checkquark(self):
        return PDG_ID.PDG_Type.isQuark(self.__PDG_ID) 
        
    def Checkmeson(self):
        return PDG_ID.PDG_Type.isMeson(self.__PDG_ID)
        
    def Checkbaryon(self):
        return PDG_ID.PDG_Type.isBaryon(self.__PDG_ID)
        
    def Checkdiquark(self):
        return PDG_ID.PDG_Type.isDiquark(self.__PDG_ID)

    def Checkgluon(self):
        return PDG_ID.PDG_Type.isGluon(self.__PDG_ID)

    def Getflavour(self):
        if self.__PDG_ID.isQuark() is True:
            return self.__flavour
        if self.__PDG_ID.isMeson() is True:
            return self.__flavour[-3:-1]
        if self.__PDG_ID.isBaryon() is True:
            return self.__flavour[-4:-1]
        if self.__PDG_ID.isDiquark() is True:
            return self.__flavour[-4:-2]
        else:
            return self.__name, "Particle does not have flavour"    
        #if self.__PDG_ID.isLepton is True:
          #  return "ERROR: Particle does not have flavour"


    def Checkanti(self):
        return PDG_ID.PDG_Type.isAnti(self.__PDG_ID)
        
    def Setmom(self,newmom = Fourvector.vector4D()):
        self.__vector4D = newmom
        self.__mom = newmom.Getvector()
        
    def SetID(self,newID):
        self.__ID = newID # added this line previous, make sure code dont shit the bed.
        self.__PDG_ID = PDG_ID.PDG_Type(newID)

        
        
