''' Massless relativistic string (MRS), written by Ryan '''

from __future__ import division
import numpy
import Fourvector
import Particle
import PDG_ID
import inspect

class mrs4D(object):
    
    def __init__(self, quark, quarkbar):    # need position vector?
        #assert isinstance(q, Particle.particle4D)
        #assert isinstance(qbar, Particle.particle4D)
        
        self.__q = quark
        self.__qname = self.__q.Getname()
        self.__qID = self.__q.Getcode()
        self.__qPDG = PDG_ID.PDG_Type(self.__qID)
        self.__qmass = self.__q.Getmass()
        self.__qcharge = self.__q.Getcharge()
        self.__qvector4D = quark.Getvector4D()
        self.__qmom = self.__q.Getvector()
        
        self.__qbar = quarkbar
        self.__qbarname = self.__qbar.Getname()
        self.__qbarID = self.__qbar.Getcode()
        self.__qbarPDG = PDG_ID.PDG_Type(self.__qbarID)
        self.__qbarmass = self.__qbar.Getmass()
        self.__qbarcharge = self.__qbar.Getcharge()
        self.__qbarvector4D = quarkbar.Getvector4D()
        self.__qbarmom = self.__qbar.Getvector()
        
        self.__totalvector4D = self.__qvector4D + self.__qbarvector4D
        #print self.__totalvector4D.Getvector()
        self.__totalmom = Fourvector.vector4D.Getvector(self.__totalvector4D)
        self.__stringmass = numpy.sqrt(Fourvector.vector4D.Dotproduct(self.__totalvector4D,self.__totalvector4D))# functions should return the sqrt(M^2)
        self.__transversemass = numpy.sqrt(self.__totalmom[0]**2 - self.__totalmom[3]**2)        
        
     
    def Getself(self):
        return self
        
    def Copy(self):
        newstring = mrs4D(self.__q, self.__qbar)
        return newstring
        
    def Getstringmass(self):
        return self.__stringmass
        
    def Gettransversemass(self):
        return self.__transversemass
        
    def Getq(self):
        return self.__q
        
    def Getqname(self):
        return self.__qname
        
    def GetqID(self):
        return self.__qID
        
    def SetqID(self,ID):
        self.__qID = ID
        
    def GetqPDG(self):
        return self.__qPDG
        
    def Getqmass(self):
        return self.__qmass
        
    def Getqcharge(self):
        return self.__qcharge
        
    def Getqvector4D(self):
        return self.__qvector4D
        
    def Getqvector(self):
        return self.__qmom
        
    def Getqbar(self):
        return self.__qbar
        
    def Getqbarname(self):
        return self.__qbarname
        
    def GetqbarID(self):
        return self.__qbarID
        
    def SetqbarID(self,ID):
        self.__qbarID = ID
        
    def GetqbarPDG(self):
        return self.__qbarPDG
        
    def Getqbarmass(self):
        return self.__qbarmass
        
    def Getqbarcharge(self):
        return self.__qbarcharge
        
    def Getqbarvector4D(self):
        return self.__qbarvector4D
        
    def Getqbarvector(self):
        return self.__qbarmom
        
    def Gettotalvector(self):
        return self.__totalmom
        
    def Gettotalvector4D(self):
        return self.__totalvector4D
        
    def Getstringcontent(self):
        return [self.__qID, self.__qbarID]
        
    def Printself(self):
        #print "Hello, I'm the Mrs (string) Class"
        print "===== PRINT QUARK DETAILS ====="
        print "q name:", self.__qname, "|", "q PDG_ID:", self.__qID, "|", "q mass:", self.__qmass, "|", 
        print "q mom:", self.__qmom
        print " "
        print "===== PRINT ANTIQUARK DETAILS ====="
        print "qbar name:", self.__qbarname, "|", "qbar PDG_ID:", self.__qID, "|", "qbar mass:", self.__qmass, "|", 
        print "qbar mom:", self.__qbarmom
        
           
                 
    def getEnergy(self):
        return self                   