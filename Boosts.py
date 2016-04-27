''' Boost class - written by Ryan '''

from __future__ import division
import numpy
import math
import Fourvector
import cmath
import math
import MRS
#import Twobodydecay

class boost4D(object):
    ''' Takes 4 momentum 4 vector as an argument '''
    def __init__(self,mom=Fourvector.vector4D()):
        #print "boost class called"
        self.__mom = mom
        self.__beta = Fourvector.vector4D.Getvector(mom)
        self.__mass = numpy.sqrt(Fourvector.vector4D.Dotproduct(mom,mom)) # dot product w/ self - sqrt # TOOK AWAY THE SQRT
        
    def Printself(self):
        print "Input: (E = ",self.__beta[0], ", px = ",self.__beta[1], ", py = ",self.__beta[2]
        print ", pz = ",self.__beta[3], ")"
       # print "( mass = ", self.__mass, "units", "vx =", self.__betax, "c", "vy =", self.__betay, "c", "vz = ", self.__betaz, "c )"
        #print "( gammax =", self.__gammax, "gammay =", self.__gammay, "gammaz =", self.__gammaz, ")"
        
    
    def GetCMF(self):
        return Fourvector.vector4D(self.__mass,0,0,0)
           
    def __mul__(self,vector):   # create the forward boost operator
        selfvec = Fourvector.vector4D()
        for i in range(4):
            selfvec[i] = self.__beta[i]
        newvec0 = Fourvector.vector4D.Dotproduct(selfvec,vector)/self.__mass
        angle = (vector[0]+newvec0)/(self.__mass+selfvec[0])
        newvec = Fourvector.vector4D()
        for i in range(4):
            newvec[i] = vector[i]
        newvec[0] = newvec0
        for i in range(1,4):
            newvec[i] -= angle*selfvec[i]
            if abs(newvec[i]/newvec0)<1.e-14:
                newvec[i] = 0
        newvec.Checknumerical()
        return newvec

    def __div__(self,vector):   # create the inverse boost operator
        selfvec = Fourvector.vector4D()
        for i in range(4):
            selfvec[i] = self.__beta[i]   # initise the boost vector
        newvec0 = 0
        for i in range(4):
            newvec0 += (self.__beta[i]*vector[i]) 
        newvec0 = newvec0 / self.__mass 
        angle = (vector[0]+newvec0)/(self.__mass+selfvec[0])
        newvec = Fourvector.vector4D()
        for i in range(4):
            newvec[i] = vector[i]
        newvec[0] = newvec0
        for i in range(1,4):
            newvec[i] += (angle*selfvec[i])
            if abs(newvec[i]/newvec0)<1.e-10: 
                newvec[i] = 0  # 
        newvec.Checknumerical()
        return newvec 
        
                
        
class rotate4D(object):
    
    def __init__(self,mom=Fourvector.vector4D()):
        self.__mom = mom
        self.__beta = Fourvector.vector4D.Getvector(mom)
        self.__values = Fourvector.vector4D.Getorientation(mom)
        self.__theta = self.__values[1]
        self.__phi = self.__values[2]
        self.__costheta = numpy.cos(self.__theta)
        self.__sintheta = numpy.sin(self.__theta)
        self.__cosphi = numpy.cos(self.__phi)
        self.__sinphi = numpy.sin(self.__phi)
        self.__mass = numpy.sqrt(Fourvector.vector4D.Dotproduct(mom,mom)) 
        
    def Getangles(self):
        return "theta =", self.__theta, "phi =", self.__phi
    
                                                     
    
    def Rotateforwards(self,vector): # takes a vector4D as an argument?
        
         #create the rotation matrix 
        M = numpy.zeros((3,3))
        M[0,0] = self.__cosphi*self.__costheta
        M[0,1] = self.__sinphi*self.__costheta
        M[0,2] = -self.__sintheta
        M[1,0] = -self.__sinphi
        M[1,1] = self.__cosphi
        M[1,2] = 0
        M[2,0] = self.__cosphi*self.__sintheta
        M[2,1] = self.__sinphi*self.__sintheta
        M[2,2] = self.__costheta
        
        vec = vector.Getvector()[1:4] # self.___beta
        boostedvec = numpy.dot(M,vec)
        finalvec = Fourvector.vector4D()
        finalvec[0] = vector.Getvector()[0] # self.__beta
        for i in range(1,4):
            finalvec[i] = boostedvec[i-1]
        finalvec.Checknumerical()    
        return finalvec
        
    def Rotatebackwards(self,vector):
        
        M = numpy.zeros((3,3))
        M[0,0] = self.__cosphi*self.__costheta
        M[0,1] = -self.__sinphi
        M[0,2] = self.__cosphi*self.__sintheta
        M[1,0] = self.__sinphi*self.__costheta 
        M[1,1] = self.__cosphi
        M[1,2] = self.__sinphi*self.__sintheta
        M[2,0] = -self.__sintheta
        M[2,1] = 0
        M[2,2] = self.__costheta
        
        vec = vector.Getvector()
        vec = vec[1:4]
        boostedvec = numpy.dot(M,vec)
        
        finalvec = Fourvector.vector4D()
        finalvec[0] = vector[0]
        for i in range(1,4):
            finalvec[i] = boostedvec[i-1]
        finalvec.Checknumerical()
        return finalvec
        
        
