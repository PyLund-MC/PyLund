from __future__ import division
import numpy
import math
import Constants



class vector4D(object):
    ''' Class written by Ryan, defines 4 vector algebra '''
    def __init__(self,t_e=0,x_px=0,y_py=0,z_pz=0): # Type has been removed
        self.__vector = [t_e, x_px, y_py, z_pz]
        
    def Printself(self):
            print "(t_e = ",self.__vector[0], ", x_px = ",self.__vector[1], ", y_py = ",self.__vector[2],
            print ", z_pz = ",self.__vector[3], ")"
    
    def __str__(self):
        return str(self.__vector)
    
    def Getvector(self):
        return self.__vector
    
    def __getitem__(self,idx):
        assert idx<4 and idx>=0
        return self.__vector[idx]

    def __setitem__(self,idx,value):
        assert idx<4 and idx>=0
        self.__vector[idx] = value
        
    def __add__(vector1,vector2):
        vec = vector4D()
        for i in range(4):
            vec[i] = vector1[i] + vector2[i]
        return vec
        
    def __sub__(vector1,vector2):
        vec = vector4D()
        for i in range(4):
            vec[i] = vector1[i] - vector2[i]
        return vec
        
    def __mul__(vector1,vector2):
        scalar = vector1[0] * vector2[0]
        for i in range(1,4):
            scalar = scalar - (vector1[i] * vector2[i])
        if scalar < Constants.Numerical_limit():
            scalar = 0
        return scalar
        
    def __iadd__(self,vector):
        for i in range(4):
            self[i] = self.__vector[i] + vector[i]
        return self
        
    def __isub__(self,vector):
        for i in range(4):
            self[i] = self.__vector[i] - vector[i]
        return self
        
    def __imul__(self,number):
        for i in range(4):
            self[i] = self.__vector[i] * number
        return self
        
    def __idiv__(self,number):
        for i in range(4):
            self[i] = self.__vector[i] / number
        return self
        
    def Checknumerical(self):
        numerical_limit = Constants.Numerical_limit()
        for i in range(4):
            if abs(self.__vector[i])<=numerical_limit:     # Set the limit for numerical accuracy
                self.__vector[i] = 0    # Number is too small? - set it to zero.
            else:
                pass
        return self    
            
    def Addvectors(vector1,vector2):
        vec = vector4D()
        for i in range(4):
            vec.__vector[i] =  vector1[i] + vector2[i]
        return vec 
        
    def Subvectors(vector1,vector2):
        vec = vector4D()
        for i in range(4):
            vec.__vector[i] = vector1[i] - vector2[i]
        return vec        
        
    def Dotproduct(vector1,vector2):
        numerical_limit = Constants.Numerical_limit()   
        scalar = vector1[0] * vector2[0]
        for i in range(1,4):
            scalar = scalar - vector1[i] * vector2[i]
        if scalar <= numerical_limit:
            scalar = 0
        return scalar
    
    def Getorientation(self):
        mag =  numpy.sqrt(self.__vector[1]**2 + self.__vector[2]**2 + self.__vector[3]**2)
        theta = numpy.arccos(self.__vector[3] / mag)
        phi = math.atan2(self.__vector[2],self.__vector[1]) # could only boost negative values for range [0,2pi]
        if phi < 0:   
            phi += 2*numpy.pi
        return [mag, theta, phi]
        
    def Getseperation(vector1,vector2):
        sep = vector4D.Dotproduct(vector1,vector1) + vector4D.Dotproduct(vector2,vector2) - 2*vector4D.Dotproduct(vector1,vector2)
        return numpy.sqrt(sep)
        
    def Copy(self):
        vec = vector4D()
        for i in range(4):
            vec[i] = self.__vector[i]
        return vec
        
''' class additions? '''
 
# Getlike? gets spacelike, timelike or lightlike type of a vector. i.e. test s**2           
            
            
        