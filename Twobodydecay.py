''' Class for the decay of a particle into two new particles '''

import numpy
import Boosts
import Fourvector
import Particle
import PDG_ID
import random
import MRS
import Decaydictionary
import Stringdecay
import sys
import copy
import math
#Decaytables = Decaydictionary.Decays()


''' Module which contains the final state decays of a string system - e.g. the standard decays and also the failsafes for insufficient string mass '''

class twobodydecay4D(object):
    def __init__(self):
        self.__classname = "twobodydecay4D"
        
        
    def Getclassname(self):
        return self.__classname
        
    def TwobodyCMF(self,decayingvector,mass1,mass2): # masses fed in are the masses of the products? decay in cmf frame of the particle?
        boost = Boosts.boost4D(decayingvector)
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0]
        mom3 = Mass/2
        mom3 *= self.Sqrtlambda1(Mass,mass1,mass2)
        randomphi = random.random()
        randomtheta = random.uniform(-1,1)
        phi = 2.*numpy.pi*randomphi
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        theta = numpy.pi*randomtheta
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        mom = Fourvector.vector4D(0.,sintheta*cosphi,sintheta*sinphi,costheta)
        mom *= mom3
        E1 = numpy.sqrt(mass1*mass1+mom3*mom3)
        E2 = numpy.sqrt(mass2*mass2+mom3*mom3)
        finalvec1 = Fourvector.vector4D()
        finalvec1[0] = E1
        finalvec1 += mom
        finalvec2 = Fourvector.vector4D()
        finalvec2[0] = E2
        finalvec2 -= mom
        finalvec1 = boost/finalvec1
        finalvec2 = boost/finalvec2
        return [finalvec1, finalvec2]
            
        
    def Failsafe2(self,string): # IN SITU. NEEDS NEW TABLES TO BE CODED.
        stringvec4D = string.Gettotalvector4D()
        boost = Boosts.boost4D(stringvec4D)
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0] # string mass
        
        ''' Get the hadron masses '''
        Decaytables = Decaydictionary.Decays()
        hadronpair = Decaytables.Failsafe2complete(string)
        
        hadron1code = hadronpair[0]
        hadron2code = hadronpair[1]
        mass1 = Decaytables.Getmass(hadron1code)
        if mass1 == 0:
            print "mass 1 == 0"
            print hadron1code, "hadron1code"
            print hadron2code, "hadron2code"
            sys.exit()
        mass2 = Decaytables.Getmass(hadron2code)
        
        mom3 = Mass/2
        mom3 *= self.Sqrtlambda1(Mass,mass1,mass2)
        randomphi = random.random()
        randomtheta = random.uniform(-1,1)
        phi = 2.*numpy.pi*randomphi
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        theta = numpy.pi*randomtheta
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        mom = Fourvector.vector4D(0.,sintheta*cosphi,sintheta*sinphi,costheta)
        mom *= mom3
        E1 = numpy.sqrt(mass1*mass1+mom3*mom3)
        E2 = numpy.sqrt(mass2*mass2+mom3*mom3)
        finalvec1 = Fourvector.vector4D()
        finalvec1[0] = E1
        finalvec1 += mom
        finalvec2 = Fourvector.vector4D()
        finalvec2[0] = E2
        finalvec2 -= mom
        finalvec1 = boost/finalvec1
        finalvec2 = boost/finalvec2
        
        finalvec1row = finalvec1.Getvector()
        finalvec2row = finalvec2.Getvector()
        for i in range(4):
            if math.isnan(finalvec1row[i]) == True:
                finalvec1 = Fourvector.vector4D(mass1,0,0,0)


            if math.isnan(finalvec2row[i]) == True:
                finalvec2 = Fourvector.vector4D(mass2,0,0,0) # if we have nan values, then produce them at reast. see how this works
                print Mass, "stringmass"
                print "reached the problem with double kaon produictiion"
                sys.exit()


        
        finalhadron1 = Particle.particle4D(hadron1code,finalvec1)
        finalhadron2 = Particle.particle4D(hadron2code,finalvec2)
        
        return [finalhadron1, finalhadron2]
        
        
    def Twobodydiquarks(self,string,diquark1, diquark2):  
        
        stringvec4D = string.Gettotalvector4D()
        boost = Boosts.boost4D(stringvec4D)
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0] # string mass
        
        
        
        
        
        ''' Get the hadron masses '''
        Decaytables = Decaydictionary.Decays()
        #hadronpair = Decaytables.Failsafe2complete(string)
        #need to get the diquarks for string ends.
        hadron1code = diquark1
        #print hadron1code, "hadron1code from diquark twobody"
        hadron2code = diquark2
        #print hadron2code, "hadron2code from diquark twobody"
        mass1 = Decaytables.Getmass(hadron1code)
        mass2 = Decaytables.Getmass(hadron2code)
        
        mom3 = Mass/2
        mom3 *= self.Sqrtlambda1(Mass,mass1,mass2)
        randomphi = random.random()
        randomtheta = random.uniform(-1,1)
        phi = 2.*numpy.pi*randomphi
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        theta = numpy.pi*randomtheta
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        mom = Fourvector.vector4D(0.,sintheta*cosphi,sintheta*sinphi,costheta)
        mom *= mom3
        E1 = numpy.sqrt(mass1*mass1+mom3*mom3)
        E2 = numpy.sqrt(mass2*mass2+mom3*mom3)
        finalvec1 = Fourvector.vector4D()
        finalvec1[0] = E1
        finalvec1 += mom
        finalvec2 = Fourvector.vector4D()
        finalvec2[0] = E2
        finalvec2 -= mom
        finalvec1 = boost/finalvec1
        finalvec2 = boost/finalvec2
        
        ''' reset the final vecs to represent massless diquarks. Just need the flavour assignment '''
        finalvec1 = copy.deepcopy(string.Getqvector4D())
        finalvec2 = copy.deepcopy(string.Getqbarvector4D())
        
        finalhadron1 = Particle.particle4D(hadron1code,finalvec1)
        finalhadron2 = Particle.particle4D(hadron2code,finalvec2)
        
        if hadron1code > 0:
            newstring = MRS.mrs4D(finalhadron1,finalhadron2)
            
        if hadron1code < 0:
            newstring = MRS.mrs4D(finalhadron2,finalhadron1) # try to keep the negative / positive ids in the relative placements (kind of important)
            
        
        #print newstring.Getstringcontent(), "CHECKING THE TWO BODY MODULE"
        return newstring 
      
        
    def Sqrtlambda1(self,M,m1,m2):     #takes the unsquard masses as arguments
        M_2 = M*M
        m1_2 = m1**2
        m2_2 = m2**2
        lamb = (M_2-m1_2-m2_2)*(M_2-m1_2-m2_2) - 4*m1_2*m2_2
        sqrtlamb = numpy.sqrt(lamb)
        sqrtlamb = sqrtlamb / M_2
        return sqrtlamb
        
                
        
    def Sqrtlambda2(self,M_2,m1_2,m2_2): # takes the masses squared as arguments
        lamb = (M_2-m1_2-m2_2)*(M_2-m1_2-m2_2) - 4*m1_2*m2_2
        sqrtlamb = numpy.sqrt(lamb)
        sqrtlamb = sqrtlamb /M_2
        return sqrtlamb
        
        
        
    def Failsafe1(self,string):
        ''' Get Produced hadron & photon masses '''
        Decaytables = Decaydictionary.Decays()
        Failsafe1hadron = Decaytables.GetFailsafe1hadron(string)
        if Failsafe1hadron == 0:
            Failsafe1hadron = 22 #  set to a photon
            print "Failsafe1 being fed an unphysical string - setting the hadron to photon to avoid crash but investigate the problem"
            print string.Getstringmass(), "stringmass"
            print string.Getstringcontent(), "stringcontent"
        #print Failsafe1hadron
        hadronmass = Decaytables.Getmass(Failsafe1hadron) # mass1
        #print hadronmass, "hadronmass"
                
        photonmass = 0 # mass 2
        ''' Get boost class which contains useful boosting functions '''
        string4vector = string.Gettotalvector4D()       
        boost = Boosts.boost4D(string4vector)
        
        ''' CMF mass of the string (mass of decaying object) '''
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0]
        mom3 = Mass/2
        mom3 *= self.Sqrtlambda1(Mass,hadronmass,photonmass)
        
        ''' Do the random angle generation '''
        randomphi = random.random()
        randomtheta = random.uniform(-1,1)
        phi = 2.*numpy.pi*randomphi
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        theta = numpy.pi*randomtheta
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        
        ''' start setting the 4 vectors of the products '''
        
        mom = Fourvector.vector4D(0.,sintheta*cosphi,sintheta*sinphi,costheta)
        mom *= mom3
        E1 = numpy.sqrt(hadronmass*hadronmass+mom3*mom3)
        E2 = numpy.sqrt(photonmass*photonmass+mom3*mom3)
        finalvec1 = Fourvector.vector4D()
        finalvec1[0] = E1
        finalvec1 += mom
        finalvec2 = Fourvector.vector4D()
        finalvec2[0] = E2
        finalvec2 -= mom
        
        ''' Boost back into the lab frame '''
        
        finalvec1 = boost/finalvec1 #hadron
        finalvec2 = boost/finalvec2 #photon
        
        ''' Create the particles and return them in a list '''
        Failsafe1hadron = self.Checkforkaons(Failsafe1hadron)
        
        finalhadron = Particle.particle4D(Failsafe1hadron,finalvec1)
        finalphoton = Particle.particle4D(22,finalvec2)
        
        if abs(finalhadron.Getcode()) == 311:
            print "FOUND IT"
            sys.exit()
        
        return [finalhadron, finalphoton]
        
        

    def threebodyparticledecay(self,decayingparticle,producthadrons): # product hadrons = list of particles (from chosendecay in hadrondecay module)
        
        decayingvector = decayingparticle.Getvector4D()
        hadron1code = producthadrons[0]
        hadron2code = producthadrons[1]
        hadron3code = producthadrons[2]
        
        Decaytables = Decaydictionary.Decays()
        m1 = Decaytables.Getmass(hadron1code)
        m2 = Decaytables.Getmass(hadron2code)
        m3 = Decaytables.Getmass(hadron3code)
                
        
        boost = Boosts.boost4D(decayingvector)
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0]
        
        ''' PLACE HOLDER HERE TO TRY AND GET A RESULT. NOT TOO CONCERNED ABOUT MOMENTA YET '''
        finalvec1 = Fourvector.vector4D(m1,0,0,0)
        finalvec2 = Fourvector.vector4D(m2,0,0,0)
        finalvec3 = Fourvector.vector4D(m3,0,0,0)
        
        product1 = Particle.particle4D(hadron1code,finalvec1)
        product2 = Particle.particle4D(hadron2code,finalvec2)
        product3 = Particle.particle4D(hadron3code,finalvec3)
        
        return [product1,product2,product3]
       
        while True:
            print "ARE WE STUCK IN THE TIME WARP LOOP. 3 - BODY DECAY"
            
            
            ''' Constrain m12 '''
            
            m12_min = m1 + m2
            m12_max = Mass - m3
            
            D_m12 = m12_max - m12_min
            rng = random.random()
            
            m12 = m12_min + rng * D_m12 # select a number randomly (uniform) between m12min and m12max.
            
            ''' Use m12 to constrain m23 '''
            
            E2 = (m12*m12 - m1*m1 + m2*m2) / (2*m12) # the energy of particle 2 in the rest frame of m12
            E3 = (Mass*Mass - m12*m12 - m3*m3) / (2*m12) # same but for 3
            
            m23_min = (E2 + E3)*(E2 + E3) - (numpy.sqrt(E2*E2 - m2*m2) + numpy.sqrt(E3*E3 - m3*m3))**2 # these two are m23 squared
            m23_max = (E2 + E3)*(E2 + E3) - (numpy.sqrt(E2*E2 - m2*m2) - numpy.sqrt(E3*E3 - m3*m3))**2
            
            D_m23 = numpy.sqrt(m23_max) - numpy.sqrt(m23_min)
            
            rng2 = random.random()
            
            m23 = numpy.sqrt(m23_min) + rng2 * D_m23
            
            ''' Calculate m13 from m12, m23 '''
            
            m13_2 = Mass*Mass + m1*m1 + m2+m2 + m3*m3 - m12*m12* - m23*m23
            m13 = numpy.sqrt(m13_2)
            
            ''' Hardcode the momenta functions ''' # Can use lambda as well to get the expressions
            
            p1 = (Mass/2) * self.Sqrtlambda1(Mass,m23,m1)
            #print p1, "p1"
            p2 = (Mass/2) * self.Sqrtlambda1(Mass,m13,m2)
            #print p2, "p2"
            p3 = (Mass/2) * self.Sqrtlambda1(Mass,m12,m3)
            #print p3, "p3"        
            
            ''' Generate the angles phi is isotropic '''
            
            randomphi = random.random()
            phi = 2.*numpy.pi*randomphi
            sinphi = numpy.sin(phi)
            cosphi = numpy.cos(phi)
            
            ''' Generate the thetas in the plane of decay '''
            
            # First momentum theta is random.
            
            randomtheta1 = random.uniform(-1,1)
            theta1 = numpy.pi*randomtheta1
            sintheta1 = numpy.sin(theta1)
            costheta1 = numpy.cos(theta1)
            
            # Second is also random, and helps determine the second. However have to generate sensible / possible angles
            for i in range(100):
                randomtheta2 = random.uniform(-1,1)
                theta2 = numpy.pi*randomtheta2
                sintheta2 = numpy.sin(theta2)
                costheta2 = numpy.cos(theta2)
                
                # Third constrained by the first two in order to balance momentum.
                
                sintheta3 = (-p1*sintheta1 - p2*sintheta2) / p3
                costheta3 = (-p1*costheta1 - p2*costheta2) / p3
                
                
                
                if abs(sintheta3) <= 1 and abs(costheta3) <= 1:
                    
                    ''' Build the 4 vectors '''
                    
                    vecE1 = numpy.sqrt(m1*m1 + p1*p1)
                    #print vecE1, "vecE1"
                    vecE2 = numpy.sqrt(m2*m2 + p2*p2)
                    #print vecE2, "vecE2"
                    vecE3 = numpy.sqrt(m3*m3 + p3*p3)
                    #print vecE3, "vecE3"
                    
                    finalvec1 = Fourvector.vector4D(vecE1,p1*sintheta1*cosphi,p1*sintheta1*sinphi,costheta1)
                    finalvec1 = boost/finalvec1
                    
                    finalvec2 = Fourvector.vector4D(vecE2,p2*sintheta2*cosphi,p2*sintheta2*sinphi,costheta2)
                    finalvec2 = boost/finalvec2
                    
                    finalvec3 = Fourvector.vector4D(vecE3,p3*sintheta3*cosphi,p3*sintheta3*sinphi,costheta3)
                    finalvec3 = boost/finalvec3
                    
                    
                    ''' Construct the particles '''
                    product1 = Particle.particle4D(hadron1code,finalvec1)
                    product2 = Particle.particle4D(hadron2code,finalvec2)
                    product3 = Particle.particle4D(hadron3code,finalvec3)
                    
                    return [product1,product2,product3]
                    
                    
                    
    def twobodyparticledecay(self,particle,products): # takes a list of two particles and decays them.
        
        ''' Get Produced hadron masses '''
        Decaytables = Decaydictionary.Decays()
        particle4vector = particle.Getvector4D()
        product1code = products[0]
        product2code = products[1]
        
        mass1 = Decaytables.Getmass(product1code)
        mass2 = Decaytables.Getmass(product2code)
        
       
        ''' Get boost class which contains useful boosting functions '''
        
        #print "boost in two body called"
        boost = Boosts.boost4D(particle4vector)
        
        ''' CMF mass of the string (mass of decaying object) '''
        cmsvector = boost.GetCMF()
        Mass = cmsvector[0]
        mom3 = Mass/2
        mom3 *= self.Sqrtlambda1(Mass,mass1,mass2)
        
        ''' Do the random angle generation '''
        randomphi = random.random()
        randomtheta = random.uniform(-1,1)
        phi = 2.*numpy.pi*randomphi
        sinphi = numpy.sin(phi)
        cosphi = numpy.cos(phi)
        theta = numpy.pi*randomtheta
        sintheta = numpy.sin(theta)
        costheta = numpy.cos(theta)
        
        ''' start setting the 4 vectors of the products '''
        
        mom = Fourvector.vector4D(0.,sintheta*cosphi,sintheta*sinphi,costheta)
        mom *= mom3
        E1 = numpy.sqrt(mass1*mass1+mom3*mom3)
        E2 = numpy.sqrt(mass2*mass2+mom3*mom3)
        finalvec1 = Fourvector.vector4D()
        finalvec1[0] = E1
        finalvec1 += mom
        finalvec2 = Fourvector.vector4D()
        finalvec2[0] = E2
        finalvec2 -= mom
        
        ''' Boost back into the lab frame '''
        print  "Is the crash here"
        finalvec1 = boost/finalvec1 #hadron
        finalvec2 = boost/finalvec2 #photon
        
        ''' Create the particles and return them in a list '''
        
        finalproduct1 = Particle.particle4D(product1code,finalvec1)
        finalproduct2 = Particle.particle4D(product2code,finalvec2)
        
        return [finalproduct1, finalproduct2]
        
                    
    def Checkforkaons(self,code):
        if abs(code) == 311:
            rng = random.getrandbits(1) # generate a number which is 0 or 1. Line here in case of method change to speed up code.
            if rng == 0:
                return 130
            if rng == 1:
                return 310
                
        return code                    
            




            
'''
twobodyclass = twobodydecay4D()

vector = Fourvector.vector4D(numpy.sqrt(101),2.,4.,5.)
decayingparticle = Particle.particle4D(333,vector)

products = [111,111,111]

testing = twobodyclass.threebodydecay(decayingparticle,products)

print testing
        
'''        
    