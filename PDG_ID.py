''' Class for storing pdg_ID values, and catalogueing particles - written by Krauss - augmented by Ryan '''

from __future__ import division
import numpy
import string


class PDG_Entry(object):
    """Particle properties: [pdg code,mass,width]"""    
    def __init__(self,name,mass=0.,width=0.,charge=0.):
        self.__name  = name
        self.__mass  = mass
        self.__width = width
        self.__charge = charge

    def getName(self):
        return self.__name

    def getMass(self):
        return self.__mass

    def getWidth(self):
        return self.__width
        
    def getCharge(self):
        return self.__charge

class PDG_ID(object):
    """A class containing particle type definitions"""
    ''' Masses are expressed in units of GeV: Width in GeV, with c=1 mu=2 '''
    ''' Particles with "infinite" lifetimes are given the value 1e-40 for decay width '''
    ''' Values taken from PDG online data sheets, as of 5/11/15 '''
    """Dictionary of particles and their names and properties"""    #other properties? charge, isospin etc??
    def __init__(self):
        self.__pdgs = {
            1:   PDG_Entry('d-quark',0.0,0.0,-1./3), # 0.0048
            2:   PDG_Entry('u-quark',0.0,0.0,2./3),  # 0.0023
            3:   PDG_Entry('s-quark',0.0,0.0,-1./3),  # 0.095
            4:   PDG_Entry('c-quark',1.275,0.0,2./3),
            5:   PDG_Entry('b-quark',4.66,0.0,-1./3),     # using the 1S decay mass here
            6:   PDG_Entry('t-quark',173.21,1.41,2./3),  # Width taken from PDG data (errors???)
            
            11:  PDG_Entry('e^-',0.000511,0.0,-1.0),
            12:  PDG_Entry(r'\nu_e',0.0,0.0,0.0),
            13:  PDG_Entry('\mu^-',0.106,0.0,-1.0),
            14:  PDG_Entry(r'\nu_{\mu}',0.0,0.0,0.0),
            15:  PDG_Entry('tau',1.777,2.36e-12,-1.0),   # Value for width here uncertain
            16:  PDG_Entry('nu_tau',0.0,0.0,0.0),
            21:  PDG_Entry('gluon',0.0,0.0,0.0),
            22:  PDG_Entry('\gamma',0.0,0.0,0.0),
            23:  PDG_Entry('Z boson',91.188,2.49,0.0),
            24:  PDG_Entry('W+ boson',80.385,2.09,1.0),
            25:  PDG_Entry('Higgs',125.1,0.00421,0.0),  # Value for width here uncertain
            #Mesons begin
           111:  PDG_Entry('\pi^0',0.135,7.7255e-9,0.0),  # value for width uncertain (999) means value to be determined.'
           113:  PDG_Entry(r'\rho^0',0.7753,0.1478,0.0),  #Width? Need to refresh memory of what to use (decay dependant)
           130:  PDG_Entry('K^0_L',0.497611,1.28657e-17,0.0),
           211:  PDG_Entry('\pi^+',0.1396,2.5284e-17,1.0), # Width
           213:  PDG_Entry(r'\rho^+',0.7751,0.1491,1.0), # AJUST
           221:  PDG_Entry('\eta',0.5479,131.e-8,0.0),     # Width quoted on PDG database?
           223:  PDG_Entry('\omega',0.7827,0.00849,0.0),   # average width quoted
           310:  PDG_Entry('K^0_S',0.497611,7.351e-15,0.0),
           311:  PDG_Entry('K^0$',0.497611,999,0.0),
           313:  PDG_Entry('K^{*0}',0.8958,0.0474,0.0),
           321:  PDG_Entry('K^+',0.493677,999,1.0),     # width to be confirmed
           323:  PDG_Entry('K^{*+}',0.8917,0.0508,0.0), # used HADROPRODUCED READINGS
           331:  PDG_Entry('\eta\prime',0.9578,0.000230,0.0),  # width average value, not their fit.
           333:  PDG_Entry('\phi',1.019,0.004266,0.0),  #av width quoted
           411:  PDG_Entry('D+-meson',1.8696,6.329e-10,1.0),
           413:  PDG_Entry('D*+-meson',2.01027,0.0834,1.0),
           421:  PDG_Entry('D0-meson',1.8648,1.605e-09,0.0),
           423:  PDG_Entry('D0*-meson',2.007,2.1,0.0),
           431:  PDG_Entry('D+_s-meson',1.969,1.316e-09,1.0),
           433:  PDG_Entry('D*+_s-meson',2.1121,1.9,1.0),
           441:  PDG_Entry('eta_c-meson(1S)',2.9836,32.0,0.0), 
           443:  PDG_Entry('J/psi-meson',3.0969,0.0929,0.0),
           511:  PDG_Entry('B0-meson',5.2796,4.330e-10,0.0),
           521:  PDG_Entry('B+-meson',5.2793,4.018e-10,1.0),
           531:  PDG_Entry('B0_s-meson',5.3668,4.359e-10,0.0),
           533:  PDG_Entry('B*0_s-meson',5.4158,999,0.0),    # this value of width unknown.
           541:  PDG_Entry('B+_c-meson',6.2751,1.298e-09,1.0),
           551:  PDG_Entry('eta_b-meson(1S)',9.398,10.8,0.0),
           553:  PDG_Entry('Epsilon(1S)-meson',9.460,0.054,0.0),
           # Baryons
          2212:  PDG_Entry('proton',0.9383,1e-40,1.0),        # width of p, n?
          2112:  PDG_Entry('neutron',0.9396,1e-40,0.0),     #
          
          2224:  PDG_Entry('Delta++',1.232,0.117,2.0),
          2214:  PDG_Entry('Delta+',1.232,0.117,1.0),
          2114:  PDG_Entry('Delta0',1.232,0.117,0.0),
          1114:  PDG_Entry('Delta-',1.232,0.117,-1.0),
          
          3122:  PDG_Entry('Lambda',1.115683,2.5008e-15,0.0),
          
          3222:  PDG_Entry('Sigma+',1.1894,8.20917e-15,1.0),
          3212:  PDG_Entry('Sigma0',1.1926,8.8948e-6,0.0),
          3112:  PDG_Entry('Sigma-',1.1974,4.45e-15,-1.0),
          3224:  PDG_Entry('Sigma*+',1.38280,0.036,1.0),
          3214:  PDG_Entry('Sigma*0',1.3837,0.036,0.0),
          3114:  PDG_Entry('Sigma*-',1.3872,0.0394,-1.0),
          
          3322:  PDG_Entry('Xi0',1.31486,2.2696e-15,0.0),
          3312:  PDG_Entry('Xi-',1.32171,4.0159e-15,-1.0),
          3324:  PDG_Entry('Xi*0',1.5318,0.0091,0.0),
          3314:  PDG_Entry('Xi*-',1.535,0.0099,-1.0),
          
          3334:  PDG_Entry('Omega-',1.67245,8.017e-15,-1.0),
        # Diquarks coded in. assume that u,d 300Mev, s 450Mev
        # Width is irrelevant. they should never be in a situation where they would have to decay into product states.
          1103:  PDG_Entry('dd(1)',0.60,1e-40,-2./3),
          2101:  PDG_Entry('ud(0)',0.60,1e-40,1./3),
          2103:  PDG_Entry('ud(1)',0.60,1e-40,1./3),
          2203:  PDG_Entry('uu(1)',0.60,1e-40,4./3),
          3101:  PDG_Entry('sd(0)',0.75,1e-40,-2./3),
          3103:  PDG_Entry('sd(1)',0.75,1e-40,-2./3),
          3201:  PDG_Entry('su(0)',0.75,1e-40,1./3),
          3203:  PDG_Entry('su(1)',0.75,1e-40,1./3),
          3303:  PDG_Entry('ss(1)',0.90,1e-40,-2./3),
          
          

  
        }

        self.__quarks   = [1,2,3,4,5,6]
        self.__leptons  = [11,12,13,14,15,16]
        self.__mesons   = [111,113,130,211,213,221,223,310,311,313,321,323,331,333,411,413,421,423,431,433,441,443,511,521,531,533,541,551,553]
        self.__baryons  = [2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334]
        self.__hadrons  = self.__mesons + self.__baryons # HADRONS IS MESONS + BARYONS
        self.__diquarks = [1103,2101,2103,2203,3101,3103,3201,3203,3303]
        self.__gluons = [21]
    """Simple access methods"""

    def name(self,code):
        if self.__pdgs.has_key(code):
            return self.__pdgs[code].getName()
        return "not known!"

    def mass(self,code):
        if self.__pdgs.has_key(code):
            return self.__pdgs[code].getMass()
        return -1.

    def width(self,code):
        if self.__pdgs.has_key(code):
            return self.__pdgs[code].getWidth()
        return -1.
        
    def charge(self,code):
        if self.__pdgs.has_key(code):
            return self.__pdgs[code].getCharge()
        return code, "is not known"

    def code(self,name):
        for i in self.__pdgs:
            if self.__pdgs[i].getName()==name:
                return i
        return name+" is not known"

    def has_key(self,code):
        if self.__pdgs.has_key(code):
            return True
        return False

    def isLepton(self,code):
        return self.__leptons.count(code)>0

    def isQuark(self,code):
        return self.__quarks.count(code)>0
        
    def isMeson(self,code):
        return self.__mesons.count(code)>0

    def isHadron(self,code):
        return self.__hadrons.count(code)>0
        
    def isBaryon(self,code):
        return self.__baryons.count(code)>0
        
    def isDiquark(self,code):
        return self.__diquarks.count(code)>0

    def isFermion(self,code):
        return self.isLepton(code) or self.isQuark(code)

    def isGluon(self,code):
        return self.__gluons.count(code)>0
        
    def Getmesonlist(self):
        return self.__mesons
        
    def Getbaryonlist(self):
        return self.__baryons
        
    def Gethadronlist(self):
        return self.__hadrons
        
    def Getdiquarklist(self):
        return self.__diquarks

PDG_IDs = PDG_ID()

class PDG_Type(object):
    """The actual type of particle - allows for anti-particles"""
    def __init__(self,pdgcode):
        assert PDG_IDs.has_key(abs(pdgcode))
        self.__anti = False
        if pdgcode<0:
            self.__anti = True
        self.__pdgcode = pdgcode
        
        
    def isAnti(self):
        return self.__anti

    def isFermion(self):
        return PDG_IDs.isFermion(abs(self.__pdgcode))

    def isQuark(self):
        return PDG_IDs.isQuark(abs(self.__pdgcode))
        
    def isMeson(self):
        return PDG_IDs.isMeson(abs(self.__pdgcode))

    def isLepton(self):
        return PDG_IDs.isLepton(abs(self.__pdgcode))

    def isHadron(self):
        return PDG_IDs.isHadron(abs(self.__pdgcode))
        
    def isBaryon(self):
        return PDG_IDs.isBaryon(abs(self.__pdgcode))
        
    def isDiquark(self):
        return PDG_IDs.isDiquark(abs(self.__pdgcode))

    def isGluon(self):
        return PDG_IDs.isGluon(abs(self.__pdgcode))

    def code(self):
        return self.__pdgcode

    def name(self):
        name = ""
        tempname = PDG_IDs.name(abs(self.__pdgcode))
        if self.__anti:
            if str(tempname) == r"\nu_{\mu}":
                return r"$\bar{\nu}_{\mu}$"
            if str(tempname) == r"\nu_e":
                return r"$\bar{\nu}_e$"
            
            if "+" in str(tempname): # charged antiparticles
                tempname1 = string.replace(tempname,"+", "-", 1)
            
                name += "$"
                name += tempname1
                name += "$"
                return name
            
            if "-" in str(tempname): # charged antiparticles
                tempname2 = string.replace(tempname,"-", "+", 1)
            
                name += "$"
                name += tempname2
                name += "$"
                return name
            name += r"$\bar{"
            name += tempname[0]
            name += "}"
            name += tempname[1:]
            name += "$"
            return name
        name += "$"    
        name.replace(" ", "")
        name += PDG_IDs.name(abs(self.__pdgcode))
        name += "$"
        
        return name

    def mass(self):
        return PDG_IDs.mass(abs(self.__pdgcode))

    def width(self):
        return PDG_IDs.width(abs(self.__pdgcode))
        
    def charge(self):
        if self.__anti:
            return PDG_IDs.charge(abs(self.__pdgcode)) * -1.
        return PDG_IDs.charge(abs(self.__pdgcode))
            
        
        
