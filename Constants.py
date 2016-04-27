''' File for storing and adjusting constants used in generation of functions '''

import numpy

''' Defines the numerical limit, to which numerical errors are tolerated '''

def Numerical_limit():
    numerical_accuracy = 1.e-12
    return numerical_accuracy
    
def detector_limit():
    limit = 1.e-1
    return limit
    
def baryonproductionthreshold():
    bpt = 3 # threshold set to which baryon production is vetoed in strings
    return bpt

''' These values are relevant for the stringdecay module - pertain to parameters or constants called in that code '''

def aLund():
    a = 2   # kinematic parameter of the lund string decay
    return a
    
def bLund():
    b = 2 # kinematic parameter of the lund string decay
    return b
        
def kappa():
    kappa = 0.2 # alternate expression of sigma. as yet unused
    return kappa
    
def kperpsigma():
    sigma = 0.300 # Sigma parameter in the code. Representative of the string colour force, and affects the average k_perp of produced hadrons as well as their multiplicities
    return sigma

def gcc(): # gluon combination constnt. referred to as delta in the thesis body, controls the gluon recombination threshold
    gcc = 4
    return gcc

def gmm(): # gluon min mass modifier. Set to 2 by default, but can change the minimum required energy of gluons to not be absorbed within the simulation
    gmm = 2
    return gmm
    
def gsigma():
    gsigma = 0.7030 # variance" of the g function for the kperp hit miss - there is a code to generate this value
    return gsigma

def g_kperpA():     #"amplitude" of the g function for kperp hit miss
    g_kperpA = 1.1
    return g_kperpA
    
def hbar_gevs():
    hbar = 6.582119514e-25
    return hbar
    
def speedoflight():
    c = 2.9979e8
    return c
    
    
''' Constants pertaining to the available decays code - determine things like the relative frequencies of the strange mesons and the eta eta' mesons '''

def decaythreshold():
    threshold = 4 # number of stringmasses required for a decay. as yet unused
    return threshold
    
def minmassmodifier(): # the minimum mass modifier when calculating the minimum pair mass of a quark, diquark string. Used as an approximater to calculate the FS2 masses of hadrons
    mmf = 0.85 # 
    return mmf
    
def doublebaryonmassmodifier():
    dbmm = 0.85 # similar in function to the minmass modifier, just for two diquark strings
    return dbmm

def strange_suppression_factor():
    SSF = 1.1 # constant to vary the relative production of strange mesons
    return SSF

def charm_scaling_factor():
    csf = 1 # as yet unused
    return csf
    
def bottom_scaling_factor():
    bsf = 1 # as yet unused
    return bsf

def eta_scaling_factor():
    ESF = 2 # change relative production of eta mesons
    return ESF
    
def eta_prime_scaling_factor():
    EPSF = 10000  # change relative production of eta prime mesons
    return EPSF
    
def vector_scaling_factor():
    vsf = 50000 # alter relative produciton of vector mesons
    return vsf
    

def eta_mixing_angle(): # "wavefunction theta"
    theta = -11.5*2*numpy.pi/360  # first value is numpy in degrees - then need to convert
    return theta
    
def decay_weight_sigma():
    sigma = kperpsigma()
    return sigma     #same value as the kperp distribution

def pi0_wf2(): # wf is essentially symmetric (equally weighted when squared up
    wf = 0.5
    return [0.5,0.5]
    
def eta_wf2(): # returns wf2 for d,u,s repspectivley
    theta = eta_mixing_angle()
    ''' uubar|eta = 1/rt(6) costheta - 1/rt(3) sintheta '''
    uquark = (1./numpy.sqrt(6)) * numpy.cos(theta) - (1./numpy.sqrt(3)) * numpy.sin(theta)
    uquark_2 = uquark*uquark
    
    ''' ddbar|eta = uubar|eta '''
    dquark = (1./numpy.sqrt(6)) * numpy.cos(theta) - (1./numpy.sqrt(3)) * numpy.sin(theta)
    dquark_2 = dquark*dquark
    
    ''' ssbar|eta = -2/rt(6) costheta - 1/rt(3) sintheta '''
    squark = (-2./numpy.sqrt(6)) * numpy.cos(theta) - (1/numpy.sqrt(3)) * numpy.sin(theta)
    squark_2 = abs(squark) * abs(squark)
    
    return [dquark_2,uquark_2,squark_2]
    
def eta_prime_wf2():
    theta = eta_mixing_angle()
    ''' uubar|etaprime = 1/rt(3) costheta + 1/rt(6) sintheta '''
    uquark = 1./numpy.sqrt(3) * numpy.cos(theta) + 1./numpy.sqrt(6) * numpy.sin(theta)
    uquark_2 = uquark*uquark
    
    ''' uubar|etaprime = ddbar|etaprime '''
    dquark = 1./numpy.sqrt(3) * numpy.cos(theta) + 1./numpy.sqrt(6) * numpy.sin(theta)
    dquark_2 = dquark*dquark
    
    ''' ssbar|etaprime = 1/rt(3) costheta - 2./rt(6) sintheta '''
    squark = 1./numpy.sqrt(3) * numpy.cos(theta) - 2./numpy.sqrt(6) * numpy.sin(theta)
    squark_2 = abs(squark) * abs(squark)
    
    return [dquark_2, uquark_2, squark_2]
    
    
    
    
''' Baryon wavefunction overlap with diquarks '''


    
    
    
''' Baryon scaling factors '''

def baryon_scaling_factor():
    bsf = 1 # as yet unused
    return bsf    
    
def vector_baryon_scaling_factor():
    vbsf = 1 # change production rate of vector baryons
    return vbsf
    

def strange_baryon_scaling_factor():
    sbsf = 1 # change production rate of strange baryons
    return sbsf
    
def baryon_exp_factor(): # Exponential factor used when getting the weights. Set to 1 by default.
    bef = 1
    return bef
    
    

''' Diquark scaling factors '''
    
def diquark_scaling_factor():
    dsf = 2# relative production of diquark in string decay and gluon splitting
    return dsf
    
def vector_diquark_scaling_factor():
    vdsf = 2 # relative production of vector diquark in string decay and gluon splitting
    return vdsf
    
        
def strange_diquark_scaling_factor():
    sdsf = 1 # relative production of strange diquark in string decay and gluon splitting
    return sdsf
    
def diquark_exp_factor():
    dief = 1 # Exponential factor used when getting the weights. Set to 1 by default.
    return dief    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    

    
    


