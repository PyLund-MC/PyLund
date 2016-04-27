''' Hadrondecaymodule - contains the tables and information pertaining to the further decay of hadrons after string decay has finished '''

import numpy
#import Boosts
import Fourvector
import Particle
import PDG_ID
import random
import MRS
#import sys
import Constants
#import Stringdecay
import copy
import Twobodydecay
import Decaydictionary
import Radioactivedecay


class Hadrondecayinput(object):
    def __init__(self, *args):
        
        self.__products = args
        
    def Getdecays(self):
        return self.__products
        
    def Printdecays(self):
        print self.__products


class hadrondecay4D(object):
    
    def __init__(self):
        
        self.__Hadrondecays = {
           111:   Hadrondecayinput(([22,22],0.99), ([11,-11,22],0.01)),   #[pi0] 2photon - 99% : e+e- photon - 1%
           113:   Hadrondecayinput(([111,111],0.50), ([211,-211],0.50)), # [rho0] 2pi0 - 50% : pi+pi- - 50%
           130:   Hadrondecayinput(([211,11,-12],0.20275), ([-211,-11,12],0.20275), ([211,13,-14],0.1352), ([-211,-23,14],0.1352), ([111,111,111],0.1952), ([211,-211,111],0.1254)),  #[K0L]   pi+ e- anti(nu_e) - 20.275% : pi- e+ nu_e - 20.275% : pi+ mu- anti(nu_mu) - 13.52% : pi- mu+ nu_mu - 13.52% : 3pi0 - 19.52% : pi+pi-pi0 - 12.54%
           
           211:   Hadrondecayinput(([-13,14],1)),        #[pi+]   mu+ nu_mu - 100%
          -211:   Hadrondecayinput(([13,-14],1)),        #[pi-]   mu- anti(nu_mu) - 100%
           213:   Hadrondecayinput(([211,111],1)),
          -213:   Hadrondecayinput(([-211,111],1)),
           221:   Hadrondecayinput(([22,22],0.3941), ([111,111,111],0.327), ([211,-211,111],0.229), ([211,-211,22],0.042)),    #[eta]   2photon - 39.41% : 3pi0 - 32.7% : pi+pi-pi0 - 22.9% : pi+pi- photon - 4.2%
           223:   Hadrondecayinput(([211,-211,111],0.892), ([111,22],0.0828), ([211,-211],0.0153)),      #[omega]   pi+pi-pi0 - 89.2% : pi0 photon - 8.28% : pi-pi+ - 1.53%
           
           #311 not coded as its the K0L and K0S mesons. You never see a 311 meson per say- only variations of short and long kaons
           310:   Hadrondecayinput(([111,111],0.3069), ([211,-211],0.6920)),     # [K0S]    2pi0 - 30.69% : pi-pi+ - 69.2%
           313:   Hadrondecayinput(([130,111], 0.16667), ([310,111], 0.16667), ([211,-321], 0.33333), ([-211,321], 0.33333)), # [K0*]  K0L pi0 - 0.167% : K0S, pi0 - 0.167% : pi+K- - 0.33% : pi-K+ - 0.33%
          -313:   Hadrondecayinput(([130,111], 0.16667), ([310,111], 0.16667), ([211,-321], 0.33333), ([-211,321], 0.33333)), # [K0*bar]  K0L pi0 - 0.167% : K0S, pi0 - 0.167% : pi+K- - 0.33% : pi-K+ - 0.33% SAME AS 313
           321:   Hadrondecayinput(([-13,14],0.6356), ([111,-11,12],0.0507), ([111,-13,14],0.03352), ([211,111],0.2067), ([211,111,111],0.0170), ([211,211,-211],0.0583)),    #[K+] mu+ nu_mu - 63.56% : pi0 e+ nu_e - 5.07% : pi0 mu+ nu_mu - 3.352% : pi+pi0 - 20.67% : 2pi+pi- - 5.583%
          -321:   Hadrondecayinput(([13,-14],0.6356), ([111,11,-12],0.0507), ([111,13,-14],0.03352), ([-211,111],0.2067), ([-211,111,111],0.0170), ([211,-211,-211],0.0583)),    #[K-] mu- anti(nu_mu) - 63.56% : pi0 e- anti(nu_e) - 5.07% : pi0 mu- anti(nu_mu) - 3.352% : pi-pi0 - 20.67% : 2pi-pi+ - 5.583%
           323:   Hadrondecayinput(([321,111],0.50), ([130,211], 0.25), ([310,211], 0.25)), # [K+*] K+pi0 - 50% : K0Lpi+ - 25% : K0Spi+ - 25%
          -323:   Hadrondecayinput(([-321,111],0.50), ([130,-211], 0.25), ([310,-211], 0.25)), # [K-*] K-pi0 - 50% : K0Lpi- - 25% : K0Spi- - 25%
          
           331:   Hadrondecayinput(([211,-211,221],0.429), ([211,-211,22],0.291), ([111,111,221],0.222), ([223,22],0.0275), ([22,22],0.022)),    #[etaprime]   pi+pi-eta - 42.9% : pi+pi-photon - 29.1% : 2pi0 eta - 22.2% : omega photon - 2.75% : 2photon - 2.2%
           # INcluding non-resonant decay argument here. Just put the final state and skipped the rho state
           333:   Hadrondecayinput(([321,-321],0.489), ([310,130],0.342), ([211,-211,111],0.1532), ([221,22],0.01309)), #[phi]   K+K- - 48.9% : K0L,K0S - 34.2% : pi0pi+pi- 15.32% NOTE, ACTUAL DECAY LISTED AS (RHOPI + PI0PI-PI+) : eta photon - 1.309%
           
           # 333 3RD DECAY CAN ALSO BE SEEN AS RHOPI + PI-PI+PI0. MAYBE ASK KRAUSS AGAIN BRIEFLY IF THIS TREATMENT IS FINE.
          2212:   Hadrondecayinput(([21,21],1)), # [proton] gg - 100% (Dummy variable. This particle should NEVER deacy 
         -2212:   Hadrondecayinput(([21,21],1)), # [antiproton] gg -100% (Again this should never occur) 
          2112:   Hadrondecayinput(([2212,11,-12],1)), # [n] p,e-,antinu_e - 100%
         -2112:   Hadrondecayinput(([-2212,-11,12],1)),  #[anti-neutron] antip,e+,nu_e - 100% (Neutrons should also nver decay) 
          
          2224:   Hadrondecayinput(([2212,211],1)), # [Delta++] p,pi+ - 100%
         -2224:   Hadrondecayinput(([-2212,-211],1)), # [anti Delta++] antip,pi- - 100%
          2214:   Hadrondecayinput(([2212,111],0.5), ([2112,211],0.5)), # [Delta+] p,pi0 - 50% : n,pi+ - 50%
         -2214:   Hadrondecayinput(([-2212,111],0.5), ([-2112,-211],0.5)), # [antidelta+] antip,pi0 - 50% : antin,pi- - 50% 
          2114:   Hadrondecayinput(([2112,111],1)), # [Delta0] n,pi0 - 100%
         -2114:   Hadrondecayinput(([-2112,111],1)), # [antidelta0] antin,pi0 - 100% 
          1114:   Hadrondecayinput(([2112,-211],1)), # [Delta-] n,pi- - 100%
         -1114:   Hadrondecayinput(([-2112,211],1)), # [antidelta-] antin,pi+  -100% 
          
          3122:   Hadrondecayinput(([2212,-211],0.639), ([2112,111],0.358)), # [lambda] p,pi- - 63.9% : n,pi0 - 35.8%
         -3122:   Hadrondecayinput(([-2212,211],0.639), ([-2112,111],0.358)), # [antilambda] antip,pi+ - 63.9% : antin,pi0 - 35.8% 
          
          3222:   Hadrondecayinput(([2212,111],0.5157), ([2112,211],0.4831)), # [Sigma+] p,pi0 - 51.57% : n,pi+ - 48.31%
         -3222:   Hadrondecayinput(([-2212,111],0.5157), ([-2112,-211],0.4831)), # [antiSigma+] antip,pi0 - 51.57% : antin,pi- - 48.31% 
          3212:   Hadrondecayinput(([3122,22],1)), # [Sigma0] Lambda,photon - 100%
         -3212:   Hadrondecayinput(([-3122,22],1)), # [antiSigma0] antilambda,photon - 100% 
          3112:   Hadrondecayinput(([2112,-211],1)), # [Sigma-] n,pi- - 100%
         -3112:   Hadrondecayinput(([-2112,211],1)), # [AntiSigma-] antin,pi+ - 100% 
          3224:   Hadrondecayinput(([3122,211],0.87), ([3212,211],0.0585), ([3222,111],0.0585)), # [Sigma*+] lambda,pi+ - 87& : Sigma0,pi+ - 5.85% : Sigma+,pi0  - 5.85% : 
         -3224:   Hadrondecayinput(([-3122,-211],0.87), ([-3212,-211],0.0585), ([-3222,111],0.585)), #[AntiSigma*+] antilambda,pi- - 87% : antiSigma0,pi- - 5.585% : AntiSigma+,pi0 - 5.85% : 
          3214:   Hadrondecayinput(([3122,111],0.87), ([3112,211],0.039), ([3212,111],0.039), ([3222,-211],0.039), ([3122,22],0.013)), # [Sigma*0] lambda,pi0 - 87% : sigma-,pi+ - 3.9% : sigma0,pi0 - 3.9% : sigma+pi- - 3.9% : lambda,photon - 1.25%
         -3214:   Hadrondecayinput(([-3122,111],0.87), ([-3112,-211],0.039), ([-3212,111],0.039), ([-3222,211],0.039), ([-3122,22],0.013)), # [antiSigma*0] antilambda,pi0 - 87% : antisigma-,pi- - 3.9% : antisigma0,pi0 - 3.9% : antisigma+pi+ - 3.9% : antilambda,photon - 1.25%
          3114:   Hadrondecayinput(([3122,-211],0.87), ([3112,111],0.0585), ([3212,-211],0.0585)), # [sigma*-] lambda,pi- - 87% : sigma-pi0 - 5.85% : sigma0pi- - 5.85%
         -3114:   Hadrondecayinput(([-3122,211],0.87), ([-3112,111],0.0585), ([-3212,211],0.0585)), # [antisigma*-] antilambda,pi+ - 87% : antisigma-pi0 - 5.85% : sigma0pi+ - 5.85%
         
          3322:   Hadrondecayinput(([3122,111],1)), # [Xi0] lambda,pi0 - 100%
         -3322:   Hadrondecayinput(([-3122,111],1)), # [antiXi0] antilambda,pi0 - 100% 
          3312:   Hadrondecayinput(([3122,-211],1)), # [Xi-] lambda,pi- - 100%
         -3312:   Hadrondecayinput(([-3122,211],1)), # [AntiXi-] antilambda,pi+ - 100%
          3324:   Hadrondecayinput(([3322,111],0.5), ([3312,211],0.5)), # [Xi*0] Xi0,pi0 - 50% : Xi-pi+ - 50%
         -3324:   Hadrondecayinput(([-3322,111],0.5), ([-3312,-211],0.5)), # [antiXi*0] antiXi0,pi0 - 50% : antiXi-pi- - 50% 
          3314:   Hadrondecayinput(([3322,-211],0.5), ([3312,111],0.5)), # [Xi*-] Xi0,pi- - 50% : Xi-,pi0 - 50%
         -3314:   Hadrondecayinput(([-3322,211],0.5), ([-3312,111],0.5)), # [antiXi*-] antiXi0,pi+ - 50% : antiXi-,pi0 - 50%
         
          3334:   Hadrondecayinput(([3122,-321],0.678), ([3322,-211],0.236), ([3312,111],0.086)), # [Omega-] lambda,K- - 67.8% : Xi0,pi- - 23.6% : Xi-,pi0 - 8.6%
         -3334:   Hadrondecayinput(([-3122,321],0.678), ([-3322,211],0.236), ([-3312,111],0.086)), # [antiOmega-] antilambda,K+ - 67.8% : antiXi0,pi+ - 23.6% : antiXi-,pi0 - 8.6% 
          
          
          
          
             
           }
           
           
    def PrintHadrons(self):
        print self.__Hadrondecays
        
    def Printdecays(self,code):
        if self.__Hadrondecays.has_key(code):
            print self.__Hadrondecays[code].Getdecays()
        return "Hadron code not found in dictionary. Check input"
        
    def Getdecays(self,code):
        if self.__Hadrondecays.has_key(code):
            return self.__Hadrondecays[code].Getdecays()
        return "Hadron code not found in dictionary. Check input"
            
    def Getproducts(self,code):
        hadronlists = []
        #print code, "CODE CHECK"
        if self.__Hadrondecays.has_key(code):
            decays = self.__Hadrondecays[code].Getdecays()
            for i in range(len(decays)):
                hadronlists.append(decays[i][0])
        return hadronlists        
                
    def Getprobabilities(self,code):
        problists = []
        if self.__Hadrondecays.has_key(code):
            decays = self.__Hadrondecays[code].Getdecays()
            for i in range(len(decays)):
                problists.append(decays[i][1])
        return problists                
        
        
    def Choosehadrondecay(self,code):
        print code, "CHECK FOR CRASH"
        problist = self.Getprobabilities(code)
        productlist = self.Getproducts(code)
        sumprobs = sum(problist)
        normalised_probs = []
        for i in range(len(problist)):
            newprob = problist[i] / sumprobs
            normalised_probs.append(newprob)
        no_probs = len(problist)
        
        rng = numpy.random.choice(no_probs,1,p=normalised_probs) # chosen index in the array
        
        chosendecay = productlist[rng]
        
        return chosendecay
        
    def Decayparticlelist(self,particlelist):
        
        stableparticles = []
        templist = particlelist # use a psudeo third list to get around pythons inability to edits lists that its iterating over.
        c = Constants.speedoflight()
        unstableparticles = particlelist # guilty until proven innocent!
        twobodyclass = Twobodydecay.twobodydecay4D()
        numberdecays = 1 # if a particle decays and forms products - then the products may decay further. this counter tells whether to keep looping or not.
        
        while numberdecays != 0:
            particlelist = templist
            templist = [] # reinitialise the templist to contain the next set of prodcut
            numberdecays = 0 # set to zero. if different that means that somethihng has decayed
            for i in range(len(particlelist)):
                
                self.Checkforkaons(particlelist[i])
                
                particlecode = particlelist[i].Getcode()
                print particlecode, "PARTICLECODE"
                
                
                
                if abs(particlecode) <= 100:
                    stableparticles.append(particlelist[i])
                    continue
                
                particlelifetime = particlelist[i].Getlifetime()
                Time = Radioactivedecay.radioactivedecay4D(particlelifetime)
                
                distance = Time * c
                #print distance, "Distacne"
                
                
                if distance >= Constants.detector_limit():  #check particle stability. if the particle gets past the detector 10CM, LIMIT, then it doesnt decay. add it to the final list.
                    stableparticles.append(particlelist[i])
                    continue
                    
                ''' If the particle decays quick enough, it will create some products. Kinematics held using standard 2/3 body decayers '''

                particledecayproducts = self.Choosehadrondecay(particlecode) # choose one of the decay routes for the unstable particle
                
                
                decayedproductlist = self.Performhadrondecay(particlelist[i],particledecayproducts)
                
                templist += decayedproductlist
                numberdecays += 1
        
        return stableparticles
                    
                    
                    
    def Performhadrondecay(self,hadron,productlist):
        twobodyclass = Twobodydecay.twobodydecay4D()
        
        if len(productlist) == 2:
            newproducts = twobodyclass.twobodyparticledecay(hadron,productlist)
            return newproducts
            
        if len(productlist) == 3:
            newproducts = twobodyclass.threebodyparticledecay(hadron,productlist)
            return newproducts
                        
    def Checkforkaons(self,particle):
        if abs(particle.Getcode()) == 311:
            rng = random.getrandbits(1) # generate a number which is 0 or 1. Line here in case of method change to speed up code.
            if rng == 0:
                particle.SetID(130)
            if rng == 1:
                particle.SetID(310)
                
            
    def Decaymanylists(self,particlelists):
        finallist = []
        for i in range(len(particlelists)):
            finalstate = self.Decayparticlelist(particlelists[i])
            finallist.append(finalstate)
            
        return finallist
        
            
