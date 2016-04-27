''' Decay strings class '''

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
import Decaydictionary
import Twobodydecay
#Decaytables = Decaydictionary.Decays()


class stringdecay4D(object):
    def __init__(self,mrs):
        #assert isinstance(mrs,MRS.mrs4D) # throwing an error for now, have a look at in the future
        self.__classname        = "I'm the string decay class!"
        self.__description      = "I contain several functions, used for the boosting of strings, their decay, and the associated kinematics!"
        self.__mrs              = mrs
        ''' Quark Properties '''
        self.__q                = mrs.Getq()
        self.__qID              = mrs.GetqID()
        self.__qPDG             = PDG_ID.PDG_Type(self.__qID)
        self.__qmom             = mrs.Getqvector()
        self.__qvector4D        = mrs.Getqvector4D()
        ''' Antiquark Properties '''
        self.__qbar             = mrs.Getqbar()
        self.__qbarID           = mrs.GetqbarID()
        self.__qbarPDG          = PDG_ID.PDG_Type(self.__qbarID)
        self.__qbarmom          = mrs.Getqbarvector()
        self.__qbarvector4D     = mrs.Getqbarvector4D()
        ''' Total string Properties '''
        self.__totalvector4D    = mrs.Gettotalvector4D()
        self.__totalmom         = mrs.Gettotalvector()
        
        ''' create the lorentz boosts '''
        self.__lorentzboost = Boosts.boost4D(self.__totalvector4D) # uses the total cmf vector4D of the string
        
        ''' create the rotational boost '''
        self.__qmomboosted = self.__lorentzboost*self.__qvector4D
        self.__qbarmomboosted = self.__lorentzboost*self.__qbarvector4D
        self.__rotateboost = Boosts.rotate4D(self.__qmomboosted)
        
        
    def Booststringforward(self,string):
        ''' Seperate out the quark and anti quark '''
        quark = string.Getq()
        antiquark = string.Getqbar()
        quarkvector4D = quark.Getvector4D()
        antiquarkvector4D = antiquark.Getvector4D()
        ''' Perform the boost - lorentz and then rotate afterwards '''
        quarkvector4Dboosted = self.__lorentzboost*quarkvector4D
        antiquarkvector4Dboosted = self.__lorentzboost*antiquarkvector4D
        quark_z_aligned = self.__rotateboost.Rotateforwards(quarkvector4Dboosted)   
        antiquark_z_aligned = self.__rotateboost.Rotateforwards(antiquarkvector4Dboosted)
        ''' Set the new quark momentums and create the new string '''
        quark.Setmom(quark_z_aligned)
        antiquark.Setmom(antiquark_z_aligned)
        finalstring = MRS.mrs4D(quark,antiquark)
        return finalstring
        
    def Booststringback(self,boostedstring):
        ''' Seperate the quarks and antiquarks from the given string '''
        quark1 = boostedstring.Getq()
        antiquark1 = boostedstring.Getqbar()
        quark1vector4D = quark1.Getvector4D()
        antiquark1vector4D = antiquark1.Getvector4D()
        ''' Perform the boosts. Unrotate and then boost back '''
        quark1_rotated = self.__rotateboost.Rotatebackwards(quark1vector4D)
        antiquark1_rotated = self.__rotateboost.Rotatebackwards(antiquark1vector4D)
        quark1_reverseboosted = self.__lorentzboost/quark1_rotated
        antiquark1_reverseboosted = self.__lorentzboost/antiquark1_rotated
        ''' Set the momentum of the quarks, and create a new string with those quarks '''
        quark1.Setmom(quark1_reverseboosted)
        antiquark1.Setmom(antiquark1_reverseboosted)
        finalstring = MRS.mrs4D(quark1,antiquark1)
        return finalstring
        
        
    def Decaystring(self,decayingstring,decayingquarknumber,hadronproducedID): #class that takes an (unboosted) mrs as an arguement and decays it, producing a particle and a new iteration of string.
        ''' NOTE - DECAYING SIDE ONLY CARES IF THE NUMBER IS POSITIVE OR NEGATIVE - I.E. PDG CODE OF QUARK / ANTIQUARK REPSECTIVELY '''
        #print "decaystring class reached"
        Decaytables = Decaydictionary.Decays()
        twobodyclass = Twobodydecay.twobodydecay4D()
        hadroncodetest = PDG_ID.PDG_Type(hadronproducedID)
        
        if hadroncodetest.isBaryon() is True:
            if self.__qPDG.isQuark() or self.__qbarPDG.isQuark() is True:
                baryonproducts = self.Baryonproduction(hadronproducedID,decayingquarknumber)
                if baryonproducts[1] == 0:
                    stringcontent = [baryonproducts[0],self.__qbarID]
                    #print stringcontent, "IS THERE SOMETHING WRONG HERE"
                if baryonproducts[1] == 1:
                    stringcontent = [self.__qID, baryonproducts[0]]
                    #print stringcontent, "HOW ABOUT HERE, IS THERE ERROR"
        
        ''' If the produced particle is a diquark i.e. producing diquarks at the ends of the string. just need to redistribute momentum via a two body decay '''
        if hadroncodetest.isDiquark() is True: # i.e. we're producing two diquarks, one at each end of the string
            recievingdiquark = self.Finddiquarkstringcontent(decayingquarknumber,hadronproducedID)
            #print recievingdiquark, "RECIEVINGDIQUARK"
            newstring = twobodyclass.Twobodydiquarks(decayingstring,hadronproducedID,recievingdiquark) # produce a new string, hadronproduced ID here is a DIQUARK.
            adc = Availabledecays.availabledecays4D(newstring)# AVAILABLE DECAY CLASS
            newproducts = adc.Finddecay()
            
            ''' Assign the new values to the parameters in the code. now producing a baryon from a diquark end '''
            decayingquarknumber = newproducts[0]
            #print decayingquarknumber, "DECAYINGQUARK NUMBER IN THE BARYONPROD"
            hadronproducedID = newproducts[1]   
            #print hadronproducedID, "HADRONPRODUCED ID IN THE BAREEREREREOWIJIW3ERJIOWEJ"     
            #print PDG_ID.PDG_Type.mass(PDG_ID.PDG_Type(hadronproducedID)), "HADRON MASS IN THE BARREEEEEE"
            #print "SUCESSFULLY DETERMINED A MASS"
            
            stringcontent = self.Getbaryonstringcontent(hadronproducedID,decayingquarknumber,recievingdiquark) # takes the (baryon, diquark) IDs. that have just been found
            decayingstring = copy.deepcopy(newstring)
            #print stringcontent, "stringcontent from the module"
            
            
            
        ''' If the produced particle is a meson, leaves behind a quark remnant at the end of the string '''
        
        if hadroncodetest.isMeson() is True:
            stringcontent = self.Findstringcontent(decayingquarknumber,hadronproducedID) # finds [quark, antiquark] of the string
            #print stringcontent, "DO WE GOT A PROBLEM HERE"
         
        ''' Get out the string information (ie the quark information) '''
        boostedstring = self.Booststringforward(decayingstring)
        quark1 = boostedstring.Getq()
        antiquark1 = boostedstring.Getqbar()
        quark1vector4D = quark1.Getvector4D()
        antiquark1vector4D = antiquark1.Getvector4D()
       
        ''' define particle produced values ''' # Initialise the particle for manipulation later. (hadron is object k)
        hadronmomentum = Fourvector.vector4D(1.,1.,1.,1.) # placeholder (used for class creation), is reset later
        hadron = Particle.particle4D(hadronproducedID,hadronmomentum)
        hadronmass = hadron.Getmass()
        
        '''Create p+, p-'''
        #print antiquark1vector4D.Getvector(), "#PRINT VECOTRS"
        #print quark1vector4D.Getvector()
        if quark1vector4D[3] < 0:
            p_plus = quark1vector4D[0] - quark1vector4D[3]
            p_minus = antiquark1vector4D[0] + antiquark1vector4D[3]  
            
        if quark1vector4D[3] > 0:
            p_plus = quark1vector4D[0] + quark1vector4D[3]
            p_minus = antiquark1vector4D[0] - antiquark1vector4D[3]  
        
        ''' Find the quark content of the string '''
        
            
        
        
        #print stringcontent[0], "stringcontent 0 "
        #print stringcontent[1], "stringcontent 1 "           
        
        ''' Try 100 hundred times to get a successful kperp and z value for the current string decay '''
        n = 0
        while True:
            n = n+1
            ##print "Are we stuck"
            #print "stringmass", boostedstring.Getstringmass()
            #print "hadronmass", hadronmass
            #print antiquark1vector4D.Getvector(), "#PRINT VECOTRS"
            #print quark1vector4D.Getvector()
        
            ''' Get k perp value '''
            kperp = self.Getkperphitmiss(boostedstring) # 0.305
            #print kperp, "kperp"
        
            ''' Get z value '''
            z_plus = self.Getzhitmiss(boostedstring,kperp,hadronmass) # 0.24
            z_minus = 1 - z_plus
            #print z_plus, "z_plus"
        
            ''' Create the k+- values '''
            k_plus = z_plus*p_plus
            k_minus = (hadronmass*hadronmass + kperp*kperp) / (z_plus*p_plus)
            ##print k_plus, "k_plus"
            ##print k_minus, "k_minus"
        
            ''' Create the angles for future resolving '''
            ''' phi is isotropic, and taken as angle to the quark. pi is added to phi to get the hadron angle '''
            q_phi = 2*numpy.pi*random.random()
            k_phi = q_phi + numpy.pi
            qnew_cosphi = numpy.cos(q_phi)
            qnew_sinphi = numpy.sin(q_phi)
            k_cosphi = numpy.cos(k_phi)
            k_sinphi = numpy.sin(k_phi)
             
            ''' Create the final vectors '''
            qnew_plus = z_minus * p_plus
            #print z_minus
            #print p_plus, "P PLUS"
            qnew_minus = (kperp*kperp) / (z_minus * p_plus)
        
            ''' Arrange arrays of [p+, p_, pcostheta, psintheta] '''
            qnew = [qnew_plus, qnew_minus, kperp*qnew_cosphi, kperp*qnew_sinphi]
            #print qnew, "QNEW"
            knew = [k_plus, k_minus, kperp*k_cosphi, kperp*k_sinphi]
            #print knew, "KNEW"
            antiqnew = [0,p_minus - k_minus - qnew_minus, 0, 0]
            #print antiqnew, "ANTIQNEW"
            
            ''' From the arrays constructed above, create the vectors E,p '''
            
            ''' Start with the particle '''
            k_E = (k_plus + k_minus) /2
            k_px = knew[2]
            k_py = knew[3]
            k_pz = (k_plus - k_minus) /2
            kmom = Fourvector.vector4D(k_E,k_px,k_py,k_pz)
            
            ''' And now for the quark in the string '''
            q_E = (qnew[0]+qnew[1])/2
            q_px = qnew[2]
            q_py = qnew[3]
            q_pz = (qnew[0]-qnew[1])/2
            qmom = Fourvector.vector4D(q_E,q_px,q_py,q_pz)
            
            ''' And now for the antiquark in the string '''
            antiq_E = (antiqnew[0] + antiqnew[1]) /2
            antiq_px = antiqnew[2]
            antiq_py = antiqnew[3]
            antiq_pz = (antiqnew[0] - antiqnew[1]) /2
            antiqmom = Fourvector.vector4D(antiq_E,antiq_px,antiq_py,antiq_pz) 
            
            ''' Testing the dynamics of the system - momentum should be conserved ''' 
            totalmom = antiqmom+qmom+kmom
            ##print p_plus*p_minus, "total cmf momentum of the system originally, i.e. 2E"
            ##print totalmom*totalmom, "Invariant mass squared of the final system"   
            
            ''' Now need to boost the system back into the lab frame '''
            ''' Adjust the quarks and create a new string, and boost to lab frame '''
            #print stringcontent, "STRINGCONTENT"
            quark2 = Particle.particle4D(float(stringcontent[0]),quark1.Getvector4D())
            #print hadronproducedID, "HADRON ID"
            ##print stringcontent[1], "ANTI QUARK NUMBER"

            antiquark2 = Particle.particle4D(float(stringcontent[1]),antiquark1.Getvector4D())
            quark2.Setmom(qmom)
            antiquark2.Setmom(antiqmom)
            quark1.Setmom(qmom)
            antiquark1.Setmom(antiqmom)
            newstring = MRS.mrs4D(quark2,antiquark2)
            newstring_labframe = self.Booststringback(newstring) # Final string, in the lab frame
            #print newstring_labframe.Getstringcontent(), "LABFRAME CONTENT"

            #testthequark = newstring_labframe.Getq()
            #print testthequark.Getinvariantmass(), "invariant mass of the quark in the string"
            
            ''' Now need to boost then adjust the produced hadron '''
            kmomrotated = self.__rotateboost.Rotatebackwards(kmom) # rotate out or rotated frame
            kmomboosted = self.__lorentzboost/kmomrotated          # boost into the lab frame (original frame)
            hadron.Setmom(kmomboosted)
           
            ##print Decaytables.Findminandmaxpairmass(newstring_labframe)[1], "testing the condition for failure here"
            ##print newstring_labframe.Getstringmass(), "stringmass in the new labframe"    
            ##print Decaytables.Stringdecaysmallesthadron(newstring_labframe), "result of the decaytables testing shit, searching for error"
            ##print newstring_labframe.Getstringcontent(), "stringcontent"
            
            #if n > 3:
                #sys.exit()
            ''' Check validity of the decay (cons of energy) '''
            if qnew_minus + k_minus >= quark1vector4D[0] + antiquark1vector4D[0]:
                #print qnew_minus
                #print k_minus
                #print quark1vector4D[0]
                #print antiquark1vector4D[0]
                #print "CONDITION 1 FAILED"
                continue
                
            if quark2.Checkquark() is True and antiquark2.Checkquark() is True:
                meson = Decaytables.GetFailsafe1hadron(newstring_labframe)
                #print meson, "STRINGDECAY ERROR MESSAGE. THIS IS WHAT TRYING TO GET MASS OF"
                #mass = Decaytables.Getmass(meson)
                if meson == 0:
                    continue
                mass = Decaytables.Getmass(meson)
                if newstring_labframe.Getstringmass() <= mass:
                    #print "CONDITION 2 FAILED"
                    continue
                break
            if Decaytables.GetFailsafe1hadron(newstring_labframe) == 0:  # if no possible stringcollapse avaiable, try and produce something with less mom frac
                #print "CONDITION 3 FAILED"
                continue # previous findminandmaxpairmass(string)[0]
            else:
                break
            
        kaoncheck = self.Checkforkaons(hadronproducedID)
        hadron.SetID(kaoncheck)
        #print "STRING SUCCESSFULLY DECAYED"
        return [newstring_labframe, hadron] # returns the new string and a hadron, both boosted into the lab frame
        
        
    def Getzhitmiss(self,string,kperp,mhad):
        ''' function used to generate a sensible z value based on a hit or miss rejection basis '''
        ''' Get the constants to be used '''
        M = string.Getstringmass()
        m_T = numpy.sqrt(kperp*kperp + mhad*mhad) #
        
        z_low = (mhad*mhad + kperp*kperp) / (M*M) # M is the string mass TRYING THE SYTEM WITH TRANS MASS
        z_high = 1 - kperp*kperp/(M*M)
        a = Constants.aLund()
        b = Constants.bLund()
        
        ''' Calculate the f_max or g(x) '''
        def f(z):
            return (1/z)*((1-z)**a)*numpy.exp(-b*m_T*m_T/z)
        z_max = scipy.optimize.fminbound(lambda z: -f(z), z_low, z_high)
        z_max = z_max*1.01
        g_z = (1/z_max)*((1-z_max)**a)*numpy.exp(-b*m_T*m_T/z_max) # this is taken as g(x) (its a constant)
        Constant = g_z #  a reminder that g(z) here used is a constant. The maximum value that f(z) can take
        G_z_max = z_high * g_z
        G_z_min = z_low * g_z
         
        ''' create the iterative loop to find a suitable z value, based on hit or miss principle '''
        while True:
            rand_no_1 = numpy.random.uniform(0,1) #(z_low,z_high) # this number used to calculate X test
            rand_no_2 = numpy.random.uniform(0,1) #(z_low,z_high) # this number used to accept / reject with given probability, as compared to a test value
            z_test = (rand_no_1*(z_high-z_low) + z_low) #/ Constant
            
            f_z_test = (1/z_test)*((1-z_test)**a)*numpy.exp(-b*m_T*m_T/z_test)
            g_z_test = Constant # constant as listed above
            testvalue = f_z_test / g_z_test
        
            if rand_no_2 <= testvalue:
                return z_test
                break
            else:
                continue
        
    
    def Getkperphitmiss(self,string): # to be used within the string decay class
        ''' function used to generate a sensible kperp value based on a hit or miss rejection basis '''
        ''' Find the suitable range (min - max) of kperp values based on the given string '''
        sigma = Constants.kperpsigma()
        g_sigma = Constants.gsigma() # this value is for g(x) such that g(x) > f(x). If edited, needs to be recalculated (c.f. finding a suitable.... generation)
        #g_sigma = Constants.kperpsigma()
        B = 1/(g_sigma*g_sigma)
        A = Constants.g_kperpA()
        M = string.Getstringmass()
        kperpmin = 0
        kperpmax = M/2 # where M is the invariant mass of the string decaying
        # g_kperp = Aexp(-kperp/(gsigma*gsigma) = Aexp(-Bkperp)            #numpy.exp(-M*M/(sigma1*sigma1)) - this factor can be omited as it cancels when creating testvalue (its essentially a constant)
        G_kperp_max = -(A/B)*numpy.exp(-B*kperpmax)
        G_kperp_min = -(A/B)*numpy.exp(-B*kperpmin)
        
        while True:      
            rand_1 = numpy.random.uniform(0,1) #(kperpmin,kperpmax)
            rand_2 = numpy.random.uniform(0,1) #(kperpmin,kperpmax)
            kperp_test = -(1/B)*numpy.log(rand_1*(numpy.exp(-kperpmax*B)-1) + 1) # working inverse function
            #kperp_test = -(1/B)*numpy.log(rand_1*(1-numpy.exp(kperpmax*B)) - 1) # bfroken inverse function
            #kperp_test = -(1/B)*numpy.log(rand_1 - 1 +numpy.exp(-kperpmax*B)) # dans inverse function

            f_kperptest = numpy.exp(-kperp_test*kperp_test/(sigma*sigma))     # system designed such that gx is always greater than fx
            g_kperptest = A*numpy.exp(-B*kperp_test)
            testvalue = f_kperptest/g_kperptest
        
            if rand_2 <= testvalue:
                return kperp_test
                break
            else:
                continue
          
          
    def Findstringcontent(self,quark_no,hadron_no): # this submodule takes the deacying quark number and meson and gives back the quark and antiquark values for that decay.
        ''' FUNCTION USED TO FIND STRING COMP BASED ON MESON ENDS - CALLED THE DECAYSTRING FUNCTION '''
        
        ''' Find hadron quark / antiquark '''
        Decaytables = Decaydictionary.Decays()
        mesoncontent = Decaytables.Getmesoncontent(quark_no,hadron_no) # returns [meson q, meson qbar]
        #print mesoncontent, "mesoncontent"
        
        if quark_no == self.__qID:   # real q of string decays - string is string qbar + meson q
            antiquark = self.__qbarID
            if quark_no < 0: # negative quark number
                return [-1 * mesoncontent[0], antiquark]
            if quark_no > 0: # positive quark number
                return [-1 * mesoncontent[1], antiquark]
        
        if quark_no == self.__qbarID:    # qbar of string decays - string is string q + meson qbar
            quark = self.__qID
            if quark_no < 0: # negative number
                return [quark, -1* mesoncontent[0]]
            if quark_no > 0:
                return [quark, -1 * mesoncontent[1]]
        
        #print "STRINGCONTENT MODULE WORKING INCORRECTLY"
        sys.exit()
            
            
    def Finddiquarkstringcontent(self,quark_no,diquark_no):
        ''' FUNCTION USED TO FIND THE FLAVOUR CONTENT OF THE DIQUARKS, AND ASSIGN THE PDGIDS PRODUCED (EACH ONE IS UNIQUE) '''
        
        if quark_no > 0:
            recievingflavour = -1
            
        if quark_no < 0:
            recievingflavour = 1 # check if reciever is positive or negative (real / anti)
                      

        
        absflavourstring = abs(diquark_no)
        flavourstring = str(absflavourstring)
        flavourstring = flavourstring[0:4]

        spin = flavourstring[-1]
        flavour = flavourstring[-4:-2]


        if abs(float(quark_no)) == abs(float(flavour[0])):
            vacumnflavour = flavour[1] # build absolute then invert if needed
            absrecievingflavour = abs(recievingflavour)
            first = max(float(vacumnflavour),float(absrecievingflavour))
            first = float(first)
            first = int(first)
            first = str(first)
            second = min(float(vacumnflavour),float(absrecievingflavour))
            second = float(second)
            second = int(second)
            second = str(second)
            zero = str(0)
            spin = str(spin)
            code = first+second+zero+spin
            
            code = float(code)
            code = round(code,0)
            code = int(code)
            if recievingflavour < 0:
                code = code *-1
                
                        
                
                
            return code
            
        if abs(float(quark_no)) == abs(float(flavour[1])):
            vacumnflavour = flavour[0]
            absrecievingflavour = abs(recievingflavour)
            first = max(float(vacumnflavour),float(absrecievingflavour))
            first = float(first)
            first = int(first)
            first = str(first)
            second = min(float(vacumnflavour),float(absrecievingflavour))
            second = float(second)
            second = int(second)
            second = str(second)
            zero = str(0)
            spin = str(spin)
            code = first+second+zero+spin

            code = float(code)
            code = round(code,0)
            code = int(code)
            if recievingflavour < 0:
                code = code *-1
            return code
            
        #print "diquark flavour association not working correctly"
        return None
        
        
    def Getbaryonstringcontent(self,baryonID,diquarkID,recievingdiquark):
        Decaytables = Decaydictionary.Decays()
        
        leftover = Decaytables.Getbaryonleftover(baryonID,diquarkID)

        if recievingdiquark < 0:
            # Get the complete anti string
            leftover = leftover *-1
            self.__qID = leftover
            self.__qbarID = recievingdiquark
            return [leftover, recievingdiquark]
            
        if recievingdiquark > 0:
            # get the complete real string
            self.__qID = recievingdiquark
            self.__qbarID = leftover
            return [recievingdiquark, leftover]
            
        
        
        
    def Checkforkaons(self,code):
        if abs(code) == 311:
            rng = random.getrandbits(1) # generate a number which is 0 or 1. Line here in case of method change to speed up code.
            if rng == 0:
                return 130
            if rng == 1:
                return 310
                
        return code
                
  
    def Baryonproduction(self,baryonID,diquarkID):
        Decaytables = Decaydictionary.Decays()
        leftover = Decaytables.Getbaryonleftover(baryonID,diquarkID)
        if self.__qPDG.isDiquark() is True:
            if diquarkID > 0:
                return leftover*-1, 0
            if diquarkID < 0:
                return leftover, 0
                
        if self.__qbarPDG.isDiquark() is True:
            if diquarkID > 0:
                return leftover*-1, 1
            if diquarkID < 0:
                return leftover,1 