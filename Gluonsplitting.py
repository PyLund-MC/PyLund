import Hepmcconverter
import numpy
import Constants
import PDG_ID
import Decaydictionary
import copy
import MRS
import scipy.optimize
import Particle
import sys
import random

def Splitgluons(lop): # takes a list of particles. Create another class to call this multiple times for ease.
	
	no_parts = len(lop)
	no_gluons = 0
	splitlist = [] # create a new list to put the completed list in.
	#splitlist.append(lop[0]) # Put the first gluon into the list.
	for i in range(no_parts):
		if lop[i].Checkgluon() == True:
			no_gluons += 1 # Get how many gluons that we have in the list. Can use this to get a finishing condition when they all split.

	index = 0 # create an index to keep track of location in the oringal particle list.

	if no_gluons == 0:
		return lop

	while no_gluons != 0:

		if lop[index].Checkgluon() == True:
			"print current particle in splitgluons function is a gluon"
			# split the gluon. append the new quarks to the list. increase the index, and subtract 1 from number of gluons.
			# Need to feed in the left, current gluon, and the right particle in the list.
			finalquarks = Getsplitquarks(splitlist[-1],lop[index],lop[index+1]) # the left side particle is the rightmost entry  in the finalparticle list.
			splitlist.append(finalquarks[0])
			splitlist.append(finalquarks[1])
			index +=1
			no_gluons -= 1
			continue

		else:
			splitlist.append(lop[index]) # If not a gluon - must be a quark. append it to the list and increase the index value
			index +=1
			continue

	splitlist.append(lop[-1]) # some bookkeeping. append the final antiquark in the chain to the list.
	

	return splitlist
		


def Gluonflavour():

	fdict = flavourdict() # fdict - flavour dictionary class
	flavours = [1,2,3,1103,2101,2103,2203,3101,3103,3201,3203,3303]
	probs = []

	for i in flavours:
		prob = fdict.Getprob(i)
		probs.append(prob)

	#print probs, "probs"
	no_probs = len(probs)
	#print no_probs, "noprobs"
	sum_probs = sum(probs)
	#print sum_probs, "sumprobs"
	normalised_probs = probs/sum_probs
	#print normalised_probs, "normed probs"

	rng = numpy.random.choice(no_probs,1,p=normalised_probs)

	chosenflavour = flavours[rng]

	return chosenflavour

def Selectflavour(left,dec_g,right): # selects a flavour based on criteria to keep the result physical! Contains veto conditions for production of various stuff
	
	leftflavour = left.Checkdiquark()
	rightflavour = right.Checkdiquark()
	mq = 0.6
	Egr = right.Getvector()[0]
	strange = [3,-3]

	while True:
		splitflavour = Gluonflavour()
		if leftflavour == True and splitflavour > 10:
			continue
		if rightflavour == True and splitflavour > 10: # if we try to produce adjacent diquarks we have a problem
			continue
			# if the adjacent are diquark and our chosen isnt then no problem!

		if splitflavour > 10: # we have a diquark produced
			rightvec4D = right.Getvector4D()
			gvec4d = dec_g.Getvector4D()
			inv = rightvec4D + gvec4d
			invmass = numpy.sqrt(inv * inv)
			if invmass < 3: # invariant mass is less than two for the decaying gluon and the right hand side particle. then veto diquark production
				continue
		
		if splitflavour in strange:
			if 2*mq > Egr: # additional condition, not sure if needed.
				continue

		



		return splitflavour
	
def Getsplitquarks(left,dec_g,right):
	#gets teh various functions and collates them and gives teh final [qq] state to append to the final list. -- the handler function
	n = 0
	while True:
		#print "getsplitquarks function stuck"
		n += 1
		if n >= 500:

			if left.Checkdiquark() == True:
				#print "diquark on the left, have 500 runs so cannot find a mathc. trying to replace left quark - We canny do that ."
				#print left.GetID()
				sys.exit()

			#print "gluon splitting has failed 500 times, reassigning the left ID to u/d"
			rng = bool(random.getrandbits(1))
			if left.GetID() < 0:
				if rng == 0:
					left.SetID(-1)
				if rng == 1:
					left.SetID(-2)

			if left.GetID() > 0:
				if rng == 0:
					left.SetID(1)
				if rng == 1:
					left.SetID(2)



			#sys.exit()

		flavour = Selectflavour(left,dec_g,right)
		z = Gluonz(left,dec_g,right,flavour)
		if z == False:
			#print "problem with the z value. too light a quark is getting through"
			continue
		
		finalstate = Getmh(left,dec_g,right,z,flavour)
		if finalstate == False:
			#print "reached finalstate, finalstate failed"
			n+=1
			continue

		else:
			return finalstate
			break

def Getmh(left,dec_g,right,z,splitflavour): # return False if conditions not met. Otherwise return normal vector output.

	Decaytables = Decaydictionary.Decays()
	gvec4Dx = dec_g.Getvector4D()
	gvec4D = copy.deepcopy(gvec4Dx)
	sgl_flavour = copy.deepcopy(splitflavour)

	if right.Checkgluon() == True: # particle on the right is a gluon. only have to consider the left hand side particle.
		#print "right particle is gluon, check that part"
		# Need the information of the quark / antiquark on the left.
		leftvector4D = left.Getvector4D()
		#print leftvector4D, "left vector"
		if left.Checkdiquark() == False: # left is a quark. Have to assign according to that idea
			if left.GetID() > 0: # if the left is the quark. need antiflavour
				if splitflavour > 10: # i.e. were producing a diquark
					sgl_flavour = copy.deepcopy(splitflavour) # we want the positive again.
				# UNLESS THE FLAVOUR CREATED IS A DIQUARK
				if splitflavour < 10: # i.e. we not producing a diquark
					sgl_flavour *= -1 # we want the negative quark on the left of the gluon split

			if left.GetID() < 0: # we have an anti quark on the left somehow??. Perhaps we've already had some diquark produciton :/
				if splitflavour > 10:
					sgl_flavour = copy.deepcopy(splitflavour) * -1

				if splitflavour < 10:
					sgl_flavour = copy.deepcopy(splitflavour)

		if left.Checkdiquark() == True:
			if left.GetID() < 0: # negative diquark
				sgl_flavour = copy.deepcopy(splitflavour) * -1
			if left.GetID() > 0: # positive diquark
				#print "we haev a positive diquark"
				#print left.GetID(), "left ID"
				sys.exit()
				sgl_flavour = copy.deepcopy(splitflavour)


		sgl_vec4D = copy.deepcopy(gvec4D)
		#print sgl_vec4D.Getvector()[0], "Energy of the gluon"
		sgl_vec4D *= z
		#print "scaled vector", sgl_vec4D

		# get the scaled momentum.
		sgl = Particle.particle4D(sgl_flavour, sgl_vec4D)



		# Create a temp string to test for FS1 hadron. I.e. if there is a small hadron that can be made.
		#print left.GetID(), "left ID"
		#print left.Getvector(), "left vector"
		#print sgl_flavour, "sgl flavour"
		#print sgl.Getvector(), "sgl vector"
		string = MRS.mrs4D(left,sgl)
		#print string.Getstringmass(), "string mass of the left,sgl pair"

		fs1hadron = Decaytables.GetFailsafe1hadron(string)
		#print fs1hadron, "fs1hadron in the right split code"

		if fs1hadron == 0:
			return False # condition failed. cant produce a hadron.


		fs1mass = PDG_ID.PDG_Type.mass(PDG_ID.PDG_Type(fs1hadron))
		#print fs1mass

		totalmom = sgl_vec4D + leftvector4D
		if totalmom*totalmom < (fs1mass*fs1mass): # Condition that the resultant qqbar pair is heavy enough to at least form a hadron.
			return False # Condition failed. Need to reloop again for different z.

		# Provided above conditions are satisfied, return the particles in correct sides based on nature of the left quark.
		# create the other particle
		sgr_flavour = copy.deepcopy(sgl_flavour) * -1
		sgr_vec4D = copy.deepcopy(gvec4D)
		sgr_vec4D *= (1-z)

		''' Create a new proxy particle to form a string with the right hand side to check cons stuff '''
		copyright = copy.deepcopy(right)

		newrightflavour = 1
		if abs(sgr_flavour) > 10: # we have a diquark
			if sgr_flavour > 0: # positive diquark
				newrightflavour = 1
			if sgr_flavour < 0:
				newrightflavour = -1

		if abs(sgr_flavour) < 10: # we ahve a regular quarks
			if sgr_flavour > 0: # positive quark
				newrightflavour = -1
			if sgr_flavour < 0:
				newrightflavour = 1


		copyright.SetID(newrightflavour)
		#print copyright.Getvector(), "vector of the right hand side particle"

		#print sgr_vec4D, "vector of the right hand side gluon after split. They must be collinear or smthing"


		sgr = Particle.particle4D(sgr_flavour,sgr_vec4D)

		
		string2 = MRS.mrs4D(sgr,copyright)
		#print string2.Getstringmass(), "stringmass of the second pseudo string"
		fs1hadron2 = Decaytables.GetFailsafe1hadron(string2)
		if fs1hadron2 == 0:
			return False # new conditional

		if abs(splitflavour) > 10: # we have a diquark produced
			rightvec4D = right.Getvector4D()
			sgr4D = sgr.Getvector4D()
			inv = rightvec4D + sgr4D
			invmass = numpy.sqrt(inv * inv)
			if invmass < 2.5: # invariant mass is less than two for the decaying gluon and the right hand side particle. then veto diquark production
				return False
		
		# Flavours already flipped from the above condition. Should always be correct relative to the left and eachother.
		return [sgl,sgr]

	if right.Checkgluon() == False: # so the right side is a quark. need to account for both sides having sufficient mass to form hadrons.
		#print "right particle is a quark, see that algorithm"

		if left.Checkdiquark() == False:
			#print "left is not diquark"
			#print splitflavour, "splitflavour in this party"
			if left.GetID() > 0: # if the left is the quark. need antiflavour
				if splitflavour > 10: # i.e. were producing a diquark
					sgl_flavour = copy.deepcopy(splitflavour) # we want the positive again.
				# UNLESS THE FLAVOUR CREATED IS A DIQUARK
				if splitflavour < 10: # i.e. we not producing a diquark
					sgl_flavour *= -1 # we want the negative quark on the left of the gluon split

			if left.GetID() < 0: # we have an anti quark on the left somehow??. Perhaps we've already had some diquark produciton :/
				#print "left ID is less than zero"
				if splitflavour > 10:
					"we've reached the left neg, g diquark part"
					sgl_flavour = copy.deepcopy(splitflavour) * -1
					#print left.GetID()
					#print sgl_flavour, "flavour of sgl"

				if splitflavour < 10:
					sgl_flavour = copy.deepcopy(splitflavour)
			
		if left.Checkdiquark() == True:
			if left.GetID() < 0: # negative diquark
				sgl_flavour = copy.deepcopy(splitflavour) * -1
			if left.GetID() > 0: # positive diquark
				sgl_flavour = copy.deepcopy(splitflavour)

		sgl_vec4D = copy.deepcopy(gvec4D)
		print sgl_vec4D, "sgl vec4D"
		print sgl_vec4D * sgl_vec4D
		sgl_vec4D *= z
		 
		sgr_vec4D = copy.deepcopy(gvec4D)
		print sgr_vec4D, "sgrvec4D"
		sgr_vec4D *= (1-z)
		sgr_flavour = copy.deepcopy(sgl_flavour) * -1 # sgr is the anti of the sgl
		#print right.GetID(), "right ID"
		#print left.GetID(), "left ID"
		#print sgr_flavour, "sgrflavour"
		#print sgl_flavour, "sglflavour"

		sgl = Particle.particle4D(sgl_flavour,sgl_vec4D)
		sgr = Particle.particle4D(sgr_flavour,sgr_vec4D)

		string1 = MRS.mrs4D(left,sgl)
		#print string1.Getstringcontent(), "stringcontent of string 1 in gluon splitting rigt hand is quark"
		string2 = MRS.mrs4D(sgr,right)
		#print string2.Getstringcontent(), "stringcontent of string 2 in gs right hand is quark"
		#print string1.Getstringmass(), string2.Getstringmass(), "string masses of 1,2 respectuively"

		fs1hadron1 = Decaytables.GetFailsafe1hadron(string1)
		#print fs1hadron1
		if fs1hadron1 == 0:
			return False
		fs1hadron2 = Decaytables.GetFailsafe1hadron(string2)

		#print fs1hadron2
		if fs1hadron2 == 0:
			return False

		return [sgl,sgr]













	





def Gluonz(left,dec_g,right,qflavour):# function that finds the z value for the gluon splitting using MC techniques.
	''' The function looks at the adjacent quark (and anti quark if final gluon in chain) \n
	and generates a z value to split the gluon into two quarks that subsequently form hadrons that \n
	are mass possible '''

	''' Need to find the variables that constrain the flavour and z values. mass, nearby stuff etc. '''
	
	qid = PDG_ID.PDG_Type(qflavour)
	if qid.isDiquark() == True:
		mq = PDG_ID.PDG_Type.mass(qid) # mass of the quark. zero for the light quarks, non zero for the diquarks.

	lights = [1,2,-1,-2]
	strange = [3,-3]

	if qid.code() in lights:
		mq = 0.3

	if qid.code() in strange:
		mq = 0.45

	''' gluon information '''
	gmom = dec_g.Getvector()
	gvec4D = dec_g.Getvector4D()
	Eg = gmom[0] # energy of the gluon.

	#print mq, "quark mass"
	#print Eg, "gluon energy"
	

	''' find some properties that are needed for the split '''
	z_low = mq/Eg
	z_high = 1 - (mq/Eg)

	if z_low > z_high:
		return False

	''' Get z value. Use a constant as an overapproximator '''
	def f(z):
		return 0.5*(z*z + (1-z)*(1-z))
	z_max = scipy.optimize.fminbound(lambda z: -f(z), z_low, z_high)
	g_z = 0.5*(z_max*z_max + (1-z_max)*(1-z_max))
	Constant = g_z

	no_tries = 0
	while True:
		#print "stuck in the z generation formula"
		no_tries += 1
		r_1 = numpy.random.uniform(0,1) #(z_low,z_high)
		r_2 = numpy.random.uniform(0,1) #(z_low,z_high)
		z_test = (r_1*(z_high-z_low) + z_low)

		f_z_test = 0.5*(z_test*z_test + (1-z_test)*(1-z_test))
		g_z_test = Constant
		testvalue = f_z_test / g_z_test

		if r_2 <= testvalue:
			#put some code in here to test if the z value leaves kinematically possible hadrons / strings.

			if z_test * Eg < mq:
				continue
			if (1-z_test) * Eg < mq:
				continue
			#print z_test, "testvale of z obtained"

			return z_test # might ahve to return some additional values
			break
		else:
			continue




class flavourdict(object):
	 #class that just containing misc stuff with flavour choice with gluon splitting 

	def __init__(self):
		self.__match = {

			1: pair(1,exp(1),-1),
			2: pair(1,exp(2),-2),
			3: pair(Constants.strange_suppression_factor(),exp(3),-3),
			1103: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor(),exp(1103),-1103),
			2101: pair(Constants.diquark_scaling_factor(),exp(2101),-2101),
			2103: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor(),exp(2103),-2103),
			2203: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor(),exp(2203),-2203),
			3101: pair(Constants.diquark_scaling_factor() * Constants.strange_diquark_scaling_factor(),exp(3101),-3101),
			3103: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor(),exp(3103),-3103),
			3201: pair(Constants.diquark_scaling_factor() * Constants.strange_diquark_scaling_factor(),exp(3201),-3201),
			3203: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor(),exp(3203),-3203),
			3303: pair(Constants.diquark_scaling_factor() * Constants.vector_diquark_scaling_factor() * Constants.strange_diquark_scaling_factor(),exp(3303),-3303)

			}

	def Getprob(self,ID):
		if self.__match.has_key(ID):
			supp = self.__match[ID].Getsf()
			mf = self.__match[ID].Getmf()
			prob = supp * mf


			return prob

		return "incorrect ID input"

	def Getanti(self,ID):
		if self.__match.has_key(ID):
			anti = self.__match(ID).Getanti()
			return anti
		return "incorrect ID input"

class pair(object):
	def __init__(self,sf,mf,anti): # supp factors, masss factors, anti particle
		self.__sf = sf
		self.__mf = mf
		self.__anti = anti


	def Getanti(self):
		return self.__anti

	def Getsf(self):
		return self.__sf

	def Getmf(self):
		return self.__mf



def exp(number):
	pdgid = PDG_ID.PDG_Type(number)
	mass = PDG_ID.PDG_Type.mass(pdgid)
	sigma = Constants.kperpsigma()
	value = numpy.exp(-1*mass*mass/(sigma*sigma))
	return value
