''' String former - forms strings from the adjacent quarks. '''
import MRS
import Decaydictionary
import sys

def Formstrings(listofquarks):
	Decaytables = Decaydictionary.Decays()
	stringlist = []
	#print len(listofquarks), "number of quarks in the stringformer module"
	for i in range(len(listofquarks)):
		if i % 2 == 0:
			#print i, "ivalue"
			string = MRS.mrs4D(listofquarks[i],listofquarks[i+1])
			#print string.Getstringmass(), "stringmass"
			#print string.Getstringcontent()
			test = Decaytables.GetFailsafe1hadron(string), "3"

			if test == 0:
				print i
				#rint "no producible hadron from this string, something has gone wrong"
				for j in range(len(listofquarks)):
					print listofquarks[j].GetID()
					print listofquarks[j].Getvector()
				sys.exit()
			stringlist.append(string)
			
	#sys.exit()
	return stringlist




