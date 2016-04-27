#! /usr/bin/env python

import sys

import numpy
import Fourvector
import MRS
#import Twobodydecay
import Particle
#import Boosts
import Decaydictionary
import Fullstringdecay
import Hadrondecay
import matplotlib.pyplot as pyplot
import Hepmcconverter




def Addlists(list1,list2,list3,weight1,weight2,weight3): # weights in decimal percentages
    totalweight = weight1+weight2+weight3
    relative1 = weight1/totalweight
    relative2 = weight2/totalweight
    relative3 = weight3/totalweight
    list4 = []
    for i in range(len(list1)):
        list4.append(list1[i]*relative1 + list2[i]*relative2 + list3[i]*relative3)
    return list4


Decaytables = Decaydictionary.Decays()

''' Create some momentum 4 vectors '''
print "Code beginning"

mom1 = Fourvector.vector4D(numpy.sqrt(41),3.,4.,4.)
mom2 = Fourvector.vector4D(numpy.sqrt(45),2.,4.,5.)
mom3 = Fourvector.vector4D(numpy.sqrt(101),9,4.,2.) 
mom4 = Fourvector.vector4D(6,2,4.,4.)
mom5 = Fourvector.vector4D(6,-2.,-4.,-4.)

''' momenta for 91.2 GeV '''


mom10 = Fourvector.vector4D(91.2/2,35.,20.,numpy.sqrt(454.36))
mom11 = Fourvector.vector4D(91.2/2,-35.,-20.,-numpy.sqrt(454.36))

''' Create some quarks '''

dquark = Particle.particle4D(1,mom10)
uquark = Particle.particle4D(2,mom10)
squark = Particle.particle4D(3,mom10)
cquark = Particle.particle4D(4,mom10)
bquark = Particle.particle4D(5,mom10)
tquark = Particle.particle4D(6,mom10)

''' Create some antiquarks '''

antidquark = Particle.particle4D(-1,mom11)
antiuquark = Particle.particle4D(-2,mom11)
antisquark = Particle.particle4D(-3,mom11)
anticquark = Particle.particle4D(-4,mom11)
antibquark = Particle.particle4D(-5,mom11)
antitquark = Particle.particle4D(-6,mom11)

''' Create some strings - string1 '''

string1 = MRS.mrs4D(dquark,antidquark)
stringmass = string1.Getstringmass()

fulldecayclass = Fullstringdecay.fulldecay4D()
hadrondecayclass = Hadrondecay.hadrondecay4D()
No_runs = 200

''' string2 '''

string2 = MRS.mrs4D(uquark,antiuquark)
stringmass2 = string2.Getstringmass()

''' string3 '''

string3 = MRS.mrs4D(squark,antisquark)


''' PERFORM THE STRING DECAY AND CREATE THE PRE / POST DECAY PARTICLE LISTS '''

''' STRING 1 '''

#Hadronlists1 = fulldecayclass.Performmanydecays(string1,No_runs)
#print "run 1 complete"

#Hadronlistspostdecay1 = hadrondecayclass.Decaymanylists(Hadronlists1)


''' STRING 2 '''

# Hadronlists2 = fulldecayclass.Performmanydecays(string2,No_runs)
# print "run 2 complete"
# #Hadronlistspostdecay2 = hadrondecayclass.Decaymanylists(Hadronlists2)

# ''' STRING 3 '''

# Hadronlists3 = fulldecayclass.Performmanydecays(string3,No_runs)
# print "run 3 complete"
#Hadronlistspostdecay3 = hadrondecayclass.Decaymanylists(Hadronlists3)


''' PLOTTING CODE '''

ddbarweight = 0.156
uubarweight = 0.116
ssbarweight = 0.156

pyplot.figure()
productionlist = []
''' Add the mesons to the code '''
productionlist += [111,211,-211,113,213,-213,223,221,331,310,130,321,-321,323,-323,313,-313,333]
''' Add the baryons to the code '''
productionlist += [2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334]
''' Add the antibaryons to the code '''
productionlist += [-2212,-2112,-2224,-2214,-2114,-1114,-3122,-3222,-3212,-3112,-3224,-3214,-3114,-3322,-3312,-3324,-3314,-3334]
''' Add the leptons to the code '''
productionlist += [22,11,-11,12,-12,13,-13,14,-14] # list of the producible particlies in my code. give the number of bars
noproducablehadrons = len(productionlist)
ind = numpy.arange(noproducablehadrons)
barwidth = 0.15

''' GET HADRON PERCENTAGES PRE DECAY '''

#hadronabundance1 = fulldecayclass.Gethadronabsoluteabundance(Hadronlists1,productionlist)
# hadronabundance2 = fulldecayclass.Gethadronabsoluteabundance(Hadronlists2,productionlist)
# hadronabundance3 = fulldecayclass.Gethadronabsoluteabundance(Hadronlists3,productionlist)



# totalweightedabundance = Addlists(hadronabundance1,hadronabundance2,hadronabundance3,ddbarweight,uubarweight,ssbarweight)




''' GET HADRON PERCENTAGES POST DECAY '''
'''
hadronabundance1postdecay = fulldecayclass.Gethadronabsoluteabundance(Hadronlistspostdecay1,productionlist)
hadronabundance2postdecay = fulldecayclass.Gethadronabsoluteabundance(Hadronlistspostdecay2,productionlist)
hadronabundance3postdecay = fulldecayclass.Gethadronabsoluteabundance(Hadronlistspostdecay3,productionlist)
'''



#y_upper = max(hadronabundance1) + max(hadronabundance1)*0.2
''' GET HADRON NAMES '''


#productionlistnames = fulldecayclass.Gethadronnames(productionlist)
# print totalweightedabundance
# print productionlistnames




# ''' AXIS SETTING '''

# fig, ax = pyplot.subplots()
# #rects1 = ax.bar(ind,hadronabundance1,barwidth, color='r')
# #rects2 = ax.bar(ind + 2*barwidth, hadronabundance2, barwidth, color='b')
# #rects3 = ax.bar(ind + 3*barwidth, hadronabundance3, barwidth, color='g')
# rects4 = ax.bar(ind, totalweightedabundance, barwidth, color='r')

# ax.set_xticks(ind+barwidth)
# ax.set_xticklabels(productionlistnames)
# ax.get_yaxis().set_visible(False)

# ''' AUTOLABEL FUNCTION '''

# def autolabel(rects):
    # # attach some text labels
    # for rect in rects:
        # height = rect.get_height()
        # #if height < 1:
            # #ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%d' % "<1", ha='center', va='bottom')
        # ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%d' % height, ha='center', va='bottom')

# #autolabel(rects1)
# #autolabel(rects2)
# #autolabel(rects3)
# autolabel(rects4)


# ''' PLOTTING CRITERIA / COSMETICS '''

# pyplot.title("10000 runs of (ddbar) string decay: Initial stringmass = 8.97 GeV")
# pyplot.xlabel("Particle produced")
# pyplot.ylabel("Percentage Frequency")
# pyplot.axis([-1,noproducablehadrons+1,0,y_upper])
# pyplot.show()


''' Test importing the PYTHIA data and running that through. see how we get along. '''

fulllistofstrings = Hepmcconverter.Convertfromhepmc() # get the list of string lists.
#print fulllistofstrings
fullfinalhadrons = [] # should be a list of lists per event.

for i in range(len(fulllistofstrings)): # for each list of strings
    if i %100 == 0:
        print "Event", i, "hadronisation complete"
    eventhadronlist = [] # create a list for the hadrons from that list
    for n in range(len(fulllistofstrings[i])): # for each string in that list

        hadronlistpythia = fulldecayclass.Performdecaygeneral(fulllistofstrings[i][n]) # perform the decay of that string and get back the hadron list
        #print hadronlistpythia, "pythia hadron list in testhepmc"
        eventhadronlist += hadronlistpythia
        #print eventhadronlist, "event hadron list"
    fullfinalhadrons.append(eventhadronlist)
    #print fullfinalhadrons, "full final hadrons"







total=fullfinalhadrons
#print total, "total final list"
# total.extend(Hadronlists2)
# total.extend(Hadronlists3)
Hepmcconverter.Hepmcconverter(total)
