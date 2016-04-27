''' Hepmc conversion module '''

import sys
import numpy
import Fourvector
import MRS
#import Twobodydecay
import Particle
import copy
import Boosts
import Decaydictionary
import Fullstringdecay
import Hadrondecay
import matplotlib.pyplot as pyplot
import hepmc 
import Gluonsplitting
import Stringformer
import Constants


INPUT = "showeroutput.hepmc" # You could also specify the absolute path of the file if it is not in the PyLund folder.
OUTPUT = "eventtest13.hepmc" # Just a file name is required in the OUTPUT field. The file is saved into the PyLund folder as the output name.

def Hepmcconverter(hadronlists):          # Takes a list of lists as the argument, and converts the products into Hepmc file format for RIVET analysis (list already weighted)
    writer = hepmc.IO_GenEvent(OUTPUT,"w")
    for i in range(len(hadronlists)):
        if i %100 == 0:
            print "EVENT ", i, "converted to hepmc"
        exec('e%i=hepmc.GenEvent()'%i)
        exec('e%i.set_event_number(i)'%i)
        exec('v%i=hepmc.GenVertex()'%i)

        eplus  = hepmc.GenParticle(hepmc.FourVector(0,0,45,45),  11) # x,y,z,E
        eminus = hepmc.GenParticle(hepmc.FourVector(0,0,-45,45),-11) # x,y,z,E
        exec('v%i.add_particle_in(eplus)'%i)
        exec('v%i.add_particle_in(eminus)'%i)




        for p in range(len(hadronlists[i])): # for each particle in each event
            fv = hadronlists[i][p].Getvector() # E,x,y,z
            ID = hadronlists[i][p].GetID()
            fvhepmc = hepmc.FourVector(fv[1],fv[2],fv[3],fv[0]) # x,y,z,E
            exec('part%i = hepmc.GenParticle(hepmc.FourVector(fv[1],fv[2],fv[3],fv[0]), int(ID))'%p)
            exec('part%i.set_status(1)'%p)
            exec('v%i.add_particle_out(part%i)'%(i,p))

        exec('e%i.add_vertex(v%i)'%(i,i))
        #exec('print e%i'%i)
        exec('e%i.set_beam_particles(eplus,eminus)'%i)
        exec('writer.write_event(e%i)'%i)




def Convertfromhepmc():
    reader = hepmc.IO_GenEvent("INPUT", "r") # ryan.hepmc is the file for PYTHIA stuff

    evtnum = 0 # should be zero
    noheavyevents = 0
    evt = hepmc.GenEvent()
    reader.fill_next_event(evt)
    totaleventlist = []
    while evt.particles_size(): # evtnum < 50000
        ''' Loop for each event. I.e. each list of particles '''
        subeventlist = []
        evtnum+=1
        #print "Event %d: " % evtnum, evt
        #print "Particles: ", evt.particles()
        #print "Num particles: ", evt.particles_size()
        #print "Num vertices: ", evt.vertices_size()
        #for p in evt.particles():
        #    m = p.momentum() # "X,Y,Z,E"
        #    print m
        fsps = evt.fsParticles()
        #print "FS particles: ", fsps
        #print "Num FS particles: ", len(fsps)
        if evtnum % 100 == 0:
            print "Event %d" % evtnum
        for p in fsps:
            ''' Build the list of final state particles in PyLund, and fill their properties '''
            #print p
            m = p.momentum() # E,x,y,z
            colour = p.flow(1)
            anticolour = p.flow(2)
            #print "1", colour, "2", anticolour

            ''' Build properties into PyLund particle class and append the particle to the sublist '''
            partfv = Fourvector.vector4D(m.e(),m.px(),m.py(),m.pz())
            ID = p.pdg_id()
            newpart = Particle.particle4D(ID,partfv)#
            #print newpart.Getinvariantmass(), "particle invariant mass"
            newpart.Setcolour(colour)
            newpart.Setanticolour(anticolour)
            subeventlist.append(newpart)
        ''' Perform some checks on the list -- current edition removes photons and vetoes heavy events '''
        if Vetoheavies(subeventlist) == True: # If event contains a heavy quark. go to next iteration before append.
            noheavyevents += 1
            #print "\n\n"
            evt.clear()
            evt = reader.get_next_event()
            continue
        subeventlist = Removephotons(subeventlist) # remove the photons from the event.
        #print "Num fs particles", len(subeventlist) # WE REMOVE THE PHOTONS. NO FS PARTICLES DOES NOT INCLUDE THE PHOTONS.

        ''' The call the colour ordering function '''

        #print subeventlist
        subeventlist = Colourorderlist(subeventlist)
        #print subeventlist, "event list after colour ordering"
        #print checkcolours(subeventlist)
        #sys.exit()
        subeventlist = masslessmomenta(subeventlist) # attempting the massels mom

        subeventlist = Fixlowenergy(subeventlist)
        #print subeventlist, "event list after fixing low nergy gluons"
        #print energies(subeventlist)
        #print checkquark(subeventlist)
        #sys.exit()
        subeventlist = Gluonsplitting.Splitgluons(subeventlist)
        #print subeventlist, "event list after splitting the gluons"
        #print len(subeventlist), "number of particles fed into string former"
        subeventlist = Stringformer.Formstrings(subeventlist) # Convert the lists to strings!
        #print subeventlist, "sublist after forming all the strings."
        #sys.exit()
        totaleventlist.append(subeventlist)


        #print "\n\n"
        evt.clear()
        evt = reader.get_next_event()


    #print noheavyevents
    return totaleventlist # returns a list of lists that contain strings ready for feeding into PyLund!

def energies(lop): # module to print energies, used for testing purposes.
    energies  = []
    for i in range(len(lop)):
        energy = lop[i].Getvector()[0]
        energies.append(energy)

    return energies

def checkquark(lop):
    energies  = []
    for i in range(len(lop)):
        energy = lop[i].Checkquark()
        energies.append(energy)

    return energies

def checkcolours(lop):
    coloursx = []
    for i in range(len(lop)):
        subcolours = []
        colours = lop[i].Getcolour()
        anticolour = lop[i].Getanticolour()
        subcolours.append(colours)
        subcolours.append(anticolour)
        coloursx.append(subcolours)
    return coloursx

def Colourorderlist(lop): #Takes a list of fs particles and orders them by quark-anticolour containing pairs.
    ''' Module that orders the list of particles in the final state by colour so that adjacent quark/antiquark pairs can be split into strings '''

    #lop = listofparticles
    no_part = len(lop)
    orderedlist = []
    no_quarks = 0
    for i in range(no_part):
        ''' find how many quark pairs we have. Need to iterate for that number of times '''
        if lop[i].Checkquark() is True:
            no_quarks +=1

    no_pairs = no_quarks / 2
    #print no_pairs, "NO PARIS"
    q_index = 0 # index of the current quark being paired up.
    sli = 0 # subloopindex. used for finding when we're done pairing up.
    firstrun = True
    for i in range(no_pairs): # for each pair of quarks
        loopcomplete = False
        #print i, "i value"

        for p in range(no_part): # NEED TO ALTER THIS LOOP FOR REPEATED RUNS
            if loopcomplete == True:
                break
            if firstrun == False:
                #print p,"in", i, "loop"
                
                p += (q_index+1)
            #print q_index, "q_index"
            #print p, "P"
            #print no_part, "no_part"

            if p > (no_part):
                #print p, "PI IN THE END CONDITION"
                #print no_part-1, "NO PART -1"
                #print len(orderedlist), "LEN orderedlist WHEN CONDITION MET"
                assert len(orderedlist) == no_part
                return orderedlist
                #break


            if lop[p].Checkquark() is True: # check for quark.
                if lop[p].Checkanti() is False: # make sure only the real quark gotten.
                    #print "QUARK FOUND, QUARK CODE CALLED"
                    q_index = p # assign the index as the quark index.
                    orderedlist.append(lop[p]) # append the quark.
                    if firstrun == False:
                        sli += 1 # need to get the index of the new quark if not the first pair.
                    firstrun = False # set the first run condition
                    

                    ''' We've found the quark -- now need to connect it via gluons to its antiquark partner '''
                    while orderedlist[sli].Getcolour() != 0: # while we do not have the antiquark at the list end. (gluons / quarks have colour).
                        current_c = orderedlist[sli].Getcolour() # current colour index value.
                        for n in range(no_part):

                            anti_c = lop[n].Getanticolour()
                            if current_c == anti_c:
                                orderedlist.append(lop[n])
                                sli += 1
                                continue # found the match, no need to continue this iteration.

                    loopcomplete = True
                                

        ''' Looping finished. Successfully found the quark--gluon--antiquark chain '''
        ''' Need to now adjust the indexing such that we dont double count lists '''

    #print len(orderedlist), "LEN OF FINAL LIST"

    return orderedlist




def Removephotons(lop):
    no_part = len(lop)
    nophotonzone = []
    for i in range(no_part):
        ID = lop[i].GetID()
        if ID != 22:
            nophotonzone.append(lop[i])

    return nophotonzone







def Vetoheavies(lop):
    ''' Function used to veto heavy containing events for current edition of PyLund '''
    ''' If Heavy is present then returns True, else False '''
    no_part = len(lop)
    heavies = [4,5,-4,-5]
    for i in range(no_part):
        ID = lop[i].GetID()
        if abs(ID) in heavies:
            return True


    return False


def Fixlowenergy(lop): # If gluon has too small an energy that it cannot be fed into PyLund after decaying, its absorbed into its neighbours.

    maxqmass = 0.6 # mass of the ud diquark pair WAS
    no_parts = len(lop)
    finished = False
    while True:
        if finished == True:
            #print "final number of particles", len(lop)
            return lop
            break
        no_parts = len(lop)
        #rint no_parts, "number particles"
        # go over the list. combine the low energy gluons absorbed into their neighbours.
        for i in range(no_parts):
            #print i, "i nubmer"
            if lop[i].Checkgluon() == True:
                vec_g = copy.deepcopy(lop[i].Getvector4D())
                vecleft = copy.deepcopy(lop[i-1].Getvector4D())
                vecright = copy.deepcopy(lop[i+1].Getvector4D())
                testleft = numpy.sqrt((vec_g + vecleft) * (vec_g + vecleft)) # invariant mass of the gluon and the partner on the left
                testright = numpy.sqrt((vec_g + vecright) * (vec_g + vecright)) # invariant mass of the gluon and the partner on the right
                Eg = vec_g.Getvector()[0]
                x = Constants.gmm() * maxqmass / Eg
                #print testleft, "testleft"
                #print testright, "testright"

                #print vec_g, "vec g into the code"
                
                if 1 < x:
                    newgluons = Absorbgluons(lop[i-1],lop[i],lop[i+1])
                    lop[i-1] = newgluons[0]
                    lop[i+1] = newgluons[1]
                    lop.pop(i) # remove the quark from the list. then need to go back to start of while to get the newest len (as it has changed)
                    #print "list popped, check if we go back to start"
                    break
                
                if testleft < (Constants.gcc()*maxqmass): # inv mass with left too small. GEt absorbed into neighbours
                    newgluons = Absorbgluons(lop[i-1],lop[i],lop[i+1])
                    lop[i-1] = newgluons[0]
                    lop[i+1] = newgluons[1]
                    lop.pop(i) # remove the quark from the list. then need to go back to start of while to get the newest len (as it has changed)
                    #print "list popped, check if we go back to start"
                    break
                if testright < (Constants.gcc()*maxqmass): # inv mass with the right too small. Get absorbed into neighbours
                    newgluons = Absorbgluons(lop[i-1],lop[i],lop[i+1])
                    lop[i-1] = newgluons[0]
                    lop[i+1] = newgluons[1]
                    lop.pop(i) # remove the quark from the list. then need to go back to start of while to get the newest len (as it has changed)
                    #print "list popped, check if we go back to start"
                    break
                else: # inv mass with neightbours is high enough. move to next gluon.
                    pass
        # get to this point once finished the lot?
            if i == (no_parts-1):
                finished = True
             

def Absorbgluons(left,current,right):

    pi = left.Getvector4D()
    veci = left.Getvector()
    pj = current.Getvector4D()
    vecj = current.Getvector()
    pk = right.Getvector4D()
    veck = right.Getvector()

    Ei = veci[0]
    Ej = vecj[0]
    Ek = veck[0]

    if Ei < Ek:
        newmoms = Combinegluons(pi,pj,pk)
        left.Setmom(newmoms[0])
        right.Setmom(newmoms[1])
        return [left,right]

    if Ei > Ek:
        newmoms = Combinegluons(pk,pj,pi)
        left.Setmom(newmoms[1])
        right.Setmom(newmoms[0])
        return [left,right]


def Combinegluons(pi,pj,pk): # takes the three fourvectors . To be used within Absorbgluons fucntion. Take pi as the smaller energy.
    #print pk, "ORIGINAL PK"
    y = (pi*pj)/(pi*pj + pi*pk + pj*pk)
    #print y, "y values"
    pkcopy1 = copy.deepcopy(pk)
    pkcopy2 = copy.deepcopy(pk)
    pkcopy1 *= (y/(1-y))
    #print pkcopy1.Getvector(), "pkcopy1"
    pkcopy2 *= (1/(1-y))
    #print pkcopy2.Getvector(), "pkcopy2"

    pinew = pi + pj - pkcopy1
    #print pi, "pi"
    #print pinew, "pinew"
    pknew = pkcopy2
    #print pk, "pk"
    #print pknew, "pknew"

    #print pinew *pinew, "square of pinew"
    #print pknew * pknew," skqure of pknew"

    return [pinew,pknew]
    

def masslessmomenta(lop):
    noparts = len(lop)
    newlist = []
    for i in range(noparts):
        if lop[i].Checkgluon() == False: # not a gluon is a quark
            ''' do momenta boosting things here '''

            mass = lop[i].Getinvariantmass()
            #print mass, "mass in invariant mass"
            vec = copy.deepcopy(lop[i].Getvector())
            #print vec, "vec  before"
            p_2 = 0
            for j in range(1,4): # just go over the momentum values, try my scaling formula

                p_2 += (vec[j]*vec[j])

            al = numpy.sqrt((mass*mass + p_2) / p_2)
            newvector = Fourvector.vector4D(vec[0], al*vec[1],al*vec[2],al*vec[3])
            ID = lop[i].GetID()
            #print ID, "id of theparticle "
            newparticle = Particle.particle4D(ID,newvector)
            newlist.append(newparticle)

            #print newparticle.Getinvariantmass(), "check inv mass"
            

        if lop[i].Checkgluon() == True:
            newlist.append(lop[i])
    return newlist
