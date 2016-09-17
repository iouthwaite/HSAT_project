############     Satellitle Analysis Program         ############
############         Ian Outhwaite, Dec 23 2015      ############
############ Professor Dawn Carone, Willaims College ############

import sys, getopt
import math
from random import randint
import os, subprocess

#imports a sequence from a text file
def getSeq(filename):
    file_in = open(filename, "r")
    data = ''
    for line in file_in:
        data = data+line
    return data

#Counts hits of HSAT with a minimum of 75% fidelity to original 24mers. Does not include indels. Counts hits in both directions, and returns the locations of hits as well as the number of missed bases in each hit, up to a maximum of 6 mismatches per 24mer (aka 75%).
#Same as above, except if the first character isn't correct it moves on to the next cycle instead of making up to six comparisons before determining that the hit is non-existent. There is a possibility for missing some hits, but it is roughly 6x faster. Better for getting a general idea of where loci might be.
def countSat2(seq,sat,startpercent,endpercent):
    l = len(seq)
    rsat = sat[::-1]
    sats = [sat,rsat]
    seqlen = len(seq)
    satlen = len(sat)
    count = 0
    results = [[],[]] #locations, number of misses per location
    for s in sats:
        i = int(startpercent*l/100)
        while i<(seqlen-satlen+1-((100-endpercent)*l/100)):
            c = 0
            misses = 0
            while c<satlen:
                querychar = s[c]
                seqchar = seq[i+c]
                if querychar != seqchar:
                    misses += 1
                if misses > 6:
                    break
                if c == satlen-1:
                    count += 1
                    results[0].append(i)
                    results[1].append(misses)
                    i=c+i
                    break
                c += 1
            i += 1
    results.sort()
    return results

hsat2A1 = 'TTGATTCCATTAGTTTCCATTGGA'
hsat2A2 = 'CATTCGATTCCATTCGATGATAAT'
hsat2B = 'TTCGATTCCATTTGATGATTCCAT'

satelliteSeqs = [hsat2A1,hsat2A2,hsat2B]

def main():
    args=(sys.argv[1:])
    try:
        options = (getopt.getopt(args,"s:e:"))[0]
    except getopt.GetoptError:
        print 'Error: missing or unknown argument input'
        sys.exit(2)
    sval = 0
    eval = 100
    for arg in options:
        if arg[0] == '-e':
            temp = arg[1]
            try:
                float(temp)
                eval = int(temp)
            except ValueError:
                print "Input argument after -e must be a whole number from 0->100"
                sys.exit(2)
        if arg[0] == '-s':
            temp = arg[1]
            try:
                float(temp)
                sval = int(temp)
            except ValueError:
                print "Input argument after -s must be a whole number from 0->100"
                sys.exit(2)

    inputdir = './Genomes'
    file_new = open('results','w')
    for GenomeFile in os.listdir(inputdir):
        if (GenomeFile != '.DS_Store'):
            file_new.write("Genome Name: " + GenomeFile + "\n")
            #print "Genome Name: " + GenomeFile
            genome = getSeq("./Genomes/"+GenomeFile)
            genlen = len(genome)
            print GenomeFile, genlen
            '''
            #find SATII Subfamily Locations, A1, A2, B
            for SATII in satelliteSeqs:
                results = countSat2(genome,SATII,sval,eval)
                print results
                i = 0
                temp = []
                for item in results:
                    if i==0:
                        temp = temp+item
                        i = 1
                    else:
                        temp = temp+item
                        file_new.write("%s\n" % temp)
                        temp = []
            '''

main()



























