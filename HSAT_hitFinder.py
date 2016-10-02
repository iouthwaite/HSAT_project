############                   Satellitle Analysis Program                    ############
############                   Ian Outhwaite, Oct 1 2016                      ############
############ Professor Dawn Carone, @Willaims College, now Swarthmore College ############

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
def countSat2(seq,sat,startpercent,endpercent,Name):
    l = len(seq)
    rsat = sat[::-1]
    sats = [sat,rsat]
    satlen = len(sat)
    results = []
    i = int(startpercent*l/100)
    start = i
    end = (((endpercent/100.0) * l) - satlen)
    while i< end:
        c = 0
        misses = [0,0]
        while c<satlen:
            querychar1 = sats[0][c]
            querychar2 = sats[1][c]
            seqchar = seq[i+c]
            if querychar1 != seqchar:
                misses[0] += 1
            if querychar2 != seqchar:
                misses[1] += 1
            if misses[0] > 6 and misses[1] > 6:
                break
            if c == satlen-1:
                if misses[0] < misses[1]:
                    results.append([i,misses[0]])
                else:
                    results.append([i,misses[1]])
                i = c+i
                break
            c += 1
        i+=1
        if i%1000000 == 0:
            percent = ((i-start)/(end-start))*100
            print (str(percent)[0:4] + '% done with ' + Name)
    results = sorted(results, key=lambda location: location[0]) #sort hits by location
    r = [[],[]]
    for entry in results:
        r[0].append(entry[0])
        r[1].append(entry[1])
    return r

satelliteSeqs = [
    ['hsat2A1','TTGATTCCATTAGTTTCCATTGGA'],
    ['hsat2A2','CATTCGATTCCATTCGATGATAAT'],
    ['hsat2B','TTCGATTCCATTTGATGATTCCAT']
]

def main():
    args=(sys.argv[1:])
    try:
        options = (getopt.getopt(args,"s:e:"))[0]
    except getopt.GetoptError:
        print ('Error: missing or unknown argument input')
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
                print ("Input argument after -e must be a whole number from 0->100")
                sys.exit(2)
        if arg[0] == '-s':
            temp = arg[1]
            try:
                float(temp)
                sval = int(temp)
            except ValueError:
                print ("Input argument after -s must be a whole number from 0->100")
                sys.exit(2)

    inputdir = './Genomes'
    file_new = open('HSAT_hit_results','w')
    for GenomeFile in os.listdir(inputdir):
        if (GenomeFile != '.DS_Store'):
            file_new.write("Genome Name: " + GenomeFile + "\n")
            #print "Genome Name: " + GenomeFile
            genome = getSeq("./Genomes/"+GenomeFile)
            genlen = len(genome)
            print (str(GenomeFile) + '     length: ' + str(genlen)+'bp')
            #find SATII Subfamily Locations, A1, A2, B
            for SATII in satelliteSeqs:
                seq = SATII[1]
                results = countSat2(genome,seq,sval,eval,SATII[0])
                file_new.write("%s\n" % results)

main()










'''
    misses = 0
    while c<satlen:
    querychar = s[c]
    seqchar = seq[i+c]
    if querychar != seqchar: #if the two characters you're looking at aren't the same, it is a miss
    misses += 1
    if misses > 6: #worse than 75% identity, break out
    break
    if c == satlen-1: #if you've found an entire 24mer
    count += 1
    results.append([i,misses])
    i=c+i
    break
    c += 1
    i += 1
    '''


















