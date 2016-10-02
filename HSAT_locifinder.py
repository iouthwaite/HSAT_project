#HSATII Loci Finder Program
#By Ian Outhwaite, 10/1/2016, iro1@williams.edu, iouthwaite1@gmail.com,
#For: Professor Dawn Carone, Swarthmore College

#REQUIRES SCIPY
#command line: first the reference dataset, and then the results you wish to analyze, and then any "-" commands
#command line input format: python HSAT_locusfinder.py NameOfCuratedLociFile.txt NameOfHSATHitsResults.txt (optional commands -h _ -u _ -c)
# -c means search for combined data (loci = of all three types, HSATIIA1, A2, B)
# -h means minimum hits per loci
# -u means stringency of loci - I found 8 worked well, but you can move it up or down as you wish
#DEFAULTS: minhits per loci is 25, and the kernel (named uk, or unlikelynumber sometimes) is set to 3, so looking at a location +/- 3 other locations around it

#Formatted for python, will run on 2.7 and 3.5

import os, subprocess
import sys, getopt
import math
from random import randint
import random
from scipy.stats import beta

#imports a sequence from a text file
def getSeq(filename):
    file_in = open(filename, "r")
    data = ''
    for line in file_in:
        data = data+line
    return data

#gerneates a list of diatances between hits
def generateDistances(locations):
    dists = []
    for i in range(0,len(locations)-1):
        if locations[i+1]>locations[i]:
            dists.append(locations[i+1]-locations[i])
    return dists

#computes the mean from a list
def computeMeanDist(data):
    temp = 0.0
    for i in range(0,len(data)):
        temp = temp+data[i]
    temp = temp/len(data)
    return temp

#gets the HSATII example distances from a list of lists text file of the locations
def getListsFromLocations(textfile):
    tempstuff = getSeq(textfile)
    temp = []
    inner = []
    temp2 = ''
    pchar = ''
    for charr in tempstuff:
        if ((charr == ',') and (temp2 != '') and (pchar != ']')):
            inner.append(int(temp2))
            temp2 = ''
        if ((charr == ']') and (temp2 != '')):
            inner.append(int(temp2))
            temp.append(inner)
            inner = []
            temp2 = ''
        if ((charr != '[') and (charr != ' ') and (charr != ',') and (charr != ']')):
            temp2 = temp2 + charr
        pchar = charr
    dists = []
    for entry in temp:
        temp = generateDistances(entry)
        for entry2 in temp:
            dists.append(entry2)
    return dists

#gets the HSATII example distances from a text file of the distances (format, list of lists)
def getListsFromDistances(textfile):
    tempstuff = getSeq(textfile)
    temp = []
    inner = []
    temp2 = ''
    pchar = ''
    for charr in tempstuff:
        if ((charr == ',') and (temp2 != '') and (pchar != ']')):
            inner.append(int(temp2))
            temp2 = ''
        if ((charr == ']') and (temp2 != '')):
            inner.append(int(temp2))
            temp.append(inner)
            inner = []
            temp2 = ''
        if ((charr != '[') and (charr != ' ') and (charr != ',') and (charr != ']')):
            temp2 = temp2 + charr
        pchar = charr
    data = []
    for entry in temp:
        for item in entry:
            data.append(item)
    return data

#computes the standard deviation from a list of values and a mean
def computeSD(data,m2):
    sd = 0.0
    for i in range(0,len(data)):
        sd = sd + ( (data[i] - m2)**2 )
    sd = (sd/(len(data)-1))**.5
    return sd

#creates a bootstrap dataset 100x larger than the input dataset
def generateBootstrapData(distances):
    dists = []
    for i in range (0,100):
        bdata = []
        for j in range(0,len(distances)):
            randintg = randint(0,len(distances)-1)
            temp = distances[randintg]
            bdata.append(temp)
        dists = dists + bdata
    return dists

#gets the entries under a certain value
def getDUnder(data,val):
    ndata = []
    for entry in data:
        if entry<val:
            ndata.append(entry)
    return ndata

#gives the entries over a certain value
def getDOver(data,val):
    ndata = []
    for entry in data:
        if entry>val:
            ndata.append(entry)
    return ndata

#returns the largest value out of a list
def findmaxvalue(data):
    temp = 0
    for entry in data:
        if entry>temp:
            temp = entry
    return temp

#determines what values to use for the beta distribution, as well as giving other stats about the distribution.
def determineDistribution2(bootstrapdatadists):
    
    overalm = computeMeanDist(bootstrapdatadists)
    overalsd = computeSD(bootstrapdatadists,overalm)
    
    modelbp = int(overalm+(3*overalsd))
    
    over = getDOver(bootstrapdatadists,modelbp)
    under = getDUnder(bootstrapdatadists,modelbp)
    
    punif = (len(over) + 0.0)/len(bootstrapdatadists)
    max = findmaxvalue(bootstrapdatadists)
    
    meanunder = computeMeanDist(under)
    sdUnder = computeSD(under,meanunder)
    
    u = overalm/(max+0.0)
    sd = overalsd/(max+0.0)
    Vr = sd*sd
    alpha = ( ((1-u)/Vr) - (1/u)) * (u*u)
    beta = alpha*((1/u)-1)
    
    print ('')
    print ("Model Information")
    print('')
    print ("Upper limit for model: " +str(max) + "bp")
    print ("Overall mean distance between hits: " + str(int(overalm))+ "bp")
    print ("Overall standard deviation for distances between hits: " +str(overalsd)+ "bp")
    print ("alpha: " + str(alpha))
    print ("beta: " + str(beta))
    print ('')
    print ("Three standard deviations from the mean: " + str(int(overalm+(3*overalsd))) + "bp")
    print ("Mean distance for distances under " + str(modelbp) + "bp: " + str(int(meanunder))+ "bp")
    print ("Standard Deviation for distances under " + str(modelbp) + "bp: " + str(sdUnder)+ "bp")
    print ("Probability a distance is more than " + str(modelbp) + "bp: " + str(punif))
    return [punif,alpha,beta,max,modelbp]

#finds loci of hits
#locations are the original data of locations
#alpha and beta1 refer to the beta distribution
#minhits are the minimum number of hits per loci in order to count one
#uk is the kernel size; if a locations =/- uk locations around it is considered too unlikely (<0.05, or <5% of what is expected in the model) then it is excluded
#maxval is the largest value considered "possible", distances >= maxval are set equal to maxval-1 since in the beta distribution the probability of 1 is 0, and we'd like to give them a chance, even if that chance is essentially 0.
def LociFinder(locations,mbp,alpha,beta1,minhits,uk,maxval,pval):
    loci = []
    dists = generateDistances(locations)
    mean = beta.mean(alpha,beta1)
    std = beta.std(alpha,beta1)
    limit = 0.05
    dists = generateDistances(locations)
    locs = [] #list of the indices of potentially good distances
    i = 0
    while i<len(dists): #denote "good" locations that might be ok given their probability of occuring
        start = i-uk
        end = i+uk
        if start < 0:
            start = 0
        if end > (len(dists)-1):
            end = len(dists)-1
        dlist = dists[start:end+1]  #list of the distances +/- "uk" from the seed, i
        plist = []
        for item in dlist:
            if (item >= mbp) and (item <= maxval):
                plist.append(pval)
            else:
                plist.append(beta.pdf( (item/(maxval + 0.0)),alpha,beta1))

        if len(plist) >= 1:
            cumulativep = sum(plist) / (float(len(plist)))
        else:
            cumulativep = 0
        if cumulativep > limit:
            locs.append(i)
        i+=1
    i = 0
    currentlocs = []
    while i<len(locs): #groups the "good" locations according to whether or not they have enough adjacent neighbors that are also "good"
        t = 0
        if i != 0:
            if locs[i-1] == locs[i] -1:
                t = 1
                currentlocs.append(locs[i])
                if (i != (len(locs)-1)):
                    if (locs[i+1] != locs[i]+1):
                        if len(currentlocs) >= minhits:
                            startloc = locations[currentlocs[0]]
                            endloc = locations[currentlocs[len(currentlocs)-1]]
                            start = currentlocs[0]
                            end = currentlocs[len(currentlocs)-1]
                            loci.append([[startloc,endloc],len(currentlocs),start,end])
                        currentlocs=[]
        if (i != (len(locs)-1)) and (t==0):
            if locs[i+1] == locs[i]+1:
                currentlocs.append(locs[i])
                t = 1
        if t != 1:
            if len(currentlocs) >= minhits:
                startloc = locations[currentlocs[0]]
                endloc = locations[currentlocs[len(currentlocs)-1]]
                start = currentlocs[0]
                end = currentlocs[len(currentlocs)-1]
                loci.append([[startloc,endloc],len(currentlocs),start,end])
            currentlocs=[]
        if i == len(locs)-1:
            if len(currentlocs) >= minhits:
                startloc = locations[currentlocs[0]]
                endloc = locations[currentlocs[len(currentlocs)-1]]
                start = currentlocs[0]
                end = currentlocs[len(currentlocs)-1]
                loci.append([[startloc,endloc],len(currentlocs),start,end])
        i+=1
        t = 0
    return loci

#method for looking at just 1 HSATII subfamilly at a time
def getAllHSATIIData(textfile,referencefile,minhits,uk):
    data = getSeq(textfile)
    referencedata = getListsFromLocations(referencefile)
    print ('')
    print ('Generating Bootstrap Data...')
    bdata = generateBootstrapData(referencedata)
    print ('')
    print ('Generating Model Parameters...')
    backgrounddata = determineDistribution2(bdata)
    modelbreakpoint = backgrounddata[4]
    probval = backgrounddata[0]
    maxval = backgrounddata[3]
    alpha = backgrounddata[1]
    beta = backgrounddata[2]
    Chromosome = ''
    c = []
    currentCHR = 'null'
    data2 = data.split('\r')
    if len(data2) == 1:
        data2 = data2[0].split('\n')
    for entry in data2:
        name = [currentCHR,'null']
        if len(entry) > 0:
            if entry[0] == 'G':
                temp = entry.split(' ')
                Chromosome = temp[2]
                hsatcounter = 0
                print ('')
                print (Chromosome)
                currentCHR = Chromosome
                name = [currentCHR,'null']
            else:
                nameNumber = hsatcounter%3
                if nameNumber == 0:
                    print ('HSATIIA1')
                    name[1] = '_HSATIIA1'
                if nameNumber == 1:
                    print ('HSATIIA2')
                    name[1] = '_HSATIIA2'
                if nameNumber == 2:
                    print ('HSATIIB')
                    name[1] = '_HSATIIB'
                hsatcounter += 1
                data3 = entry[2:len(entry)-2].split('], [')
                templocs = data3[0].split(', ')
                temphits = data3[1].split(', ')
                locations = []
                j = 0
                while j<len(temphits):
                    locations.append([int(templocs[j]),int(temphits[j])])
                    j += 1
                if len(locations) > 0:
                    c.append( [getLocations(locations,modelbreakpoint,alpha,beta,minhits,uk,maxval,probval,name[1]),(name[0]+name[1])] )
                else:
                    c.append([[[]],name[0]+name[1]])
    return c

#For getting data from all HSAT loci together
def getCombinedLociHSATIIData(textfile,referencefile,minhits,uk):
    data = getSeq(textfile)
    referencedata = getListsFromLocations(referencefile)
    print ('')
    print ('Generating Bootstrap Data...')
    bdata = generateBootstrapData(referencedata)
    print ('')
    print ('Generating Model Parameters...')
    backgrounddata = determineDistribution2(bdata)
    modelbreakpoint = backgrounddata[4]
    probval = backgrounddata[0]
    maxval = backgrounddata[3]
    alpha = backgrounddata[1]
    beta = backgrounddata[2]
    data2 = data.split('\r')
    if len(data2) == 1:
        data2 = data2[0].split('\n')
    Chromosome = ''
    locations = []
    c = []
    if len(data2) > 0:
        for entry in data2:
            if len(entry) > 0:
                if entry[0] == 'G':
                    temp = entry.split(' ')
                    oldc = Chromosome
                    Chromosome = temp[2]
                    c.append([getLocations(locations,modelbreakpoint,alpha,beta,minhits,uk,maxval,probval,'all'),oldc])
                    locations = []
                    print ('')
                    print (Chromosome)
                else:
                    data3 = entry[2:len(entry)-2].split('], [')
                    templocs = data3[0].split(', ')
                    temphits = data3[1].split(', ')
                    j = 0
                    while j<len(temphits):
                        locations.append([int(templocs[j]),int(temphits[j])])
                        j += 1
        temp = [getLocations(locations,modelbreakpoint,alpha,beta,minhits,uk,maxval,probval,'all'),Chromosome]
        c.append(temp)
    else:
        c.append([[[]],Chromosome])
    return c
        
def getLocations(locations,modelbreakpoint,alpha,beta,minhits,uk,maxval,probval,subfamily):
    c = []
    l = sorted(locations,key=lambda tup: tup[0])
    locations = []
    hittypes = []
    for entry in l:
        locations.append(entry[0])
        hittypes.append(entry[1])
    lociLocs = LociFinder(locations,modelbreakpoint,alpha,beta,minhits,uk,maxval,probval)
    j = 0
    while j<len(lociLocs):
        loc = lociLocs[j]
        start = loc[2]
        end = loc[3]
        mval = 0.0
        for i in range(start,end):
            mval = mval + hittypes[i]
        mval = mval/(end-start)
        temp = []
        print ('Loci location: ' + str(loc[0]))
        temp.append(loc[0][0])
        temp.append(loc[0][1])
        print ('Number of Hits in Loci: ' + str(loc[1]))
        temp.append(loc[1])
        print ('Mean Value of Hits in Loci: ' + str(int((((24.0-mval)/24)*100))) + "% identity")
        temp.append(str(int((((24.0-mval)/24)*100))))
        temp.append(subfamily)
        c.append(temp)
        j+= 1
    return c

def writeToSheet(c):
    import xlwt
    book = xlwt.Workbook(encoding="utf-8")
    sheet1 = book.add_sheet("Sheet 1")
    sheet1.write(0, 0, "Name")
    sheet1.write(0, 1, "start value")
    sheet1.write(0, 2, "end value")
    sheet1.write(0, 3, "number of hits")
    sheet1.write(0, 4, "percent identity of hits in loci")
    sheet1.write(0, 5, "HSATII subfamily")
    x = 1
    for entry1 in c:
        name = entry1[1]
        t = 0
        for entry in entry1[0]:
            if len (entry) > 0:
                sheet1.write(x, 0, name)
                sheet1.write(x, 1, str(entry[0]))
                sheet1.write(x, 2, str(entry[1]))
                sheet1.write(x, 3, str(entry[2]))
                sheet1.write(x, 4, str(entry[3]))
                sheet1.write(x, 5, str(entry[4]))
                x += 1
                t = 1
        if t==0:
            sheet1.write(x, 0, name)
            x += 1
    book.save("HSATfileforplotproduction.xls")

def main():
    args=(sys.argv[1:])
    opts = args[2:]
    combined = 0
    c=[]
    try:
        options = (getopt.getopt(opts,"ch:u:"))[0]
    except getopt.GetoptError:
        print ('Error: missing or unknown argument input')
        sys.exit(2)
    referencehits = str(args[0])
    name = str(args[1])
    minhits = 25
    unlikelynumber = 2
    for arg in options:
        if arg[0] == '-h':
            temp = arg[1]
            try:
                float(temp)
                temp2 = int(temp)
                if temp2 > 1:
                    minhits = temp2
            except ValueError:
                print ("Input argument after -h, minimum hits per loci, must be a whole number greater than 1")
                sys.exit(2)
        if arg[0] == '-u':
            temp = arg[1]
            try:
                float(temp)
                temp2 = int(temp)
                if temp2 > 2:
                    unlikelynumber = temp2
            except ValueError:
                print ("Input argument after -u must be a whole number greater than 2")
                sys.exit(2)
        if arg[0] == '-c':
            combined = 1
    if combined == 0:
        c = getAllHSATIIData(name,referencehits,minhits,unlikelynumber)
    else:
        c = getCombinedLociHSATIIData(name,referencehits,minhits,unlikelynumber)
    print ('')
    writeToSheet(c)



main()




