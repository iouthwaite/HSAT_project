import datetime
import numpy as np
import plotly.plotly as py
py.sign_in('', '')
import plotly.graph_objs as go
import xlrd

workbook = xlrd.open_workbook('HSATfileforplotproduction_c2.xls')
sheet = workbook.sheet_by_name('Sheet 1')


def importGenomeLengths():
    wb = xlrd.open_workbook('Genome Lengths.xlsx')
    sheet2 = wb.sheet_by_name('Lengths')
    dict = {}
    x = 0
    while x<sheet2.nrows:
        n = str(sheet2.cell(x,0).value)
        val = int(sheet2.cell(x,1).value)
        dict[n] = val
        x += 1
    return dict

lengths = importGenomeLengths()

names = []

def getData():
    data = []
    cname = ''
    x = 1
    j = -1
    lim = sheet.nrows
    while x<sheet.nrows:
        n = str(sheet.cell(x,0).value)
        if n != cname:
            cname = n
            data.append([])
            j += 1
        start = (str(sheet.cell(x,1).value))
        dat = []
        if start != '':
            start = int(start)
            end = int((str(sheet.cell(x,2).value)))
            for item in lengths:
                if item in cname:
                    k = len(item)
                    if k == len(cname):
                        tlen = lengths[item]
                        start = int(((int(start)+0.0)/tlen) * 100)
                        end = int(((int((str(sheet.cell(x,2).value)))+0.0)/tlen) * 100)
                        #if len(cname) > k:
                    elif cname[k] == '_':
                        tlen = lengths[item]
                        start = int(((int(start)+0.0)/tlen) * 100)
                        end = int(((int((str(sheet.cell(x,2).value)))+0.0)/tlen) * 100)
            #print str(start) + '  ' + str(end)
            hits = int((str(sheet.cell(x,3).value)))
            dat=dat + [start,end,hits]
        names.append(cname)

        l = len(data)
        if l>=1:
            data[l-1].append(dat)
        x+=1
    return data

data = getData()

numchrs = len(data)
programmers = range(1,numchrs+1)
      
z = []
x = 0
while x<numchrs:
#for chr in data:
    new_row = []
    t = 0
    while t < 101:
        new_row.append(0)
        t+=1
    z.append(list(new_row))
    x += 1
  
counter = 0
x = 0
while x < numchrs:
#for chromosome in data:
    chromosome = data[x]
    if len(chromosome[0]) > 0:
        for entry in chromosome:
            start = entry[0]
            end = entry[1]
            zval = entry[2]
            while start <= end:
                print start
                z[counter][start] = zval
                start += 1
    counter += 1
    x += 1

xcord = [x for x in range(0,101)]
data = [
    go.Heatmap(
        z=z,
        x=xcord,
        y=names,
        colorscale='Viridis',
    )
]

layout = go.Layout(
    autosize =False,
    height = 5000,
    width = 5000,
    margin=go.Margin(l=200),
    title='Number of HSATII hits inside HSATII Loci Accross Primate Chromosomes',
    xaxis = dict(title='Percent Distance Along Chromosome',ticks='', nticks=20),
    yaxis = dict(title='Chromosome',ticks='', nticks=numchrs*2 )
)

fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename='datetime-heatmap')
