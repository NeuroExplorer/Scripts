# Finds peaks in numerical results and adds graphics commands to show peaks.
# Also adds columns to summary of numerical results with peak info.
import nex
import math
import sys
import json

def PeakDetect(x, y, delta):
#    Converted from MATLAB script at http://billauer.co.il/peakdet.html
#    Finds the local maxima in the vector y.    
#    A point is considered a maximum peak if it has the maximal
#    value, and was preceded (to the left) by a value lower by delta.

    maximums = []
    minimums = []
       
    if len(y) != len(x):
        raise ValueError('Input vectors v and x must have same length')
        
    if delta <= 0:
        raise ValueError('Input argument delta must be positive')
    
    theMinimum = float('inf')
    theMaximum = float('-inf')
    minPos = float('nan')
    maxPos = float('nan')
    
    lookformax = True
    
    for i in range(len(y)):
        this = y[i]
        if this > theMaximum:
            theMaximum = this
            maxPos = x[i]
        if this < theMinimum:
            theMinimum = this
            minPos = x[i]
        
        if lookformax:
            if this < theMaximum-delta:
                maximums.append((maxPos, theMaximum))
                theMinimum = this
                minPos = x[i]
                lookformax = False
        else:
            if this > theMinimum+delta:
                minimums.append((minPos, theMinimum))
                theMaximum = this
                maxPos = x[i]
                lookformax = True

    return maximums
    
def MakeEmptyColumn(name, numValues):
    # creates empty column object that can be used in
    # doc.SetProperty('AdditionalNumResColumns', columns)
    col = {}
    col['name'] = name
    col['values'] = []
    for row in range(numValues):
        col['values'].append('')
    return col
    
def MakeLine(graphIndex,xFrom,yFrom,xTo,yTo,thickness=0.5,color='(0;0;0)'):
    # creates line drawing command that can be used in
    # doc.SetProperty('AdditionalGraphicsCommands', graphicsCommands)
    line = {}
    line['type'] = 'line'
    line['graphNumber'] = graphIndex
    line['thickness'] = thickness
    line['color'] = color
    line['xFrom'] = xFrom
    line['xTo'] = xTo
    line['yFrom'] = yFrom
    line['yTo'] = yTo
    return line
    
def MakeText(graphIndex,x,y,text,fontsize=8,color='(0;0;0)'):
    # creates text command that can be used in
    # doc.SetProperty('AdditionalGraphicsCommands', graphicsCommands)
    textCommand = {}
    textCommand['type'] = 'text'
    textCommand['graphNumber'] = graphIndex
    # font size is in points
    textCommand['fontsize'] = fontsize
    textCommand['color'] = color
    # x and y are the coordinates of the lower-left corner of text
    textCommand['x'] = x
    textCommand['y'] = y
    textCommand['text'] = text
    return textCommand


def FindPeaksAndAddColumnsAndGraphicsCommands():
    doc = nex.GetActiveDocument()
    
    deltaPercentString = doc.GetPostProcessingScriptParameter()
    try:
        deltaPercent = float(deltaPercentString)
    except:
        deltaPercent = 25

    res = doc.GetAllNumericalResults()
    resNames = doc.GetNumResColumnNames()

    if len(resNames) < 2:
        raise ValueError('Numerical Results table need to have at least 2 columns (x an y)')
            
    # the first column should contain X axis values
    if not resNames[0] in ['Bin Left', 'Bin Middle', 'Bin Right', 'Frequency Value', 'Time']:
        raise ValueError('No x data in Numerical Results')
            
    # select only columns with real graph names 
    # (we may have columns with confidence, standard deviation, etc.)
    sumNumRes = doc.GetAllNumResSummaryData()    
    # real graph names are specified in the first column 
    # of summary of numerical results
    realGraphNames = sumNumRes[0]
    graphIndexes = []
    for i in range(len(resNames)):
        if resNames[i] in realGraphNames:
            graphIndexes.append(i)

    numResults = len(graphIndexes)        
    # first column of numerical results contains x values
    x = res[0]
    xRange = max(x) - min(x)
    # we may need to adjust x to draw peak center at bin center
    step = x[1]-x[0]
    xAdjustment = 0
    if resNames[0] == 'Bin Left':
        xAdjustment = step*0.5
    if resNames[0] == 'Bin Right':
        xAdjustment = -step*0.5
    
    # prepare empty columns with peak info
    columns = []        
    # we will report up to 5 peaks
    maxNumPeaks = 5
    for i in range(maxNumPeaks):
        columns.append(MakeEmptyColumn('PeakValue' + str(i+1),numResults))
        columns.append(MakeEmptyColumn('PeakPosition' + str(i+1),numResults))
        
    graphicsCommands = []
    # for each graph, find peaks and add columns and graphical commands
    for igraph in range(len(graphIndexes)):
        y = res[graphIndexes[igraph]]
        yRange = max(y) - min(y)
        if yRange > 0 :
            delta = yRange*deltaPercent/100.0
            peaks = PeakDetect(x, y, delta)
            for ip in range(len(peaks)):
                xp = peaks[ip][0] + xAdjustment
                yp = peaks[ip][1]
                if ip < maxNumPeaks:
                    # fill extra columns
                    columns[ip*2]['values'][igraph] = '%.6f' % (yp)
                    columns[ip*2+1]['values'][igraph] = '%.6f' % (xp)
                # add graphics commands
                dx = xRange*0.010
                dy = yRange*0.015
                # add graphics commands to draw 'x' (red diagonal lines) at each peak
                graphicsCommands.append(MakeLine(igraph,xp-dx,yp-dy,xp+dx,yp+dy,2.0,'(255;0;0)'))
                graphicsCommands.append(MakeLine(igraph,xp-dx,yp+dy,xp+dx,yp-dy,2.0,'(255;0;0)'))
                # add text command showing x value of the peak
                # if you want both x and y, use
                # peakText = '(%.3f,%.3f)'%(xp,yp)
                peakText = '%.3f'%(xp) 
                graphicsCommands.append(MakeText(igraph,xp,yp+dy,peakText+ ' s',7,'(0;0;255)'))

        
    doc.SetProperty('AdditionalNumResColumns', columns)
    doc.SetProperty('AdditionalGraphicsCommands', graphicsCommands)


FindPeaksAndAddColumnsAndGraphicsCommands()
