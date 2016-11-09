#Calculates the spectrum of the selected signal range via DFTs (no zero-padding).
#    Range selection is done using Select Data, From, To in Data Selection tab.
#    Assumes that less than 100,000 data points are selected.
#    Returns spectrum values as % of total spectrum.
#Parameter:name=Show Frequency From;default=0
#Parameter:name=Show Frequency To;default=100000
import nex
import json
import math
import cProfile

# helper function to select data from specified time range
def IntervalFilterVarIfNeeded( doc, var, pars ):
    selData = pars['Select Data']
    theFilter = pars['Interval Filter']
    if selData == 'All' and theFilter == 'None':
        # no need to filter
        return var
        
    mainFrom = 0.0
    mainTo = nex.GetDocEndTime(doc)
    if selData == 'From Time Range':
        mainFrom = float(pars['Select Data From (sec)'])
        mainTo = float(pars['Select Data To (sec)'])
        
    tempInt = nex.NewIntEvent(doc, 0)
    nex.AddInterval(tempInt, mainFrom, mainTo)
    
    if theFilter == 'None':
        return nex.IntervalFilter(var, tempInt)
    else:
        tempInt1 = nex.IntAnd(tempInt, doc[theFilter])
        return nex.IntervalFilter(var, tempInt1)        
        
# main calculate function of the script 
# called in the last line of the script below        
def Calculate():
    doc = nex.GetActiveDocument()
    # get the variable info and values of script parameters
    inputJson = doc.GetPythonAnalysisInput()
    inputPars = json.loads(inputJson)
    
    # get all analysis parameters (we need to get data selection options)
    parsJson = doc.GetAllAnalysisParameters()
    pars = json.loads(parsJson)
    
    variable = inputPars['Variable']

    if variable['Type'] != 'Continuous':
        raise ValueError('This analysis can only analyze Continuous variables')
        
    filteredVar = IntervalFilterVarIfNeeded(doc, doc[variable['Name']], pars) 

    v = filteredVar.ContinuousValues()
    numberOfContValues = len(v)
    if numberOfContValues > 100000:
        raise ValueError('Too many Continuous values')
   
    freqFrom = float(inputPars['ScriptParameters']['Show Frequency From'])
    freqTo = float(inputPars['ScriptParameters']['Show Frequency To'])
    samplingRate = float(variable['SamplingRate'])
            
    spectrumStep = samplingRate/numberOfContValues
    spectrum = nex.Spectrum(v)
    s = sum(spectrum)
        
    result = {}
    result['XAxisLabel'] = 'Frequency (Hz)'
    result['YAxisLabel'] = 'Spectrum (%)'
    result['XValues'] = []
    result['YValues'] = []
    for i in range(len(spectrum)):
        freq = i*spectrumStep
        if freq >= freqFrom and freq <= freqTo:
            result['XValues'].append(freq)
            result['YValues'].append(100.0*spectrum[i]/s)
            
    if len(result['XValues']) == 0:
        raise ValueError('Frequency range is too narrow: frequency step = ' + str(spectrumStep))
        
    doc.SetPythonAnalysisOutput(json.dumps(result))
    
cProfile.run('Calculate()')

#Calculate()
