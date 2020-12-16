# Calculates global synchronization index for a list of neurons
#Parameter:name=PhaseMethod;default=1
#Parameter:name=Bin (Seconds);default=0.001
#Parameter:name=Number of Shuffles;default=15
#Parameter:name=K threshold;default=3
#Option:callOnce
import nex
import numpy as np
import math
import random
import json

# The following papers describe the calculation of the global synchronization index:
    
# Li [2007]: Li X, Cui D, Jiruska P, Fox JE, Yao X, Jefferys JG. Synchronization
# measurement of multiple neuronal populations. J Neurophysiol 98: 3341–3348, 2007.

# Patel [2012]: Patel TP, Ventre SC, Meaney DF. Dynamic changes in neural circuit
# topology following mild mechanical injury in vitro. Ann Biomed Eng 40: 23–36, 2012.

# Eisenman [2015] Eisenman LN, Emnett CM, Mohan J, Zorumski CF, Mennerick S. Quantification of bursting
# and synchrony in cultured hippocampal neurons. J Neurophysiol 114: 1059–1071, 2015.

# the Python code from 
# https://github.com/lneisenman/aaft 
# is used for AAFT randomization
# the code in this script was also verified against the code in
# https://github.com/lneisenman/burst_sync/blob/master/burst_sync/global_sync/global_sync.py



def AAFTsur(xV):
    """ 
    Amplitude Adjusted Fourier Transform Surrogates
    
    https://github.com/lneisenman/aaft
    
    This code is based on the following MATLAB code

    <AAFTsur.m>, v 1.0 2010/02/11 22:09:14  Kugiumtzis & Tsimpiris
    This is part of the MATS-Toolkit http://eeganalysis.web.auth.gr/

    Copyright (C) 2010 by Dimitris Kugiumtzis and Alkiviadis Tsimpiris
                           <dkugiu@gen.auth.gr>

    Reference : D. Kugiumtzis and A. Tsimpiris, "Measures of Analysis of Time
                Series (MATS): A Matlab  Toolkit for Computation of Multiple
                Measures on Time Series Data Bases", Journal of Statistical
                Software, 2010
    """
    
    n = len(xV)
    zM = np.empty(n)
    T = np.argsort(xV)
    oxV = np.sort(xV)
    ixV = np.argsort(T)

    # Rank order a white noise time series 'wV' to match the ranks of 'xV'
    wV = np.random.randn(n) * np.std(xV, ddof=1)    # match Matlab std
    owV = np.sort(wV)
    yV = owV[ixV].copy()

    # Fourier transform, phase randomization, inverse Fourier transform
    n2 = n//2
    tmpV = np.fft.fft(yV, 2*n2)
    magnV = np.abs(tmpV)
    fiV = np.angle(tmpV)
    rfiV = np.random.rand(n2-1) * 2 * np.pi
    nfiV = np.append([0], rfiV)
    nfiV = np.append(nfiV, fiV[n2+1])
    nfiV = np.append(nfiV, -rfiV[::-1])
    tmpV = np.append(magnV[:n2+1], magnV[n2-1:0:-1])
    tmpV = tmpV * np.exp(nfiV * 1j)
    yftV = np.real(np.fft.ifft(tmpV, n))  # Transform back to time domain

    # Rank order the 'xV' to match the ranks of the phase randomized
    # time series
    T2 = np.argsort(yftV)
    iyftV = np.argsort(T2)
    zM = oxV[iyftV]  # the AAFT surrogate of xV
    return zM
    

def Histograms(doc, neurons, binSize):
    '''
    Calculates normalized rate histograms
    formula (1) in Li [2007], where Z is the rate histogram
    '''
    
    theStart = nex.GetDocStartTime(doc)
    theEnd = nex.GetDocEndTime(doc)
    theBins = np.arange(theStart, theEnd, binSize)
    hists = []
    for nr in neurons:
        ts = nr.Timestamps()
        hist, b = np.histogram(ts, bins=theBins)
        if len(ts) > 0: 
            hist = (hist - np.mean(hist))/np.std(hist)
        hists.append(hist)
    return hists
    

def PhaseArray(spikeIndexes, numBins):
    '''
    Calculates instantaneous phase arrays
    (Patel [2012] and Eisenman [2015], second formula on page 1061)
    
    uses spike bin indexes (with 1 ms time bins)
    
    non-normalized phase (numerator in formula) is 0 at spike, then 1, 2, ... , isi-1
    we use ramp array [0,1,2,...]
    
    phantom spike at time zero (we already have this spike in spikeIndexes parameter)
    
    random values between 0 and 2*pi for the times after the last spike
    '''
    
    if len(spikeIndexes) == 0:
        return 2*np.pi*np.random.random(numBins)

    phase = np.zeros(numBins)
    ramp = np.arange(0, numBins+1, 1)
    
    # isis are interspike intervals
    isis = np.diff(spikeIndexes)
    
    # start with second spike
    for i in range(1, len(spikeIndexes)):
        isi = isis[i-1]
        ind = spikeIndexes[i]
        if ind >= numBins:
            break
        if isi == 0:
            continue
        phase[ind-isi:ind] = (i-1) + ramp[0:isi]/isi
        
    # check if the last spike is less than numBins
    last = spikeIndexes[-1]
    if last < numBins:
        # adding random inst. phases
        phase[last:] = len(spikeIndexes) - 1 + np.random.random(numBins - last)
    
    # multiply by 2*pi
    return phase*2*np.pi


def PhaseVectors(doc, neurons, binSize):
    '''
    calculates instantaneous phase vectors for the given array of neurons
    see PhaseArray above
    '''
    
    theStart = nex.GetDocStartTime(doc)
    theEnd = nex.GetDocEndTime(doc)
    numBins = np.int64(np.floor((theEnd - theStart)/binSize))
    phaseVectors = []
    
    for nr in neurons:
        ts = nr.Timestamps()
        if len(ts) == 0:
            # no spikes, fill with random phases
            phaseVectors.append(2*np.pi*np.random.random(numBins))
            continue
        
        ts = np.array(ts) - theStart
        
        # check if we need to add phantom spike at zero
        if ts[0] > 0:
            ts = np.insert(ts, 0, 0)
            
        binIndexes = np.int64(np.floor(ts/binSize));
        phaseVectors.append(PhaseArray(binIndexes, numBins))
        
    return phaseVectors


def PhaseEigValues(allPhaseVectors):
    '''
    eigenvalues for phase vectors:
    Patel [2012], eigenvalues of Sjk -- phase synchronization matrix
    '''
    
    numNeurons = len(allPhaseVectors)
    S = np.zeros((numNeurons, numNeurons))
    
    for j in range(numNeurons):
        S[j][j] = 1.0
        for k in range(j+1, numNeurons):
            angleDiff = allPhaseVectors[j] - allPhaseVectors[k]
            angleDiff = np.mod(angleDiff, 2*np.pi)
            S[j, k] = np.abs(np.mean(np.exp(angleDiff*1j)))
            S[k, j] = S[j, k]
    
    # calculate eigenvalues of S       
    (eigenValues, eigenVectors) = np.linalg.eig(S)
    # sort eigenvalues
    return np.sort(np.real(eigenValues))
    

def EigValues(allHists):
    '''
    calculate eigenvalues of the equal-time correlation matrix C
    Li [2007], formula (3)
    
    we assume that allHists are already normalized as in (1) in Li [2007]
    '''
    
    numNeurons = len(allHists)
    C = np.zeros((numNeurons, numNeurons))
    for i in range(numNeurons):
        C[i][i] = 1.0
        for j in range(i+1, numNeurons):
            C[i, j] = np.mean(allHists[i]*allHists[j])
            C[j, i] = C[i, j]
            
    # calculate eigenvalues of C       
    (eigenValues, eigenVectors) = np.linalg.eig(C)
    # sort eigenvalues
    return np.sort(np.real(eigenValues))


def Reshuffle(allHists):
    '''
    reshuffle using amplitude-adjusted Fourier transform (AAFT)
    '''
    
    hsurrogates = []
    for neuronHist in allHists:
        hsurrogates.append(AAFTsur(neuronHist))
    return hsurrogates
    

def EigsOfShuffles(rateHists, numShuffles):
    '''
    Reshuffle histograms (rateHists) numShuffles times  
    and calculate eigenvalues for each shuffle
    '''
    
    numNeurons = len(rateHists)
    eigs = np.zeros((numShuffles, numNeurons))
    for i in range(numShuffles):
        print('shuffle {} of {}'.format(i+1, numShuffles))
        shuffled = Reshuffle(rateHists)
        eigs[i] = EigValues(shuffled)
    
    return eigs
    
    
def EigsOfPhaseShuffles(phaseArrays, numShuffles):
    '''
    Reshuffle instantaneous phase arrays (phaseArrays) numShuffles times  
    and calculate eigenvalues for each shuffle
    '''

    numNeurons = len(phaseArrays)
    eigs = np.zeros((numShuffles, numNeurons))
    for i in range(numShuffles):
        print('shuffle {} of {}'.format(i+1, numShuffles))
        shuffled = Reshuffle(phaseArrays)
        eigs[i] = PhaseEigValues(shuffled)
    
    return eigs


def SyncIndexPhase(doc, neurons, numShuffles, binSize, K):
    '''
    calculates global synchronization index using instantaneous phases
    the index is calculated according to Li [2007], formula (6)
    '''
    
    phaseArrays = PhaseVectors(doc, neurons, binSize)
    M = len(neurons)
    
    nonShuffledEigs = PhaseEigValues(phaseArrays)
    
    eigs = EigsOfPhaseShuffles(phaseArrays, numShuffles)
    
    means = np.mean(eigs, axis=0)
    std = np.std(eigs, axis=0)
    
    # indexes where lambdas should be zeros
    zeroLambda = ((nonShuffledEigs - (means + K*std))<=0)
    
    # normalize lambdas
    normLambda = (nonShuffledEigs - means)/(M - means)
    
    # assign zeros if lambda is less than mean_lambda + K*std
    normLambda[zeroLambda] = 0
    return normLambda


def SyncIndex(doc, neurons, numShuffles, binSize, K):
    '''
    calculates global synchronization index using normalized rate histograms
    the index is calculated according to Li [2007], formula (6)
    '''

    histograms = Histograms(doc, neurons, binSize)
    M = len(neurons)
    
    nonShuffledEigs = EigValues(histograms)
    
    eigs = EigsOfShuffles(histograms, numShuffles)
    
    means = np.mean(eigs, axis=0)
    std = np.std(eigs, axis=0)
    zeroLambda = ((nonShuffledEigs - (means + K*std))<=0)
    normLambda = (nonShuffledEigs - means)/(M - means)
    normLambda[zeroLambda] = 0
    return normLambda
    

def MakeText(graphIndex, x, y, text, fontsize=8, color='(0;0;0)'):
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


def Calculate():
    '''
    the main calculate function for Python analysis in NeuroExplorer
    - unpacks parameters
    - calls SyncIndex or SyncIndexPhase
    - generate results for NeuroExplorer
    '''
    
    doc = nex.GetActiveDocument()
    if not doc:
        raise ValueError('open data file first')
    inputJson = doc.GetPythonAnalysisInput()
    if not inputJson:
        raise ValueError('this script should be run from Python analysis only')
        
    inputPars = json.loads(inputJson)
    
    phaseMethod = int(inputPars['ScriptParameters']['PhaseMethod'])   
    binSize = float(inputPars['ScriptParameters']['Bin (Seconds)'])   
    numShuffles = int(inputPars['ScriptParameters']['Number of Shuffles'])   
    K = float(inputPars['ScriptParameters']['K threshold']) 
    
    variableNames = inputPars['Variables']
    if not variableNames:
        raise ValueError('no variables selected')
        
    theVars = []
    for v in variableNames:
        theVars.append(doc[v['Name']])

    if not theVars:
        raise ValueError('no variables selected')
    if len(theVars) < 2:
        raise ValueError('only one variable selected')

    if phaseMethod == 0:
        synIndex = SyncIndex(doc, theVars, numShuffles, binSize, K)
    else:
        synIndex = SyncIndexPhase(doc, theVars, numShuffles, binSize, K)

    globalSI = max(synIndex)
    result = {}
    result['XAxisLabel'] = 'Sync. Index Number'
    result['YAxisLabel'] = 'Sync. Index Value'
    xvalues = np.arange(1, len(synIndex)+1, 1).tolist()
    result['XValues'] = xvalues
    result['YValues'] = synIndex.tolist()
    
    graphicsCommands = []
    text = 'Global Synch. Index = {:.5f}'.format(globalSI)
    middleX = min(xvalues) + (max(xvalues)-min(xvalues))/2 - 0.5
    textY = globalSI
    if textY == 0:
        textY = 1
    graphicsCommands.append(MakeText(0, middleX, textY, text, 10))

    doc.SetPythonAnalysisOutput(json.dumps(result))
    doc.SetProperty('AdditionalGraphicsCommands', graphicsCommands)

Calculate()


