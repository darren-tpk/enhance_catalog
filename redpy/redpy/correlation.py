# REDPy - Repeating Earthquake Detector in Python
# Copyright (C) 2016-2020  Alicia Hotovec-Ellis (ahotovec-ellis@usgs.gov)
# Licensed under GNU GPLv3 (see LICENSE.txt)

import numpy as np
import obspy.core.trace as trace
import redpy.table
import redpy.cluster
import datetime
import matplotlib
from scipy.fftpack import fft, ifft

def calcWindow(waveform, windowStart, opt, winlen=1):

    """
    Calculates the amplitude coefficient, FFT, and frequency index for a window of data.
    
    waveform: numpy array of waveform data
    windowStart: starting sample of window
    opt: Options object describing station/run parameters
    winlen: Fraction of window to use (optional)
    
    Returns windowCoeff, windowFFT, and 
    """
    
    # Shift window left by 10% of winlen
    windowStart = windowStart - opt.winlen/10
    windowCoeff = []
    windowFFT = np.zeros(opt.winlen*opt.nsta,).astype(np.complex64)
    windowFI = []
    
    for n in range(opt.nsta):
        winstart = int(n*opt.wshape + windowStart)
        winend = int(n*opt.wshape + windowStart + opt.winlen*winlen)
        fftwin = np.reshape(fft(waveform[winstart:winend]),(opt.winlen*winlen,))
        if np.median(np.abs(waveform[winstart:winend]))==0:
            windowCoeff.append(0)
            windowFI.append(np.nan) #?
        else:
            windowCoeff.append(1/np.sqrt(sum(
                waveform[winstart:winend] * waveform[winstart:winend])))
            windowFI.append(np.log10(np.mean(np.abs(np.real(
                fftwin[int(opt.fiupmin*opt.winlen*winlen/opt.samprate):int(
                opt.fiupmax*opt.winlen*winlen/opt.samprate)])))/np.mean(np.abs(np.real(
                fftwin[int(opt.filomin*opt.winlen*winlen/opt.samprate):int(
                opt.filomax*opt.winlen*winlen/opt.samprate)])))))
        
        windowFFT[n*opt.winlen:(n+1)*opt.winlen] = fftwin
        
    return windowCoeff, windowFFT, windowFI


def xcorr1x1(windowFFT1, windowFFT2, windowCoeff1, windowCoeff2, opt):

    """
    Calculates the cross-correlation coefficient and lag for two windows.
    
    windowFFT1: FFT of first window
    windowFFT2: FFT of second window
    windowCoeff1: amplitude coefficient of first window
    windowCoeff2: amplitude coefficient of second window
    Order matters for sign of lag, but not CCC.
    
    Returns maximum cross-correlation and optimal lag (in samples)
    """

    M = opt.winlen
    
    cor = []
    lag = []
    for n in range(opt.nsta):
        coeff = windowCoeff1[n] * windowCoeff2[n]
        cors = np.real(ifft(windowFFT1[n*M:(n+1)*M] * np.conj(
            windowFFT2[n*M:(n+1)*M]))) * coeff

        indx = np.argmax(cors)
        cor.append(cors[indx])

        if indx <= M/2:
            lag.append(indx)
        else:
            lag.append(indx-M)
    
    nthcor = np.sort(np.array(cor))[::-1][opt.ncor-1]
    maxcor = np.amax(cor)
    
    if nthcor >= opt.cmin:
        maxlag = np.median(np.array(lag)[np.argsort(cor)[::-1][0:opt.ncor]])
    else:
        maxlag = lag[np.argmax(cor)]

    return maxcor, maxlag, nthcor
    

def xcorr1xtable(coeffi, ffti, subtable, opt):

    """
    Correlates a new event with all events in a subtable.
    
    coeffi: amplitude coefficient of new event
    ffti: FFT of new event
    subtable: a table of either repeaters or orphans to compare to the new event
    opt: Options object describing station/run parameters
    
    Returns correlation and lag arrays
    
    The 'subtable' can be a full table (the full orphan table) or a selection of
    rows (cluster centers from repeaters, or a full family)
    
    Contemplating figuring how to run this in parallel...
    """
    
    cor = np.zeros((len(subtable),))
    lag = np.zeros((len(subtable),))
    nthcor = np.zeros((len(subtable),))

    j = -1
    for rj in subtable:
        j = j+1
        cor[j], lag[j], nthcor[j] = xcorr1x1(
            ffti, rj['windowFFT'], coeffi, rj['windowCoeff'], opt)
        
    return cor, lag, nthcor


def compare2Family(rtable, ctable, ftable, rnumber, cnum, opt):

    """
    Correlates a known repeater with all events in a family except the core.
    
    rtable: Repeater table
    ctable: Correlation matrix table
    rnumber: Row of repeater in rtable
    cnum: Cluster or family number
    opt: Options object describing station/run parameters
    
    Writes correlations to ctable
    """
    
    members = np.fromstring(ftable[cnum]['members'], dtype=int, sep=' ')
    core = ftable[cnum]['core']

    famtable = rtable[np.setdiff1d(members, core)]
    ids = famtable['id']
    rid = rtable[rnumber]['id']
        
    cor, lag, nthcor = xcorr1xtable(rtable[rnumber]['windowCoeff'],
        rtable[rnumber]['windowFFT'], famtable, opt)
        
    if len(np.where(nthcor>=opt.cmin)[0]) > 0:
        for j in range(len(cor)):
            if (nthcor[j] >= opt.cmin):
                redpy.table.appendCorrelation(ctable, rid, ids[j], cor[j], opt)


def compareDeleted(trigs, dtable, opt):

    """
    Compares trigger against deleted events
    
    trigs: Triggers to be checked
    dtable: Deleted table (manually removed from rtable)
    opt: Options object describing station/run parameters
    
    Returns trigs that do not match deleted events
    """
    
    for t in trigs:
    
        coeffi, ffti, fii = calcWindow(t.data, int(opt.ptrig*opt.samprate), opt)  
        cor, lag, nthcor = xcorr1xtable(coeffi, ffti, dtable, opt)
        
        if np.where(cor >= opt.cmin - 0.05)[0].any():
            trigs.remove(t)
    
    return trigs


def compareGoodOrphans(rtable, otable, ctable, ftable, trig, id, coeffi, ffti, cor, lag,
    nthcor, opt):

    """
    Goes and finds the matches of the new event in the orphan table, appends them to
    the repeater table, and then compares to cores
    
    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    trig: New trigger to compare
    id: Unique ID of new trigger
    coeffi: Scaling coefficient for trigger
    ffti: FFT of trigger
    cor: Correlation of trigger to orphans
    lag: Lag ot trigger to orphans
    opt: Options object describing station/run parameters
    
    """
    
    # Loop through potential matches
    written = 0
    while len(cor[cor >= opt.cmin - 0.05]) > 0:
        
        # If not written to rtable yet, realign new event
        if written == 0:
            lagmax = lag[np.argmax(cor)]
            coeffi2, ffti2, fii2 = calcWindow(trig.data, int(
                opt.ptrig*opt.samprate + lagmax), opt)
            coeffj2 = otable[np.argmax(cor)]['windowCoeff']
            fftj2 = otable[np.argmax(cor)]['windowFFT']
            fij2 = otable[np.argmax(cor)]['FI']
        # If written already, realign older orphan to new event
        else:
            coeffj2, fftj2, fij2 = calcWindow(otable[np.argmax(cor)]['waveform'],
                int(opt.ptrig*opt.samprate + lagmax - lag[np.argmax(cor)]), opt)
            
        cor2, lag2, nthcor2 = xcorr1x1(ffti2, fftj2, coeffi2, coeffj2, opt)
        
        # If actually matches...
        if nthcor2 >= opt.cmin:
            # Move both the orphans to the repeater table
            if written == 0:
                redpy.table.populateRepeater(rtable, ftable, id, trig, opt,
                    int(opt.ptrig*opt.samprate + lagmax))
                redpy.table.moveOrphan(rtable, otable, ftable, np.argmax(cor), opt)
                redpy.table.appendCorrelation(ctable, id, rtable.cols.id[-1], cor2, opt)
                written = 2
            # Update the table to reflect the new window, then move it
            else:
                otable.cols.windowFFT[np.argmax(cor)] = fftj2
                otable.cols.windowCoeff[np.argmax(cor)] = coeffj2
                otable.cols.FI[np.argmax(cor)] = fij2
                otable.cols.windowStart[np.argmax(cor)] = int(opt.ptrig*opt.samprate +
                    lagmax - lag[np.argmax(cor)])
                redpy.table.moveOrphan(rtable, otable, ftable, np.argmax(cor), opt)
                redpy.table.appendCorrelation(ctable, id, rtable.cols.id[-1], cor2, opt)
                written = written+1
                
        lag = np.delete(lag, np.argmax(cor))
        nthcor = np.delete(nthcor, np.argmax(cor))
        cor = np.delete(cor, np.argmax(cor))
    
    # If there are no actual matches in the orphans, check new event with cores
    if written == 0:
        if len(rtable) > 0:
            compareSingleOrphan2Cores(rtable, otable, ctable, ftable, trig, id, coeffi,
                ffti, opt)
        else:
            redpy.table.populateOrphan(otable, id, trig, opt)
    # If there is a match, check new event and its matches with cores
    else:
        compareMultipleOrphans2Cores(rtable, ctable, ftable, written, opt)


def compareMultipleOrphans2Cores(rtable, ctable, ftable, written, opt):

    """
    Compares multiple orphans that have already been written to the end of the repeater
    table to the other repeaters
    
    rtable: Repeater table
    ctable: Correlation matrix table
    written: Number of new repeaters written to rtable 
    opt: Options object describing station/run parameters
    
    Note: Currently only runs clustering if there are no matches to cores, and this
    is the ONLY case where full clustering is run
    """
    
    wfam = []
    wlag = []
    
    found = 0
    if len(ftable) >= 1:
        cores = rtable[ftable.cols.core[:]]
        fftjs = cores['windowFFT']
        coeffjs = cores['windowCoeff']
        ids = cores['id']
        coresNum = range(ftable.attrs.nClust) 
        
        coeffi = rtable.cols.windowCoeff[-written]
        ffti = rtable.cols.windowFFT[-written]
        cor, lag, nthcor = xcorr1xtable(coeffi, ffti, cores, opt)
        
        # Loop through families that match
        while len(cor[cor >= opt.cmin - 0.05]) > 0:    
        
            if found == 0:
                lagmax2 = lag[np.argmax(cor)]
                coeffi2, ffti2, fii2 = calcWindow(rtable[-written]['waveform'],
                    int(rtable[-written]['windowStart'] + lagmax2), opt)
                           
            cor2, lag2, nthcor2 = xcorr1x1(ffti2, fftjs[np.argmax(cor)], coeffi2,
                coeffjs[np.argmax(cor)], opt)
                
            if nthcor2 >= opt.cmin:
                if found == 0:
                    found = 1
                    wlag.append(0)
                    # Realign all new events in the repeater catalog to the matched family
                    for i in range(-written,0):
                        rtable.cols.windowCoeff[i], rtable.cols.windowFFT[i], rtable.cols.FI[i] = calcWindow(
                            rtable.cols.waveform[i], int(rtable.cols.windowStart[i] +
                            lagmax2), opt)
                        rtable.cols.windowStart[i] = int(rtable.cols.windowStart[i] + lagmax2)
                        rtable.flush()
                        ftable.cols.members[coresNum[np.argmax(cor)]] = ftable.cols.members[
                            coresNum[np.argmax(cor)]].decode('utf-8')+' {}'.format(
                            len(rtable)+i)
                        ftable.cols.printme[coresNum[np.argmax(cor)]] = 1
                        ftable.flush()
                else:
                    wlag.append(lag2)
                
                wfam.append(coresNum[np.argmax(cor)])
                
                # Compare to full family, write to correlation table
                for i in range(-written,0):
                    cor3, lag3, nthcor3 = xcorr1x1(rtable[i]['windowFFT'],
                        fftjs[np.argmax(cor)], rtable[i]['windowCoeff'],
                        coeffjs[np.argmax(cor)], opt)
                    if nthcor3 >= opt.cmin:
                        redpy.table.appendCorrelation(ctable, rtable[i]['id'],
                            ids[np.argmax(cor)], cor3, opt)
                    compare2Family(rtable, ctable, ftable, i,
                        coresNum[np.argmax(cor)], opt)
            
            else:
                members = np.fromstring(ftable[coresNum[np.argmax(cor)]]['members'],
                    dtype=int, sep=' ')
                famtable = rtable[members]
                corx, lagx, nthcorx = xcorr1xtable(coeffi2, ffti2, famtable, opt)
                if max(nthcorx) >= opt.cmin:
                    if found == 0:
                        found = 1
                        wlag.append(0)
                        # Realign all new events in the repeater catalog to the matched family
                        for i in range(-written,0):
                            rtable.cols.windowCoeff[i], rtable.cols.windowFFT[i], rtable.cols.FI[i] = calcWindow(
                                rtable.cols.waveform[i], int(rtable.cols.windowStart[i] +
                                lagmax2 + lagx[np.argmax(corx)]), opt)
                            rtable.cols.windowStart[i] = int(
                                rtable.cols.windowStart[i] + lagmax2 + lagx[np.argmax(corx)])
                            rtable.flush()
                            ftable.cols.members[coresNum[np.argmax(cor)]] = ftable.cols.members[
                                coresNum[np.argmax(cor)]].decode('utf-8')+' {}'.format(
                                len(rtable)+i)
                            ftable.cols.printme[coresNum[np.argmax(cor)]] = 1
                            ftable.flush()
                        found = 1
                    else:
                        wlag.append(lagx[np.argmax(corx)])
                        
                    wfam.append(coresNum[np.argmax(cor)])
                    
                    for x in range(len(corx)):
                        if nthcorx[x] >= opt.cmin:
                            redpy.table.appendCorrelation(ctable, rtable[-written]['id'],
                                famtable[x]['id'], corx[x], opt)
                            
            coresNum = np.delete(coresNum, np.argmax(cor))
            fftjs = np.delete(fftjs, np.argmax(cor), axis=0)
            coeffjs = np.delete(coeffjs, np.argmax(cor), axis=0)
            ids = np.delete(ids, np.argmax(cor))
            lag = np.delete(lag, np.argmax(cor))
            cor = np.delete(cor, np.argmax(cor))
        
    # Make sure to save correlation of new events with each other
    for i in range(-written+1,0):
        cor4, lag4, nthcor4 = xcorr1x1(rtable[i]['windowFFT'],
            rtable[-written]['windowFFT'], rtable[i]['windowCoeff'],
            rtable[-written]['windowCoeff'], opt)
        if nthcor4 >= opt.cmin:
            redpy.table.appendCorrelation(ctable, rtable[-written]['id'], rtable[i]['id'],
                cor4, opt)
    
    if found == 0:
        members = np.arange(len(rtable)-written,len(rtable)).astype(int)
        core = len(rtable)-written
        redpy.table.createNewFamily(rtable, ftable, members, core, opt)
    else:
        if len(wfam) == 1:
            redpy.cluster.runFamOPTICS(rtable, ctable, ftable, wfam[0], opt)
        else:
            redpy.table.mergeFamilies(rtable, ctable, ftable, wfam, wlag, opt)
        

def compareSingleOrphan2Cores(rtable, otable, ctable, ftable, trig, id, coeffi, ffti, opt):

    """
    Compares a single orphan to the cluster cores, adds the orphan to the best cluster
    if it matches, else appends to the orphan table
    
    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    trig: New trigger to compare
    id: Unique ID of new trigger
    coeffi: Scaling coefficient for trigger
    ffti: FFT of trigger
    opt: Options object describing station/run parameters
    
    """
    
    cores = rtable[ftable.cols.core[:]]
    fftjs = cores['windowFFT']
    coeffjs = cores['windowCoeff']
    ids = cores['id']
    coresNum = np.arange(ftable.attrs.nClust)
         
    cor, lag, nthcor = xcorr1xtable(coeffi, ffti, cores, opt)
    wfam = []
    wlag = []
    
    written = 0
    # Loop through potential matching families
    while len(cor[cor >= opt.cmin - 0.05]) > 0:    

        if written == 0:
            lagmax = lag[np.argmax(cor)]
            coeffi2, ffti2, fii2 = calcWindow(trig.data, int(opt.ptrig*opt.samprate + lagmax),
                opt)
        
        cor2, lag2, nthcor2 = xcorr1x1(ffti2, fftjs[np.argmax(cor)], coeffi2,
                coeffjs[np.argmax(cor)], opt)
        
        # If it definitely matches a family...
        if nthcor2 >= opt.cmin:
            if written == 0:
                # Move the orphan to the repeater table
                redpy.table.populateRepeater(rtable, ftable, id, trig, opt,
                    int(opt.ptrig*opt.samprate + lagmax))
                ftable.cols.members[coresNum[np.argmax(cor)]] = ftable.cols.members[
                    coresNum[np.argmax(cor)]].decode('utf-8')+' {}'.format(
                    len(rtable)-1)
                ftable.flush()
                wlag.append(0)
                written = 1
            else:
                wlag.append(lag2) 
            
            # Append to family list that needs to be ordered/merged
            wfam.append(coresNum[np.argmax(cor)])
            
            # Correlate with other members of the family
            redpy.table.appendCorrelation(ctable, id, ids[np.argmax(cor)], cor2, opt)
            compare2Family(rtable, ctable, ftable, -1, coresNum[np.argmax(cor)],
                opt)
            
        # Check to make sure...    
        else:
            members = np.fromstring(ftable[coresNum[np.argmax(cor)]]['members'],
                dtype=int, sep=' ')
            famtable = rtable[members]
            idf = famtable['id']
            corx, lagx, nthcorx = xcorr1xtable(coeffi2, ffti2, famtable, opt)
            if max(nthcorx) >= opt.cmin:
                if written == 0:
                    # Move the orphan to the repeater table
                    redpy.table.populateRepeater(rtable, ftable, id, trig, opt,
                        int(opt.ptrig*opt.samprate + lagmax + lagx[np.argmax(corx)]))
                    ftable.cols.members[coresNum[np.argmax(cor)]] = ftable.cols.members[
                        coresNum[np.argmax(cor)]].decode('utf-8')+' {}'.format(
                        len(rtable)-1)
                    ftable.flush()
                    wlag.append(0)
                    written = 1
                else:
                    wlag.append(lagx[np.argmax(corx)])
                wfam.append(coresNum[np.argmax(cor)])
                for x in range(len(corx)):
                    if nthcorx[x] >= opt.cmin:
                        redpy.table.appendCorrelation(ctable, id, idf[x], corx[x], opt)
                
        coresNum = np.delete(coresNum, np.argmax(cor))
        fftjs = np.delete(fftjs, np.argmax(cor), axis=0)
        coeffjs = np.delete(coeffjs, np.argmax(cor), axis=0)
        ids = np.delete(ids, np.argmax(cor))    
        lag = np.delete(lag, np.argmax(cor))
        cor = np.delete(cor, np.argmax(cor))
        
    # If doesn't match anything, append as orphan   
    if written == 0:
        redpy.table.populateOrphan(otable, id, trig, opt)
    else:
        if len(wfam) == 1:
            redpy.cluster.runFamOPTICS(rtable, ctable, ftable, wfam[0], opt)
        else:
            redpy.table.mergeFamilies(rtable, ctable, ftable, wfam, wlag, opt)


def runCorrelation(rtable, otable, ctable, ftable, ttimes, trig, id, opt):

    """
    Adds a new trigger to the correct table, runs the correlations and clustering
    
    rtable: Repeater table
    otable: Orphan table
    ctable: Correlation matrix table
    ftable: Families table
    ttimes: Trigger times
    trig: New trigger to compare
    id: Unique ID of new trigger
    opt: Options object describing station/run parameters
    
    This is the top-level logic for processing; detailed logic is within the two compare
    functions.
    """
    
    # Check to ensure this isn't a duplicate in either rtable or otable
    try:
        stime = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S.%f'))
    except ValueError:
        stime = matplotlib.dates.date2num(datetime.datetime.strptime(
            trig.stats.starttime.isoformat(), '%Y-%m-%dT%H:%M:%S'))
    
    if not len(np.intersect1d(np.where(ttimes > stime - opt.mintrig/86400), np.where(
        ttimes < stime + opt.mintrig/86400))):

        coeffi, ffti, fii = calcWindow(trig.data, int(opt.ptrig*opt.samprate), opt)
        
        # Correlate with the new event with all the orphans
        cor, lag, nthcor = xcorr1xtable(coeffi, ffti, otable, opt)
        
        try:
            # If there's a match, run the most complex function
            if max(cor) >= opt.cmin - 0.05:
                compareGoodOrphans(rtable, otable, ctable, ftable, trig, id, coeffi, ffti,
                    cor, lag, nthcor, opt)
            else:
                # Compare that orphan to the cores in the repeater table
                if len(rtable) > 0:
                    compareSingleOrphan2Cores(rtable, otable, ctable, ftable, trig, id,
                        coeffi, ffti, opt)
                # Populate as an orphan if there are no repeaters yet
                else:
                    redpy.table.populateOrphan(otable, id, trig, opt)
        except ValueError:
            print('Could not properly correlate, moving on...')
            redpy.table.populateOrphan(otable, id, trig, opt)
