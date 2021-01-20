#!/usr/bin/python

import numpy as np
from collections import defaultdict
import inspect
from lib import generalUtil as gUtil


def conFunc(x, a):
    return a


def fVar(flux, fluxError):
    # Taken from Eq.1 of the Mrk501 paper draft
    fVar = np.std(flux)**2
    fVar -= np.average(fluxError**2)
    fVar /= np.average(flux)**2
    if fVar < 0:
        fVar = 0
    return np.sqrt(fVar)


def deltaFVar(flux, fluxError):
    # Taken from Eq.2 of the Mrk501 paper draft
    fVarNow = fVar(flux, fluxError)
    sigErr = np.average(fluxError**2)
    mean = np.average(flux)
    N = len(flux)

    if mean > 0 and N > 0:
        sigNXS = (np.sqrt(2/N)*(sigErr/(mean**2)))**2
        sigNXS += (np.sqrt(sigErr/N)*(2*fVarNow/(mean)))**2
        sigNXS = np.sqrt(sigNXS)
    else:
        sigNXS = 0

    return (np.sqrt(fVarNow**2 + sigNXS) - fVarNow)


def calcBlockFlux(time, values, uncert, edges):
    """
       Calculate mean flux and width for each block.
       Each data point is weighted by its uncertainty
    """

    fluxMean, fluxErrorMean = list(), list()

    for i in range(len(edges)-1):
        maskNow = (time >= edges[i]) & (time <= edges[i+1])
        fluxMean = np.append(fluxMean,
                             np.average(values[maskNow], weights=1/uncert[maskNow]))
        uncertSumQuad = np.sum(np.power(uncert[maskNow], 2))
        fluxErrorMean = np.append(fluxErrorMean, np.sqrt(uncertSumQuad/len(uncert[maskNow])))

    # Last block
    maskNow = time >= edges[-1]
    fluxMean = np.append(fluxMean,
                         np.average(values[maskNow], weights=1/uncert[maskNow]))
    uncertSumQuad = np.sum(np.power(uncert[maskNow], 2))
    fluxErrorMean = np.append(fluxErrorMean, np.sqrt(uncertSumQuad/len(uncert[maskNow])))

    return(fluxMean, fluxErrorMean)


def readGammaFromFile(logger, directory, priorFile, sourceNow, prob):

    gammaPrefix = '{}/determinePrior/{}/'.format(directory, sourceNow)
    fileName = gammaPrefix + priorFile

    if gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3], False):

        headersType = {'names': ('ncp', 'fp'),
                       'formats': ('f8', 'f8')}
        gammaData = np.loadtxt(fileName, dtype=headersType)

        ncp = gammaData[(np.abs(gammaData['fp'] - prob)).argmin()]['ncp']
        gamma = 10**(-ncp)

        return gamma

    else:
        gamma = 0.0002
        logger.warn([['r', 'Could not find '],
                     ['r', fileName],
                     [' using gamma -'],
                     ['b', gamma]])

        return gamma


def readVeritasLC(logger, fileName):
    """
    Read a VERITAS lightcurve from a text file.
    The structure of the file is given in headersType below.
    The fluxes are returned, those are where the flux error
    is different to zero.
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    headersType = {'names': ('DateMJD', 'Date Error',
                             'Flux', 'Flux Error'),
                   'formats': ('f8', 'f8',
                               'f8', 'f8')}

    veritasData = np.loadtxt(fileName, dtype=headersType)
    veritasFluxes = veritasData[veritasData['Flux Error'] > 0]
    veritasUpperLimits = veritasData[veritasData['Flux Error'] == 0]

    return veritasFluxes


def readVeritasObservations(logger, fileName):

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    headersType = {'names': ('run', 'date', 'flux', 'fluxError',
                             'significance', 'ze'),
                   'formats': ('f8', 'f8', 'f8', 'f8',
                               'f8', 'f8')}
    veritasObs = np.loadtxt(fileName, dtype=headersType)

    return veritasObs


def readFermiLC(logger, fileName):
    """
    Read the  Fermi lightcurve from an npy file saved by FermiPy.
    Convert the information to a structured numpy array to make it easier
    to use later.
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    fermiLC = np.load(fileName, encoding='latin1', allow_pickle=True).flat[0]
    fermiLC = convertFermiDictToNumpyArray(fermiLC)

    return fermiLC


def readSwiftLC(logger, fileName):
    """
    Read swift lightcurve obtained from the online tool.
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    logger.info([['wb', 'Reading Swift LC from'],
                 ['y', fileName]])

    swiftFile = open(fileName, 'r')
    date, dateErrUp, dateErrDn = list(), list(), list()
    rate, rateErrUp, rateErrDn = list(), list(), list()
    rateErrAve, mode = list(), list()
    for line in swiftFile:
        if '!' in line:
            if 'WT data' in line:
                modeNow = 'WT'
                continue
            if 'PC data' in line:
                modeNow = 'PC'
                continue
            if 'Upper limit' in line:
                break
        if '!' not in line and len(line) > 1 and 'NO' not in line and 'READ' not in line:
            date.append(float(line.split()[0].strip()))
            dateErrUp.append(abs(float(line.split()[1].strip())))
            dateErrDn.append(abs(float(line.split()[2].strip())))
            rate.append(float(line.split()[3].strip()))
            rateErrUp.append(abs(float(line.split()[4].strip())))
            rateErrDn.append(abs(float(line.split()[5].strip())))
            rateErrAve.append((rateErrUp[-1] + rateErrDn[-1])/2.)
            mode.append(modeNow)

    swiftData = np.c_[date, dateErrDn, dateErrUp,
                      rate, rateErrDn, rateErrUp,
                      rateErrAve, mode]

    headersType = {'names': ('Date', 'Date error down', 'Date error up',
                             'Rate', 'Rate error down', 'Rate error up',
                             'Rate error average', 'mode'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8',
                               'f8', 'U40')}

    swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    return swiftData


def readSwiftLC_fromSte():
    """
    Read swift lightcurve analysed by Ste.
    """

    swiftData = np.load(swiftDatFileName, encoding='latin1').flat[0]

    headersType = {'names': ('MJD', 'Flux_logpara', 'Flux_ErrL_logpara',
                             'Flux_ErrU_logpara', 'Obs Mode'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'U10')}
    swiftData = np.c_[swiftData['MJD'], swiftData['Flux_logpara'],
                      swiftData['Flux_ErrL_logpara'], swiftData['Flux_ErrU_logpara'],
                      swiftData['Obs Mode']]
    swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    return swiftData


def convertFermiDictToNumpyArray(fermiLC):

    fermiLC = np.c_[fermiLC['tmax_mjd'], fermiLC['tmin_mjd'],
                    fermiLC['flux'], fermiLC['flux_err'],
                    fermiLC['flux_ul95'], fermiLC['ts']]
    headersType = {'names': ('tmax_mjd', 'tmin_mjd', 'flux', 'flux_err',
                             'flux_ul95', 'ts'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
    fermiLC = np.core.records.fromarrays(fermiLC.transpose(), dtype=headersType)

    return fermiLC


def rebinFermi(fermiLC, veritasObs, dates):
    """
    Rebin the Fermi lightcurve given in fermiLC to yearly bins.
    Only contemporaneous observations with VERITAS are taken.
    The VERITAS observation dates are given in veritasObs.
    The range of time for the observations is given in dates,
    a list with two entires [mjd_min, mjd_max].
    """

    fermiMask = list()
    for fermiDataPoint in fermiLC['tmax_mjd']:
        keepFermi = False
        for veritasObsNow in veritasObs:
            if abs(fermiDataPoint - veritasObsNow) < 28:
                keepFermi = True
        fermiMask.append(keepFermi)
    fermiLC = fermiLC[fermiMask]

    nDays = 365
    mjd_min = dates[0]
    mjd_max = dates[1]
    nBins = int((mjd_max - mjd_min)/nDays)
    timeBins = np.linspace(mjd_min, mjd_max, nBins, False)

    rebinnedFermi = defaultdict(list)

    for i_bin, edgeDn in enumerate(timeBins):
        edgeUp = 1e6
        if i_bin < len(timeBins) - 1:
            edgeUp = timeBins[i_bin+1]

        tempFermiData = fermiLC[(edgeDn <= fermiLC['tmax_mjd']) & (fermiLC['tmax_mjd'] < edgeUp)]
        if len(tempFermiData) > 0:
            rebinnedFermi['tmax_mjd'].append(np.average(tempFermiData['tmax_mjd']))
            rebinnedFermi['tmin_mjd'].append(np.average(tempFermiData['tmin_mjd']))
            rebinnedFermi['flux'].append(np.average(tempFermiData['flux'],
                                                    weights=1./tempFermiData['flux_err']))
            rebinnedFermi['flux_err'].append(np.sqrt(np.sum(np.power(tempFermiData['flux_err'],
                                                                     2))))
            rebinnedFermi['flux_ul95'].append(np.average(tempFermiData['flux_ul95']))
            rebinnedFermi['ts'].append(np.average(tempFermiData['ts']))

    fermiLC = np.c_[rebinnedFermi['tmax_mjd'], rebinnedFermi['tmin_mjd'],
                    rebinnedFermi['flux'], rebinnedFermi['flux_err'],
                    rebinnedFermi['flux_ul95'], rebinnedFermi['ts']]
    headersType = {'names': ('tmax_mjd', 'tmin_mjd', 'flux', 'flux_err',
                             'flux_ul95', 'ts'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
    fermiLC = np.core.records.fromarrays(fermiLC.transpose(), dtype=headersType)

    return fermiLC


def rebinSwiftLC(logger, swiftData, binning, veritasObs, dates):

    """
    Rebin the swift lightcurve to the binning set in 'binning'.
    The only allowed values for rebinning are 'weekly', 'monthly' or 'yearly'.
    In the case of weekly and monthly, the observations in the same
    week/month are simply combined. In case of yearly binning, only
    contemporaneous observations with VERITAS are taken. The VERITAS
    observation dates are given in veritasObs.
    The range of time for the observations is given in dates,
    a list with two entires [mjd_min, mjd_max].
    """

    if binning not in ['nightly', 'weekly', 'monthly', 'yearly']:
        logger.critical([['rb', 'I do not know how to rebin',
                                'Swift observations to '],
                         ['c', binning],
                         ['rb', ' binnings, only weekly, monthly or yearly']])

    if binning == 'yearly':
        # Take only contemporaneous observations (in this case within a day)
        swiftMask = list()
        for swiftObsNow in swiftData['Date']:
            keepSwift = False
            for veritasObsNow in veritasObs:
                if abs(swiftObsNow - veritasObsNow) < 1:
                    keepSwift = True
            swiftMask.append(keepSwift)
        swiftData = swiftData[swiftMask]

    nDays = 28
    if binning == 'yearly':
        nDays = 365
    if binning == 'weekly':
        nDays = 7
    if binning == 'nightly':
        nDays = 1
    mjd_min = dates[0]
    mjd_max = dates[1]
    nBins = int((mjd_max - mjd_min)/nDays)
    timeBins = np.linspace(mjd_min, mjd_max, nBins, False)

    date, dateErrDn, dateErrUp = list(), list(), list()
    rate, rateErrDn, rateErrUp = list(), list(), list()
    rateErrAve, mode = list(), list()

    for i_bin, edgeDn in enumerate(timeBins):
        edgeUp = 1e6
        if i_bin < len(timeBins) - 1:
            edgeUp = timeBins[i_bin+1]

        # TODO - should we divide into the different modes?
        tempSwiftData = swiftData[(edgeDn <= swiftData['Date']) & (swiftData['Date'] < edgeUp)]
        if len(tempSwiftData) > 0:
            date.append(np.average(tempSwiftData['Date']))
            dateErrDn.append(date[-1] - np.min(tempSwiftData['Date']))
            dateErrUp.append(np.max(tempSwiftData['Date'] - date[-1]))
            totalError = tempSwiftData['Rate error down'] + tempSwiftData['Rate error up']
            rate.append(np.average(tempSwiftData['Rate'], weights=1./totalError))
            rateErrDn.append(np.sqrt(np.sum(np.power(tempSwiftData['Rate error down'], 2))))
            rateErrUp.append(np.sqrt(np.sum(np.power(tempSwiftData['Rate error up'], 2))))
            rateErrAve.append(np.sqrt(np.sum(np.power(tempSwiftData['Rate error average'], 2))))
            mode.append('Combined')

    swiftData = np.c_[date, dateErrDn, dateErrUp,
                      rate, rateErrDn, rateErrUp,
                      rateErrAve, mode]
    headersType = {'names': ('Date', 'Date error down', 'Date error up',
                             'Rate', 'Rate error down', 'Rate error up',
                             'Rate error average', 'mode'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8',
                               'f8', 'U40')}
    swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    return swiftData


def readCorrTable(corrTableFile):
    """
    Read the table of corrections derived from the Crab lightcurve
    """

    headersType = {'names': ('Left edges', 'Right edges',
                             'Correction factor', 'CorrFactorError',
                             'CorrFactorErrorCons'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8')}

    return np.loadtxt(corrTableFile, dtype=headersType)


def correctFluxesFromCrabLC(origLC, corrTable):
    """
    Scale points in the lighcurve based on the table of
    corrections derived from the Crab lightcurve.
    """

    corrLC = np.copy(origLC)

    for i_point, dateNow in enumerate(corrLC['DateMJD']):
        corrBin = np.argmax(dateNow < corrTable['Right edges'])
        if corrTable['Correction factor'][corrBin] != 1:
            corrLC['Flux'][i_point] = (corrLC['Flux'][i_point] /
                                       corrTable['Correction factor'][corrBin])
            corrLC['Flux Error'][i_point] = np.sqrt(np.power(corrLC['Flux Error'][i_point], 2) +
                                                    np.power(corrTable['CorrFactorError'][corrBin] *
                                                             corrLC['Flux'][i_point], 2))

    return corrLC


def lcBetweenDates(logger, lc, dates, experiment):

    if dates[0] >= dates[1]:
        logger.critical([['rb', 'Dates to choose the LC in must make sense']])

    dateFormat = 'DateMJD'
    if experiment.lower() == 'veritas':
        logger.info([['wb', 'Extracting VERITAS LC between'],
                     ['c', dates[0]],
                     ['wb', 'and'],
                     ['c', dates[1]]])
        return lc[(lc[dateFormat] >= dates[0]) & (lc[dateFormat] <= dates[1])]
    elif experiment.lower() == 'fermi':
        logger.info([['wb', 'Extracting Fermi LC between'],
                     ['c', dates[0]],
                     ['wb', 'and'],
                     ['c', dates[1]]])
        dateFormat = 'tmax_mjd'
        return lc[(lc[dateFormat] >= dates[0]) & (lc[dateFormat] <= dates[1])]
    elif experiment.lower() == 'swift':
        logger.info([['wb', 'Extracting Swift LC between'],
                     ['c', dates[0]],
                     ['wb', 'and'],
                     ['c', dates[1]]])
        dateFormat = 'Date'
        return lc[(lc[dateFormat] >= dates[0]) & (lc[dateFormat] <= dates[1])]
    else:
        logger.critical([['rb', 'Can choose lc between dates only for VERITAS, Fermi and Swift']])
