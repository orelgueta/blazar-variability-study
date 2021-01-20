#!/usr/bin/python

import numpy as np
from astropy import units as u
from collections import defaultdict
import inspect
import yaml
from lib import generalUtil as gUtil


def readFranceschini(logger, fileName):

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    eblDict = dict()
    modelFile = open(fileName, 'r')
    for line in modelFile:
        if 'redshift' in line:
            redshiftNow = round(float(line.split()[3]), 3)
            eblDict[redshiftNow] = defaultdict(list)
        elif 'TeV' not in line:
            eblDict[redshiftNow]['E [TeV]'].append(float(line.split()[0]))
            eblDict[redshiftNow]['E [eV]'].append(float(line.split()[1]))
            eblDict[redshiftNow]['tau'].append(float(line.split()[2]))
            eblDict[redshiftNow]['e^tau'].append(float(line.split()[3]))

    return eblDict


def readVeritasSpectrumFromLog(logger, fileName):
    """
    Read a VERITAS SED from a log file.
    Very strange that we have to do that,
    but it's the easiest way I found at the moment to get all of the information.
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    diffFluxTempDict = defaultdict(list)
    logFile = open(fileName, 'r')
    readNow = False
    for line in logFile:
        if line.startswith('# TeV'):
            readNow = True
            continue
        if readNow and line.startswith('#---------------------'):
            readNow = False
            continue
        if readNow:
            diffFluxData = line.strip().split()
            diffFluxTempDict['E'].append(float(diffFluxData[0]))
            diffFluxTempDict['E_min'].append(float(diffFluxData[1]))
            diffFluxTempDict['E_max'].append(float(diffFluxData[2]))
            diffFluxTempDict['dE'].append(float(diffFluxData[3])/2.)
            diffFluxTempDict['dN/dE'].append(float(diffFluxData[4]))
            diffFluxTempDict['dN/dE error'].append(float(diffFluxData[5]))
            if '(' in diffFluxData[6]:
                diffFluxTempDict['dN/dE error down'].append(float(diffFluxData[6].
                                                                  replace('(', '').
                                                                  replace(',', '')))
                diffFluxTempDict['dN/dE error up'].append(float(diffFluxData[7].
                                                                replace(')', '')))
                columnNow = 8
            else:
                diffFluxTempDict['dN/dE error down'].append(0.)
                diffFluxTempDict['dN/dE error up'].append(0.)
                columnNow = 7

            diffFluxTempDict['nOn'].append(float(diffFluxData[columnNow]))
            diffFluxTempDict['nOn error'].append(float(diffFluxData[columnNow + 1]))
            diffFluxTempDict['nOff'].append(float(diffFluxData[columnNow + 2]))
            diffFluxTempDict['nOff error'].append(float(diffFluxData[columnNow + 3]))
            diffFluxTempDict['alpha'].append(float(diffFluxData[columnNow + 4]))
            diffFluxTempDict['sigma'].append(float(diffFluxData[columnNow + 5]))
            diffFluxTempDict['ts'].append(float(diffFluxData[columnNow + 6]))

    veritasData = np.c_[diffFluxTempDict['E'], diffFluxTempDict['E_min'],
                        diffFluxTempDict['E_max'], diffFluxTempDict['dE'],
                        diffFluxTempDict['dN/dE'], diffFluxTempDict['dN/dE error'],
                        diffFluxTempDict['dN/dE error down'], diffFluxTempDict['dN/dE error up'],
                        diffFluxTempDict['nOn'], diffFluxTempDict['nOn error'],
                        diffFluxTempDict['nOff'], diffFluxTempDict['nOff error'],
                        diffFluxTempDict['alpha'], diffFluxTempDict['sigma'],
                        diffFluxTempDict['ts']]

    headersType = {'names': ('E', 'E_min', 'E_max', 'dE',
                             'dN/dE', 'dN/dE error',
                             'dN/dE error down', 'dN/dE error up',
                             'nOn', 'nOn error', 'nOff', 'nOff error',
                             'alpha', 'sigma', 'ts'),
                   'formats': ('f8', 'f8', 'f8', 'f8',
                               'f8', 'f8', 'f8', 'f8',
                               'f8', 'f8', 'f8', 'f8',
                               'f8', 'f8', 'f8')}

    veritasData = np.core.records.fromarrays(veritasData.transpose(), dtype=headersType)

    return veritasData


def readVeritasSpectrumFromTable(logger, fileName):
    """
    Read a VERITAS SED from a file containing the table (prepared manually).
    The structure of the table is seen below.
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    headersType = {'names': ('E', 'E_min', 'E_max', 'dE',
                             'dN/dE', 'dN/dE error',
                             'dN/dE error down', 'dN/dE error up'),
                   'formats': ('f8', 'f8', 'f8', 'f8',
                               'f8', 'f8',
                               'U40', 'U40')}

    veritasData = np.loadtxt(veritasFile, dtype=headersType)

    veritasData['dN/dE error down'] = [entry.replace('(', '').replace(',', '')
                                       for entry in veritasData['dN/dE error down']]
    veritasData['dN/dE error up'] = [entry.replace(')', '')
                                     for entry in veritasData['dN/dE error up']]
    return veritasData


def readVeritasData(logger, fileName, fromLog=True):
    """
    Read the VERITAS data either from a log file or a table (prepared manually).
    When reading from a log file, the spectral points to use are selected
    based on the logic explained in selectSpectralPointsVeritas.
    Returns the VERITAS spectrum after selecting the appropriate points.
    """

    if fromLog:
        veritasData = readVeritasSpectrumFromLog(logger, fileName)
        veritasData = selectSpectralPointsVeritas(logger, veritasData)
        return veritasData
    else:
        return readVeritasSpectrumFromTable(logger, fileName)


def selectSpectralPointsVeritas(logger, veritasData):
    """
    Choose which flux points to use based on the following argument from the EBL wiki:

        The question of which spectral points to include has haunted VERITAS for generations.
        We adopt the following prescription:
        include one bin past the last significant (2 sigma) bin OR
        the last bin with ON>0 by 2 sigma, whichever has the higher energy.
        This prevents a bias towards harder spectral fits,
        and ensures that well-measured spectral cutoffs are accounted for.
    """

    lastPointSigma = 0
    lastPointNon = 0
    for i_point, sigmaNow in enumerate(veritasData['sigma']):
        # Check the significance of this point
        if sigmaNow < 2 and lastPointSigma == 0:
            lastPointSigma = i_point + 1

        # Check the number of ON counts of this point
        nOnNow = veritasData['nOn'][i_point] - 2*veritasData['nOn error'][i_point]
        # If this is the first time we find ON - ON_err < 0, take the previous bin
        # to be the one fulfilling the condition stated above
        if nOnNow < 0 and lastPointNon == 0:
            lastPointNon = i_point - 1

        # Break when find both options
        if lastPointNon > 0 and lastPointSigma > 0:
            break

    lastPoint = lastPointSigma if lastPointSigma > lastPointNon else lastPointNon
    pointsToPlot = [True if i_point < lastPoint else False for i_point in range(len(veritasData))]
    # Test also if the first point has sigma > 2 and ON>0 by 2 sigma
    if veritasData['sigma'][0] < 2 or (veritasData['nOn'][0] - 2*veritasData['nOn error'][0]) < 0:
        pointsToPlot[0] = False

    veritasData = veritasData[pointsToPlot]

    return veritasData


def readFermiData(logger, fileName):

    """
    Read a Fermi SED from a numpy file saved by FermiPy
    """

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    sedFermi = np.load(fileName, encoding='latin1', allow_pickle=True).flat[0]

    return convertFermiDictToNumpyArray(sedFermi)


def convertFermiDictToNumpyArray(sedFermi):

    sedFermi = np.c_[sedFermi['e_ctr'], sedFermi['dnde'], sedFermi['e2dnde'], sedFermi['dnde_err'],
                     sedFermi['dnde_err_lo'], sedFermi['dnde_err_hi'],
                     sedFermi['e2dnde_err_lo'], sedFermi['e2dnde_err_hi']]
    headersType = {'names': ('e_ctr', 'dnde', 'e2dnde', 'dnde_err',
                             'dnde_err_lo', 'dnde_err_hi',
                             'e2dnde_err_lo', 'e2dnde_err_hi'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
    sedFermi = np.core.records.fromarrays(sedFermi.transpose(), dtype=headersType)

    return sedFermi


def readSwiftData(logger, fileName, swiftType='online'):

    gUtil.checkFileExists(logger, fileName, inspect.stack()[0][3])

    if swiftType == 'online':
        return readSwiftDataOnline(logger, fileName)
    elif swiftType == 'fromSte':
        return readSwiftDataSte(logger, fileName)
    elif swiftType == 'mireia':
        return readSwiftDataMireia(logger, fileName)
    else:
        logger.critical([['rb', 'I do not know what type of swift data is'],
                         ['c', swiftType]])
        return


def readSwiftDataOnline(logger, fileName):

    return readXrtOnline(logger, fileName), list()


def readXrtOnline(logger, fileName):

    headersType = {'names': ('E', 'dE',
                             'E^2 dN/dE', 'E^2 dN/dE error',
                             'E^2 dN/dE model'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8')}

    xrtData = np.loadtxt(fileName, dtype=headersType, skiprows=3)

    tempXrtData = dict()
    kevToGeV = (1*u.keV).to('GeV')

    tempXrtData['E'] = xrtData['E']*kevToGeV
    tempXrtData['dE down'] = xrtData['dE']*kevToGeV
    tempXrtData['dE up'] = xrtData['dE']*kevToGeV
    tempXrtData['E^2 dN/dE'] = xrtData['E^2 dN/dE']*kevToGeV
    tempXrtData['E^2 dN/dE error down'] = xrtData['E^2 dN/dE error']*kevToGeV
    tempXrtData['E^2 dN/dE error up'] = xrtData['E^2 dN/dE error']*kevToGeV

    xrtData = np.c_[tempXrtData['E'], tempXrtData['dE down'], tempXrtData['dE up'],
                    tempXrtData['E^2 dN/dE'], tempXrtData['E^2 dN/dE error down'],
                    tempXrtData['E^2 dN/dE error up']]

    headersType = {'names': ('E', 'dE down', 'dE up',
                             'E^2 dN/dE', 'E^2 dN/dE error down', 'E^2 dN/dE error up'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8')}

    xrtData = np.core.records.fromarrays(xrtData.transpose(), dtype=headersType)

    return xrtData


def readSwiftDataSte(logger, fileName):

    headersType = {'names': ('E', 'dE',
                             'E^2 dN/dE', 'E^2 dN/dE error'),
                   'formats': ('f8', 'f8', 'f8', 'f8')}

    xrtData = np.loadtxt(fileName, dtype=headersType)

    tempXrtData = dict()
    tempXrtData['dE down'] = xrtData['dE']
    tempXrtData['dE up'] = xrtData['dE']
    tempXrtData['E^2 dN/dE error down'] = xrtData['E^2 dN/dE error']
    tempXrtData['E^2 dN/dE error up'] = xrtData['E^2 dN/dE error']

    xrtData = np.c_[xrtData['E'], tempXrtData['dE down'], tempXrtData['dE up'],
                    xrtData['E^2 dN/dE'], tempXrtData['E^2 dN/dE error down'],
                    tempXrtData['E^2 dN/dE error up']]

    headersType = {'names': ('E', 'dE down', 'dE up',
                             'E^2 dN/dE', 'E^2 dN/dE error down', 'E^2 dN/dE error up'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8')}

    xrtData = np.core.records.fromarrays(xrtData.transpose(), dtype=headersType)

    return xrtData, list()


def readSwiftDataMireia(logger, fileName):

    with open(fileName, 'r') as fSwift:
        swiftDataTemp = yaml.safe_load(fSwift)

    swiftData = np.c_[swiftDataTemp['x'], swiftDataTemp['xel'], swiftDataTemp['xeh'],
                      swiftDataTemp['y'], swiftDataTemp['yel'], swiftDataTemp['yeh'],
                      swiftDataTemp['y_obs'], swiftDataTemp['yel_obs'],
                      swiftDataTemp['yeh_obs'], swiftDataTemp['inst']]

    headersType = {'names': ('E', 'dE down', 'dE up',
                             'E^2 dN/dE', 'E^2 dN/dE error down', 'E^2 dN/dE error up',
                             'obs E^2 dN/dE', 'obs E^2 dN/dE error down',
                             'obs E^2 dN/dE error up', 'inst'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8',
                               'f8', 'f8', 'f8', 'U40')}
    swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    xrtData = swiftData[swiftData['inst'] == 'XRT']
    ergToMeV = (1*u.erg).to('GeV')
    xrtData['E^2 dN/dE'] = xrtData['E^2 dN/dE']*ergToMeV
    xrtData['E^2 dN/dE error down'] = xrtData['E^2 dN/dE error down']*ergToMeV
    xrtData['E^2 dN/dE error up'] = xrtData['E^2 dN/dE error up']*ergToMeV

    uvotData = swiftData[swiftData['inst'] == 'UVOT']
    ergToMeV = (1*u.erg).to('GeV')
    uvotData['E^2 dN/dE'] = uvotData['E^2 dN/dE']*ergToMeV
    uvotData['E^2 dN/dE error down'] = uvotData['E^2 dN/dE error down']*ergToMeV
    uvotData['E^2 dN/dE error up'] = uvotData['E^2 dN/dE error up']*ergToMeV

    return xrtData, uvotData
