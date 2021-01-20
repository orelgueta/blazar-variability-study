#!/usr/bin/python

import os.path
import os
import glob
import subprocess
import numpy as np
import numpy.lib.recfunctions as rfn
from astropy.io import fits
from astropy.stats import bayesian_blocks
import argparse
from collections import defaultdict


def checkDatFile(datFileName):
    if not os.path.isfile(datFileName):
        print(datFileName + ' is not a file!\n')
        return False
    return True


def extractErrors(errorStr):
    errors = errorStr.replace('(', '').replace(')', '').split(' - ')
    return float(errors[0]), float(errors[1])


def calibrate_ncp_prior(flux=None, fluxerr=None, time=None, timebin=None,
                        p_0=[0.05], n_sims=1000, min_prior=0.2, max_prior=4,
                        n_steps=20, outPrefix=None):
                        # path='./gammas/', exp='VERITAS', source=''):

    # Calibration of ncp_prior:
    # input:
    #        flux, fluxerr, time, timebin : Lightcurve in format numpy.ndarray or pandas.Series
    #        p_0 : FPR input array
    #        n_sims : float
    #        min_prior : float/int
    #        max_prior : float/int
    #        n_stepts : number of steps in [min_prior, max_prior]

    sourceNow = outPrefix.split('/')[0]

    falsecount = np.zeros(n_steps)
    ncp_priors = np.linspace(min_prior, max_prior, n_steps)
    result = {}
    best = {}

    # distance between points not relevant but should be ordered
    x = np.arange(len(flux))
    average = np.average(flux, weights=fluxerr)

    # simulating lightcurves for n_sims times and applying algorithem
    # in n_steps steps between min_prior and max_prior. Afterwards
    # false positive rate is calculated if a block was detected.
    for k in range(n_sims):
        if k % 10 == 0:
            print(sourceNow, 'current simulation: {}'.format(k))

        # simulate the flux values
        datapoints = np.random.normal(average, fluxerr, len(fluxerr))

        # aply bayesian block and count fpr
        for l, ncp_prior in enumerate(ncp_priors):
            gamma = 10**(-ncp_prior)
            bb = bayesian_blocks(x, datapoints, fluxerr, fitness='measures', gamma=gamma)
            if len(bb) > 2:
                falsecount[l] += 1

    fp_rate = falsecount/n_sims

    # Final result of FPR in dependency of ncp_prior
    result = np.core.records.fromarrays([ncp_priors, fp_rate], names='ncp, fp')

    # Calculation of best results for the values in p_0
    for p0 in p_0:
        best[str(p0)] = result[(np.abs(result.fp - p0)).argmin()]

    # Saving result and best to txt file
    with open(outPrefix + '_result.txt', 'wb') as fOut:
        np.savetxt(fOut, result)

    # with open(outPrefix + '_results_best.txt', 'wb') as fOut:
    #     np.savetxt(fOut, [best])

    return(result, best)


def readSwiftLC(swiftFileName, rebin, veritasObs):

    swiftFile = open(swiftFileName, 'r')
    date, dateErrUp, dateErrDn = list(), list(), list()
    rate, rateErrUp, rateErrDn = list(), list(), list()
    mode = list()
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
            mode.append(modeNow)

    swiftData = np.c_[date, dateErrDn, dateErrUp,
                      rate, rateErrDn, rateErrUp,
                      mode]

    headersType = {'names': ('Date', 'Date error down', 'Date error up',
                             'Rate', 'Rate error down', 'Rate error up',
                             'mode'),
                   'formats': ('f8', 'f8', 'f8',
                               'f8', 'f8', 'f8',
                               'U40')}

    swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    if rebin == 'monthly' or rebin == 'weekly' or rebin == 'yearly':

        if rebin == 'yearly':
            # Take only contemporaneous observations
            swiftMask = list()
            for swiftObsNow in swiftData['Date']:
                keepSwift = False
                for veritasObsNow in veritasObs:
                    if abs(swiftObsNow - veritasObsNow) < 1:
                        keepSwift = True
                swiftMask.append(keepSwift)
            swiftData = swiftData[swiftMask]

        nDays = 28
        if rebin == 'yearly':
            nDays = 365
        if rebin == 'weekly':
            nDays = 7
        mjd_min = 53423  # This is exactly 147 weeks before the start day of Fermi
        mjd_max = 58465  # This is ~today
        nBins = int((mjd_max - mjd_min)/nDays)
        timeBins = np.linspace(mjd_min, mjd_max, nBins, False)

        date, dateErrDn, dateErrUp = list(), list(), list()
        rate, rateErrDn, rateErrUp = list(), list(), list()
        mode = list()

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
                mode.append('Combined')

        swiftData = np.c_[date, dateErrDn, dateErrUp,
                          rate, rateErrDn, rateErrUp,
                          mode]
        swiftData = np.core.records.fromarrays(swiftData.transpose(), dtype=headersType)

    return swiftData


def rebinFermi(fermiLC, veritasObs):

    # First convert to numpy array to make it easier
    fermiLC = np.c_[fermiLC['tmax_mjd'], fermiLC['tmin_mjd'], fermiLC['flux'], fermiLC['flux_err']]
    headersType = {'names': ('tmax_mjd', 'tmin_mjd', 'flux', 'flux_err'),
                   'formats': ('f8', 'f8', 'f8', 'f8')}
    fermiLC = np.core.records.fromarrays(fermiLC.transpose(), dtype=headersType)

    # Take only contemporaneous observations (in this case, within a month)
    # fermiBlocks = bayesian_blocks(fermiLC['tmax_mjd'], fermiLC['flux'], fermiLC['flux_err']
    fermiMask = list()
    for fermiDataPoint in fermiLC['tmax_mjd']:
        keepFermi = False
        for veritasObsNow in veritasObs:
            if abs(fermiDataPoint - veritasObsNow) < 28:
                keepFermi = True
        fermiMask.append(keepFermi)
    fermiLC = fermiLC[fermiMask]

    nDays = 365
    mjd_min = 53423  # This is exactly 147 weeks before the start day of Fermi
    mjd_max = 58465  # This is ~today
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

    fermiLC = np.c_[rebinnedFermi['tmax_mjd'], rebinnedFermi['tmin_mjd'],
                    rebinnedFermi['flux'], rebinnedFermi['flux_err']]
    fermiLC = np.core.records.fromarrays(fermiLC.transpose(), dtype=headersType)

    return fermiLC


def readCorrTable(corrTableFile):

    headersType = {'names': ('Left edges', 'Right edges',
                             'Correction factor', 'CorrFactorError',
                             'CorrFactorErrorCons'),
                   'formats': ('f8', 'f8', 'f8', 'f8', 'f8')}

    return np.loadtxt(corrTableFile, dtype=headersType)


def correctFluxesFromCrabLC(origLC, corrTable):

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


def correctFluxes(origLC, corrTable):

    corrLC = correctFluxesFromCrabLC(origLC, corrTable)
    # We increased the threshold, so no need to add a systematic uncertainty anymore

    return corrLC


def determinePriors(veritasDatFileName, fermiFile, swiftFullFileName, corrTable,
                    veritasObsFile, sourceNow, binning):

    for fileNow in [veritasDatFileName, fermiFile, swiftFullFileName]:
        if not checkDatFile(fileNow):
            return

    veritasDatFile = open(veritasDatFileName, 'r')
    headersType = {'names': ('DateMJD', 'Date Error',
                             'Flux', 'Flux Error'),
                   'formats': ('f8', 'f8',
                               'f8', 'f8')}

    veritasData = np.loadtxt(veritasDatFile, dtype=headersType)
    veritasFluxes = veritasData[veritasData['Flux Error'] > 0]

    veritasFluxes = correctFluxes(veritasFluxes, corrTable)

    nsims = 15000
    n_steps = 40
    experiment = 'veritas'
    outPrefix = '{}/{}_{}_{}'.format(sourceNow,
                                     experiment,
                                     binning,
                                     str(nsims))

    result, best = calibrate_ncp_prior(flux=veritasFluxes['Flux'],
                                       fluxerr=veritasFluxes['Flux Error'],
                                       time=veritasFluxes['DateMJD'],
                                       timebin=veritasFluxes['Date Error'],
                                       p_0=[0.01, 0.05], n_sims=nsims,
                                       min_prior=0.2, max_prior=4,
                                       n_steps=n_steps, outPrefix=outPrefix)
    gamma = 10**(- best[str(0.01)].ncp)
    print(sourceNow, 'VERITAS', 'gamma - ', gamma)

    if binning == 'yearly' or binning == 'monthly':

        fermiDatFile = open(fermiFile, 'rb')
        fermiLC = np.load(fermiDatFile, encoding='latin1').flat[0]

        if binning == 'yearly':
            headersType = {'names': ('run', 'date', 'flux', 'fluxError',
                                     'significance', 'ze'),
                           'formats': ('f8', 'f8', 'f8', 'f8',
                                       'f8', 'f8')}
            veritasObs = np.loadtxt(veritasObsFile, dtype=headersType)
            fermiLC = rebinFermi(fermiLC, veritasObs['date'])

        experiment = 'fermi'
        outPrefix = '{}/{}_{}_{}'.format(sourceNow,
                                         experiment,
                                         binning,
                                         str(nsims))
        result, best = calibrate_ncp_prior(flux=fermiLC['flux'],
                                           fluxerr=fermiLC['flux_err'],
                                           time=fermiLC['tmax_mjd'],
                                           timebin=fermiLC['tmax_mjd'] - fermiLC['tmin_mjd'],
                                           p_0=[0.01, 0.05], n_sims=nsims,
                                           min_prior=0.2, max_prior=4,
                                           n_steps=n_steps, outPrefix=outPrefix)
        gamma = 10**(- best[str(0.01)].ncp)
        print(sourceNow, 'Fermi', 'gamma - ', gamma)

    swiftBinnings = [binning]
    if binning == 'yearly':  # run also the daily for Swift in this case
        swiftBinnings = ['daily', 'yearly']

    for swiftBinNow in swiftBinnings:

        if swiftBinNow == 'yearly':
            veritasObsDates = veritasObs['date']
        else:
            veritasObsDates = list()

        swiftData = readSwiftLC(swiftFile, swiftBinNow, veritasObsDates)

        experiment = 'swift'
        outPrefix = '{}/{}_{}_{}'.format(sourceNow,
                                         experiment,
                                         swiftBinNow,
                                         str(nsims))
        swiftRateErrorAverage = (swiftData['Rate error down'] + swiftData['Rate error up'])/2.
        result, best = calibrate_ncp_prior(flux=swiftData['Rate'],
                                           fluxerr=swiftRateErrorAverage,
                                           time=swiftData['Date'],
                                           timebin=(swiftData['Date error down'] +
                                                    swiftData['Date error up']),
                                           p_0=[0.01, 0.05], n_sims=nsims,
                                           min_prior=0.2, max_prior=4,
                                           n_steps=n_steps, outPrefix=outPrefix)
        gamma = 10**(- best[str(0.01)].ncp)
        print(sourceNow, 'Swift', swiftBinNow, 'gamma - ', gamma)

    return


if __name__ == '__main__':

    np.random.seed(1234)

    parser = argparse.ArgumentParser(description=('Calculate optimal '
                                                  'priors for Bayesian blocks.'))
    parser.add_argument('source')
    parser.add_argument('binning')

    args = parser.parse_args()

    sources = {'1ES0033': '1ES 0033+595',
               '1ES0502': '1ES 0502+675',
               '1ES1011': '1ES 1011+496',
               '1ES1218': '1ES 1218+304',
               '1ES0229': '1ES 0229+200',
               'RGBJ0710': 'RGB J0710+591',
               'PG1553':  'PG 1553+113',
               'PKS1424': 'PKS 1424+240'
               }

    if args.source not in sources:
        print('Source', args.source, 'not known')

    hdulist = fits.open(('/afs/ifh.de/group/cta/scratch/ogueta/sw/anaconda/envs/fermi/'
                         'lib/python2.7/site-packages/fermipy/data/catalogs/gll_psc_8year_v5.fit'))
    sourceCatalog = hdulist[1].data

    workDir = os.getcwd() + '/'
    fermiPrefix = '/lustre/fs19/group/cta/users/ogueta/fermi/variabilityStudy/'
    veritasPrefix = '/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/makeLC/'
    swiftPrefix = '/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/swift/onlineTool/'

    corrTableFile = ('/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/'
                     'crabStability/plotLC/correctionFactors.txt')

    veritasObsPrefix = '/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/spectra/'

    for i_src, sourceTeV in enumerate(sourceCatalog['ASSOC_TEV']):
        if sources[args.source] in sourceTeV:
            fermiLC = sourceCatalog['Source_Name'][i_src].replace(' ', '_').lower()
    fermiLC += '_lightcurve.npy'

    fermiBinning = args.binning
    if fermiBinning != 'monthly':
        fermiBinning = 'monthly'
    fermiFile = os.path.join(fermiPrefix, '{}LightCurves'.format(fermiBinning),
                             args.source, args.source, fermiLC)

    veritasDirectory = os.path.join(veritasPrefix, args.source)
    veritasLC = glob.glob(os.path.join(veritasDirectory,
                          '{}*fullEnergyRange*.txt'.format(args.binning)))[0]
    veritasFile = os.path.join(veritasDirectory, veritasLC)
    corrTable = readCorrTable(corrTableFile)
    veritasObsFile = os.path.join(os.path.join(veritasObsPrefix, args.source), 'fluxPerRun.txt')

    swiftFile = os.path.join(swiftPrefix, args.source,
                             'dailyBins', '{}_lightcurve.qdp'.format(args.source))

    try:
        subprocess.check_call(['mkdir', '-p', args.source])
    except subprocess.CalledProcessError as e:
        print('Could not create output directory')
        sys.exit(1)

    determinePriors(veritasFile, fermiFile, swiftFile, corrTable,
                    veritasObsFile, args.source, args.binning)
