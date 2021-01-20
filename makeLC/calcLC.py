#!/usr/bin/python

import os
import os.path
import numpy as np
import ROOT
import argparse
import subprocess
import resource
import glob
import sys
sys.path.append("..")
from lib import colourLogger
from lib import generalUtil as gUtil


def calcLC(source, configAnal):

    logStdout.info([['wb', 'Producing a lightcurve for'],
                    ['p', source]])

    anasumFile = configAnal['anasumFiles'][source]
    sourceSpectralParam = configAnal['spectralParameters'][source]
    dirNow = os.path.join(os.getcwd(), source)

    minEnergyFluxCalc = sourceSpectralParam[0]
    normalizationEnergy = sourceSpectralParam[1]
    spectralIndex = sourceSpectralParam[2]
    maxEnergyFluxCalc = sourceSpectralParam[3]

    # Read energy threshold from file
    headersType = {'names': ('Source', 'z', 'tau = 1',
                             'tau = 2', 'tau = 3'),
                   'formats': ('U20', 'f8',
                               'f8', 'f8', 'f8')}

    sourceEnergyThresholds = np.loadtxt('../../{}_sourcesThresholds.txt'.
                                        format(configAnal['EBL']['model']), dtype=headersType)
    sourceIndex = np.where(sourceEnergyThresholds['Source'] == source)
    # first two entries are source and z
    thresholdsForThisSource = list(sourceEnergyThresholds[sourceIndex[0][0]])[2:4]
    thresholdsForThisSource.insert(0, minEnergyFluxCalc)

    binningDict = {'nightly': 1, 'weekly': 7, 'monthly': 28, 'yearly': 365}
    mjd_min = configAnal['dates']['veritas'][0]
    mjd_max = configAnal['dates']['veritas'][1]

    thresholds = list()
    # Build a list of pairs of min/max energies to run over
    for i_thr, thresholdNow in enumerate(thresholdsForThisSource):

        minEnergy = thresholdNow
        if i_thr == len(thresholdsForThisSource) - 1:
            maxEnergy = 30.
        else:
            maxEnergy = thresholdsForThisSource[i_thr + 1]

        prefix = 'LightCurve_'
        interfix = 'tau-{}_'.format(i_thr + 1)
        suffix = 'range_{0:1.0f}_{1:1.0f}_GeV.txt'.format(minEnergy*1000., maxEnergy*1000.)
        fileName = prefix + interfix + suffix
        thresholds.append({'minEnergy': thresholdNow,
                           'maxEnergy': maxEnergy,
                           'fileName': fileName})

    # Make two versions of the lightcurves, one with negative fluxes
    # and one without. The former is only for the luminosity function study.
    suffixNegativeFlux = ['', '_unbound']

    for suffixNegativeFluxNow in suffixNegativeFlux:

        prefix = 'LightCurve_'
        interfix = 'fullEnergyRange_'
        suffix = 'range_{0:1.0f}_{1:1.0f}_GeV{2}.txt'.format(minEnergyFluxCalc*1000.,
                                                             30*1000.,
                                                             suffixNegativeFluxNow)
        fileName = prefix + interfix + suffix
        thresholds.append({'minEnergy': minEnergyFluxCalc,
                           'maxEnergy': maxEnergy,
                           'fileName': fileName})

    for thresholdsNow in thresholds:

        for binning, timeBinLength_days in binningDict.items():

            logStdout.info([['p', source],
                            ['bb', binning],
                            ['wb', 'lightcurve,'],
                            ['g', '{} < E < {} TeV'.format(thresholdsNow['minEnergy'],
                                                           thresholdsNow['maxEnergy'])]])

            iLightCurve = ROOT.VLightCurve()

            iLightCurve.initializeTeVLightCurve(anasumFile, timeBinLength_days,
                                                mjd_min, mjd_max)

            # set spectral parameters (start of flux calculation [TeV],
            # normalization at energy E0 [TeV], spectral index, end of flux calculation [TeV])
            iLightCurve.setSpectralParameters(minEnergyFluxCalc, normalizationEnergy,
                                              spectralIndex, maxEnergyFluxCalc)

            # plot ULs for poins with <2 sigma significance; to get rid of
            # significance limits for upper limits of flux use (-999, 999);
            iLightCurve.setSignificanceParameters(-999, -999)
            # Avoid negative fluxes (if True it avoids them)
            iLightCurve.setFluxCalculationMethod('unbound' not in thresholdsNow['fileName'])
            # calculate fluxes and upper flux limits for energies in the range given
            iLightCurve.fill(thresholdsNow['minEnergy'], thresholdsNow['maxEnergy'])
            iLightCurve.writeASCIIFile(os.path.join(dirNow,
                                                    binning +
                                                    thresholdsNow['fileName']))
            del iLightCurve


def moveFileToLumiFunctionStudyDir(source):

    files = glob.glob('{}/*fullEnergyRange*'.format(source))
    for fileNow in files:
        if 'unbound' in fileNow:
            command = 'mv'
            fileDest = fileNow.replace('_unbound', '')
        else:
            command = 'cp'
            fileDest = fileNow.replace('.txt', '_bounded.txt')

        fileDest = fileDest.replace('fullEnergyRange_range', 'energyRange')
        # The 195 GeV threshold in ED is essentially equivalent to 200 GeV.
        fileDest = fileDest.replace('195', '200')

        subprocess.call('mkdir -p lightcurvesForLumiFunctionStudy/{}'.format(source), shell=True)
        fullCommand = '/bin/{} {} lightcurvesForLumiFunctionStudy/{}'.format(command,
                                                                             fileNow,
                                                                             fileDest)
        subprocess.call(fullCommand, shell=True)

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Produce light-curve.')
    parser.add_argument('source', action='store',
                        help='Source to produce the light-curve for')
    parser.add_argument('configAnalFile', action='store',
                        help='YAML file with the analysis configuration')

    args = parser.parse_args()

    logStdout = colourLogger.initStdOutLogger()

    # A horrible hack to avoid the following ROOT error
    # SysError in <TFile::TFile>: file  can not be opened for reading (Too many open files)
    logStdout.info([['wb', 'getrlimit before:'],
                    ['bb', resource.getrlimit(resource.RLIMIT_NOFILE)]])
    resource.setrlimit(resource.RLIMIT_NOFILE, (4096, 4096))
    logStdout.info([['wb', 'getrlimit after:'],
                    ['bb', resource.getrlimit(resource.RLIMIT_NOFILE)]])

    ROOT.gROOT.SetBatch(True)
    ROOT.TFile.Open._creates = True

    # load shared library
    ROOT.gSystem.Load("$EVNDISPSYS/lib/libVAnaSum.so")

    configAnal = gUtil.readYamlFile(logStdout, args.configAnalFile)
    calcLC(args.source, configAnal)

    moveFileToLumiFunctionStudyDir(args.source)
