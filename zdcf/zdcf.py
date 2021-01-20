#!/usr/bin/python

import os.path
import os
import sys
import glob
import subprocess
import platform
import numpy as np
import matplotlib as mlp
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import periodogram
from scipy import stats
from scipy.stats import gaussian_kde
from astropy.time import Time
from astropy.stats import bayesian_blocks
from astroML.fourier import PSD_continuous
from astroML.time_series import generate_power_law
from collections import defaultdict
from collections import OrderedDict
import argparse
import time
import datetime


sys.path.append("..")
from lib import colourLogger
from lib import generalUtil as gUtil
from lib import lcUtil
from lib import plotUtil as pUtil


def psd(x, *p):
    norm, beta = p
    return norm/np.power(x, beta)


def readLCs(file1, file2):

    headersType = {'names': ('mjd', 'flux', 'fluxErr'),
                   'formats': ('f8', 'f8', 'f8')}

    lc1 = np.loadtxt(file1, delimiter=',', dtype=headersType)
    lc2 = np.loadtxt(file2, delimiter=',', dtype=headersType)

    return lc1, lc2


def buildLC(lc, fluxes, sampling, roll):

    headersType = {'names': ('mjd', 'flux', 'fluxErr'),
                   'formats': ('f8', 'f8', 'f8')}

    # Deal with unreasonable uncertainties
    relUncert = (lc['fluxErr']/lc['flux'])
    reasonUncert = relUncert[relUncert < 1]
    meanUncert = np.mean(reasonUncert)
    relUncert[relUncert > 1] = meanUncert

    # To avoid using always the same MJD values and
    # as a result having empty lag values in the zDCF calculation
    sInDay = 60*60*24
    lag = int(np.round(np.random.normal(0, roll*sampling/sInDay, 1)[0]))
    newTime = lc['mjd'] - lag
    newTime = np.roll(newTime, lag)
    newFlux = np.roll(fluxes, lag)
    newFluxErr = np.roll(relUncert*fluxes, lag)

    newLC = np.c_[newTime, newFlux, newFluxErr]
    # newLC = np.c_[lc['mjd'], fluxes, relUncert*fluxes]

    return np.core.records.fromarrays(newLC.transpose(), dtype=headersType)


def writeLCsForZDCF(lc1, lc2, file1, file2):

    np.savetxt(file1, lc1, delimiter=',')
    np.savetxt(file2, lc2, delimiter=',')

    return


def runCrossCorrZDCF(platformNow, outPrefix, nSims, file1, file2, logOutput):

    zdcfCmd = './zdcf{}.out'.format(platformNow)
    userInput = '$"2\n./{}\nn\n0\nn\n{}\n{}\n{}"'.format(outPrefix, nSims, file1, file2)
    subprocess.call('{} <<< {} >& {}'.format(zdcfCmd, userInput, logOutput), shell=True)

    return


def calcLag(platformNow, zdcfFile, lowBound, upperBound, logOutput):

    plikeCmd = './plike{}.out'.format(platformNow)
    userInput = '$"{}\n{}\n{}"'.format(zdcfFile, lowBound, upperBound)
    subprocess.call('{} <<< {} >& {}'.format(plikeCmd, userInput, logOutput), shell=True)

    plikeFile = open(logOutput, 'r')

    lag, lagDn, lagUp = None, None, None

    for line in plikeFile:
        if '1 sigma ML' in line:
            lag = float(line.split()[7])
        if '= (' in line:
            lagDn = float(line.split()[3])
            lagUp = float(line.split()[5].split(')')[0])

    if lag is not None and lagDn is not None and lagUp is not None:
        return lag, lagDn, lagUp
    else:
        print('could not find lag')
        print(lag, lagDn, lagUp)
        sys.exit(1)


def calcPSD(lc, sampling):

    freq, power = periodogram(lc['flux'], 1/(sampling))
    power = power[freq > 0]
    freq = freq[freq > 0]
    # if np.mean(lc['flux']) < 0.01:
    #     np.save('power.npy', power)
    #     np.save('freq.npy', freq)
    # if np.max(power)/np.mean(power) > 30:  # Help the fit to converge
    #     power = np.delete(power, 0)
    #     freq = np.delete(freq, 0)
    pars, varMatrix = curve_fit(psd, freq, power, p0=[np.max(power), 0.5], maxfev=5000)

    return freq, power, pars


def simLC(origLC, n, sampling, beta, roll):

    simFluxes = generate_power_law(n, sampling, beta)
    simLC = buildLC(origLC, simFluxes, sampling, roll)

    return simLC


def getSimCorrFactors(prefix, lc1, lc2, n1, n2, sampling1, sampling2,
                      beta1, beta2, runSim, nSim):

    if runSim:
        # First run the actual LCs to get the correct ZDCF length
        writeLCsForZDCF(lc1, lc2, prefix + 'lc1.txt', prefix + 'lc2.txt')
        runCrossCorrZDCF(platformNow, prefix + 'ccf', 1000,
                         prefix + 'lc1.txt', prefix + 'lc2.txt', prefix + 'outZDCF.log')
        nomZdcf = np.loadtxt(prefix + 'ccf.dcf')

        simLC1 = simLC(lc1, n1, sampling1, beta1, 3)
        simLC2 = simLC(lc2, n2, sampling2, beta2, 3)
        plotLC(simLC1, simLC2, prefix + 'simLC.pdf')

        writeLCsForZDCF(simLC1, simLC2, prefix + 'simLC1.txt', prefix + 'simLC2.txt')
        runCrossCorrZDCF(platformNow, prefix + 'simCCF', 1000,
                         prefix + 'simLC1.txt', prefix + 'simLC2.txt', prefix + 'simZDCF.log')

        lag, lagDn, lagUp = calcLag(platformNow, prefix + 'simCCF.dcf',
                                    -200, 200, prefix + 'plikeSimCCF.log')
        zdcf = np.loadtxt(prefix + 'simCCF.dcf')

        # Make sure to get the correct ZDCF length in the simulation
        i_iter = 1
        while len(nomZdcf) != len(zdcf):
            print('rerunning first simulation (iter-{})'.format(i_iter))
            print('nominal zDCF - {}, sim zDCF - {}'.format(len(nomZdcf), len(zdcf)))
            i_iter += 1
            simLC1 = simLC(lc1, n1, sampling1, beta1, 0)
            simLC2 = simLC(lc2, n2, sampling2, beta2, 0)
            plotLC(simLC1, simLC2, prefix + 'simLC.pdf')

            writeLCsForZDCF(simLC1, simLC2, prefix + 'simLC1.txt', prefix + 'simLC2.txt')
            runCrossCorrZDCF(platformNow, prefix + 'simCCF', 1000,
                             prefix + 'simLC1.txt', prefix + 'simLC2.txt', prefix + 'simZDCF.log')

            lag, lagDn, lagUp = calcLag(platformNow, prefix + 'simCCF.dcf',
                                        -200, 200, prefix + 'plikeSimCCF.log')
            zdcf = np.loadtxt(prefix + 'simCCF.dcf')

        plotZDCF(zdcf, lag, lagDn, lagUp, prefix + 'simZdcf.pdf')

        distZDCF = zdcf[:, 3]  # To calculate quantiles
        lags, corrFactors = zdcf[:, 0], zdcf[:, 3]

        roll = 3
        for i_sim in range(nSim):
            if i_sim + 1 % 10 == 0:
                print('simulated iteration - {}'.format(i_sim + 1))
            simLC1 = simLC(lc1, n1, sampling1, beta1, roll)
            simLC2 = simLC(lc2, n2, sampling2, beta2, roll)

            writeLCsForZDCF(simLC1, simLC2, prefix + 'tempSimLC1.txt', prefix + 'tempSimLC2.txt')
            runCrossCorrZDCF(platformNow, prefix + 'tempSimCCF', 1000,
                             prefix + 'tempSimLC1.txt', prefix + 'tempSimLC2.txt',
                             prefix + 'tempSimZDCF.log')

            zdcf = np.loadtxt(prefix + 'tempSimCCF.dcf')
            # Sometimes the zDCF binning in the simulation is different, skip those cases
            if len(zdcf[:, 3]) != len(distZDCF):
                roll = 0
                continue
            lags = np.append(lags, zdcf[:, 0])
            corrFactors = np.append(corrFactors, zdcf[:, 3])
            distZDCF = np.c_[distZDCF, zdcf[:, 3]]

        np.save(prefix + 'simResults_lags.npy', lags)
        np.save(prefix + 'simResults_corrFactors.npy', corrFactors)
        np.save(prefix + 'simResults_distZDCF.npy', distZDCF)
    else:
        distZDCF = np.load(prefix + 'simResults_distZDCF.npy')
        lags = np.load(prefix + 'simResults_lags.npy')
        corrFactors = np.load(prefix + 'simResults_corrFactors.npy')

    return distZDCF, lags, corrFactors


def getFullRangeSimLC(noFlareLC, flareLC, simNoFlareLC):

    simNoFlareFlux = simNoFlareLC['flux'] - np.min(simNoFlareLC['flux'])
    simNoFlareFlux = (np.mean(noFlareLC['flux'])/np.mean(simNoFlareFlux))*simNoFlareFlux
    simNoFlareLC['flux'] = simNoFlareFlux
    simNoFlareLC['fluxErr'] = np.random.normal(np.mean(noFlareLC['fluxErr']),
                                               np.std(noFlareLC['fluxErr']),
                                               len(noFlareLC['flux']))
    simFlareFlux = np.random.normal(np.mean(simNoFlareFlux),
                                    np.mean(noFlareLC['fluxErr']),
                                    len(flareLC['flux']))
    simFlareFluxErr = np.random.normal(np.mean(noFlareLC['fluxErr']),
                                       np.std(noFlareLC['fluxErr']),
                                       len(flareLC['flux']))
    if np.any(simFlareFlux < 0):
        simFlareFlux = simFlareFlux - np.min(simFlareFlux)
    headersType = {'names': ('mjd', 'flux', 'fluxErr'),
                   'formats': ('f8', 'f8', 'f8')}

    simFlarePeriod = np.c_[flareLC['mjd'], simFlareFlux, simFlareFluxErr]
    simFlarePeriod = np.core.records.fromarrays(simFlarePeriod.transpose(), dtype=headersType)

    simFullRange = np.concatenate((simNoFlareLC, simFlarePeriod), axis=0)
    simFullRange.sort(order='mjd')

    return simFullRange


def addSimLCtoLC(lc, simFullRange, tau):

    opacity = np.exp(-tau)
    secondaryLC = np.copy(lc)
    secondaryLC['flux'] = secondaryLC['flux']*opacity
    secondaryLC['fluxErr'] = secondaryLC['fluxErr']*opacity
    secondaryLC['flux'] = secondaryLC['flux'] + simFullRange['flux']*(1 - opacity)
    secondaryLC['fluxErr'] = secondaryLC['fluxErr'] + simFullRange['fluxErr']*(1 - opacity)

    return secondaryLC


def plotPSD(freq, power, pars, experiment, outputFile):

    pu = pUtil.plotUtil()
    colors = pu.getColors()
    markers = pu.getMarkers()
    lines = pu.getLines()
    fontsize = pu.getFontsize(bigPlot=True)
    markersize = pu.getMarkersize(bigPlot=True)

    fig = plt.figure(figsize=(15, 10))

    xTitle = 'frequency [Hz]'
    yTitle = 'PSD [V**2/Hz]'

    plt.plot(freq, power, markers[0],
             color=colors[0],
             markersize=markersize,
             label='{}'.format(experiment),
             zorder=10)
    plt.plot(freq, psd(freq, *pars), '--',
             color=colors[1],
             label='Fit',
             zorder=10)
    plt.gca().set_xscale("log", nonposx='clip')
    plt.gca().set_yscale("log", nonposy='clip')

    plt.xlabel(xTitle, fontsize=20, labelpad=0)
    plt.ylabel(yTitle, fontsize=20)
    plt.gca().get_yaxis().get_offset_text().set_fontsize(15)

    plt.legend(loc='upper right', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.grid(True)
    plt.gca().set_axisbelow(True)
    plt.tight_layout()

    plt.savefig(outputFile)

    plt.close('all')

    return


def plotLC(lc1, lc2, outputFile):

    pu = pUtil.plotUtil()
    colors = pu.getColors()
    markers = pu.getMarkers()
    lines = pu.getLines()
    fontsize = pu.getFontsize(bigPlot=True)
    markersize = pu.getMarkersize(bigPlot=True)
    elinewidth = pu.getElinewidth(bigPlot=True)

    fig = plt.figure(figsize=(25, 10))
    ax1 = fig.add_subplot(211)

    xTitle = 'MJD'
    yTitle = 'Flux (a.u.)'

    plt.errorbar(lc1['mjd'],
                 lc1['flux'],
                 yerr=lc1['fluxErr'],
                 fmt=markers[0], color=colors[0],
                 markersize=markersize, elinewidth=elinewidth,
                 zorder=4, label=r'VERITAS')
    plt.xlabel(xTitle, fontsize=18, labelpad=0)
    plt.ylabel(yTitle, fontsize=18)
    plt.gca().get_yaxis().get_offset_text().set_fontsize(18)

    plt.legend(loc='upper left', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=18)

    plt.grid(True)
    plt.gca().set_axisbelow(True)
    plt.tight_layout()

    fig.add_subplot(212, sharex=ax1)

    plt.errorbar(lc2['mjd'],
                 lc2['flux'],
                 yerr=lc2['fluxErr'],
                 fmt=markers[1], color=colors[1],
                 markersize=markersize, elinewidth=elinewidth,
                 zorder=4, label=r'Fermi-LAT')

    plt.xlabel(xTitle, fontsize=18, labelpad=0)
    plt.ylabel(yTitle, fontsize=18)
    plt.gca().get_yaxis().get_offset_text().set_fontsize(18)
    plt.legend(loc='upper left', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=18)

    plt.grid(True)
    plt.gca().set_axisbelow(True)
    plt.tight_layout()

    ax2 = ax1.twiny()
    mjds = Time(ax1.get_xticks(), format='mjd', scale='utc', out_subfmt='date')
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xticklabels(mjds.iso, fontdict={'fontsize': 18})
    ax2.set_xlim(ax1.get_xlim())

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)

    plt.savefig(outputFile)

    plt.close('all')

    return


def plotZDCF(zdcf, lag, lagDn, lagUp, outputFile, title='', titleLCs=['', ''],
             legendLoc='upper right', textLoc=[-99, -99],
             simDistZDCF=list(), simLags=list(), simCorrFactors=list()):

    pu = pUtil.plotUtil()
    colors = pu.getColors()
    markers = pu.getMarkers()
    lines = pu.getLines()
    fontsize = pu.getFontsize(bigPlot=True) + 2
    markersize = pu.getMarkersize(bigPlot=True)
    elinewidth = pu.getElinewidth(bigPlot=True) + 1

    fig = plt.figure(figsize=(20, 8))

    xTitle = 'Lag (days)'
    yTitle = 'zDCF'

    plt.errorbar(zdcf[:, 0],
                 zdcf[:, 3],
                 xerr=[zdcf[:, 1], zdcf[:, 2]],
                 yerr=[zdcf[:, 4], zdcf[:, 5]],
                 fmt=markers[0], color=colors[0],
                 markersize=markersize, elinewidth=elinewidth,
                 zorder=4, label='t({}) - t({})'.format(titleLCs[1], titleLCs[0]))

    if len(simDistZDCF) > 0 and len(simLags) > 0 and len(simCorrFactors) > 0:

        oneSigma = 0.317/2
        twoSigma = 0.045/2
        threeSigma = 0.0027/2
        fourSigma = 0.000063/2

        quantiles = np.quantile(distZDCF, [fourSigma, 0.999936 + fourSigma,
                                           threeSigma, 0.997 + threeSigma,
                                           twoSigma, 0.954 + twoSigma,
                                           oneSigma, 0.682 + oneSigma], 1)

        plt.hist2d(lags, corrFactors, bins=(len(zdcf[:, 0]), np.linspace(-1, 1, 100)),
                   cmap='Blues', normed=True)
        sigmaTitles = [r'$4\sigma$', r'$3\sigma$', r'$2\sigma$', r'$1\sigma$']
        sigmaLines = [1, 2, 5, 3]
        sigmaColors = pu.getColors('greens')
        for i_quant in range(quantiles.shape[0]):
            i_att = int(i_quant/2.)
            plt.plot(zdcf[:, 0], quantiles[i_quant], lw=5,
                     color=sigmaColors[i_att], linestyle=lines[sigmaLines[i_att]], alpha=1)
        for i_att in range(len(sigmaTitles)):
            i_quant = i_att*2 + 1
            quantForLabel = quantiles[i_quant][(np.abs(zdcf[:, 0] - (-850))).argmin()]
            plt.text(-850, quantForLabel,
                     sigmaTitles[i_att],
                     horizontalalignment='left',
                     verticalalignment='bottom',
                     fontsize=30,
                     color=sigmaColors[i_att])

        if np.any(zdcf[:, 3] < quantiles[2]) or np.any(zdcf[:, 3] > quantiles[3]):
            ymin, ymax = plt.gca().get_ylim()
            plt.plot([lag, lag], [ymin, ymax],
                     label=r'Lag (${0:1.0f}^{{+{1:1.0f}}}_{{-{2:1.0f}}}$)'.format(lag,
                                                                                  lagUp - lag,
                                                                                  lag - lagDn),
                     lw=4, linestyle=lines[1], color=colors[0])
            plt.fill_between([lagDn, lagUp], ymin, ymax, facecolor=colors[0], alpha=0.2)

    plt.gca().set_xlim([-1000, 1000])

    plt.title(title, fontsize=30, y=1.03)
    if textLoc[0] >= 0:
        plt.text(textLoc[0], textLoc[1],
                 'ICRC 2019', transform=plt.gca().transAxes,
                 horizontalalignment='left',
                 verticalalignment='bottom',
                 fontsize=30)

    plt.xlabel(xTitle, fontsize=fontsize, labelpad=0)
    plt.ylabel(yTitle, fontsize=fontsize)
    plt.gca().get_yaxis().get_offset_text().set_fontsize(fontsize)

    plt.legend(loc=legendLoc, fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.grid(True)
    plt.gca().set_axisbelow(True)
    plt.tight_layout()

    plt.savefig(outputFile)

    plt.close('all')


if __name__ == '__main__':

    start = time.time()
    logStdout = colourLogger.initStdOutLogger()

    platformNow = platform.system()
    if platformNow == 'Darwin':
        platformNow = 'Mac'

    config = gUtil.readYamlFile(logStdout, 'config.yaml')
    for sourceShort, sourceLong in config['sources'].items():
        for configNameNow in config['config'][sourceShort]:

            configNow = config['config'][sourceShort][configNameNow]

            logStdout.info([['wb', 'Running'],
                            ['p', ' {}'.format(sourceLong)],
                            ['wb', 'for'],
                            ['b', '({}, {})'.format(configNow['titles'][0],
                                                    configNow['titles'][1])]])

            prefix = '{}_{}_{}_{}/'.format(sourceShort, configNow['names'][0],
                                           configNow['names'][1], configNow['binning'])
            title = '{} ({}, {}) correlation'.format(sourceLong, configNow['titles'][0],
                                                     configNow['titles'][1])

            lc1, lc2 = readLCs(prefix + '{}_{}_{}.txt'.format(sourceShort, configNow['binning'],
                                                              configNow['names'][0].lower()),
                               prefix + '{}_{}_{}.txt'.format(sourceShort, configNow['binning'],
                                                              configNow['names'][1].lower()))

            plotLC(lc1, lc2, prefix + 'lightcurve.pdf')

            sInDay = 60*60*24
            sampling1 = float(configNow['sampling'][0])*sInDay
            freq, power, pars = calcPSD(lc1, sampling1)
            beta1 = pars[1]
            plotPSD(freq, power, pars, configNow['titles'][0],
                    prefix + 'psd{}.pdf'.format(configNow['names'][0]))

            sampling2 = float(configNow['sampling'][1])*sInDay
            freq, power, pars = calcPSD(lc2, sampling2)
            beta2 = pars[1]
            plotPSD(freq, power, pars, configNow['titles'][1],
                    prefix + 'psd{}.pdf'.format(configNow['names'][1]))

            np.random.seed(4242)
            distZDCF, lags, corrFactors = getSimCorrFactors(prefix, lc1, lc2,
                                                            len(lc1['flux']), len(lc2['flux']),
                                                            sampling1, sampling2,
                                                            beta1, beta2,
                                                            runSim=configNow['runSim'],
                                                            nSim=int(float(configNow['nSim'])))

            writeLCsForZDCF(lc1, lc2, prefix + 'lc1.txt', prefix + 'lc2.txt')
            runCrossCorrZDCF(platformNow, prefix + 'ccf', 1000,
                             prefix + 'lc1.txt', prefix + 'lc2.txt', prefix + 'outZDCF.log')

            lag, lagDn, lagUp = calcLag(platformNow, prefix + 'ccf.dcf',
                                        float(configNow['lagRange'][0]),
                                        float(configNow['lagRange'][1]),
                                        prefix + 'plikeCCF.log')
            zdcf = np.loadtxt(prefix + 'ccf.dcf')
            plotZDCF(zdcf, lag, lagDn, lagUp, prefix + 'zdcf.pdf',
                     title, configNow['titles'], configNow['legend'],
                     [configNow['xInfo'], configNow['yInfo']],
                     distZDCF, lags, corrFactors)

            if sourceShort == '1ES1011' and 'tau' in configNow['names'][0]:

                headersType = {'names': ('mjd', 'flux', 'fluxErr'),
                               'formats': ('f8', 'f8', 'f8')}

                fullLC = np.loadtxt('1ES1011_tau-1_tau-3_nightly/1ES1011_nightly_tau-1.txt',
                                    delimiter=',', dtype=headersType)

                tau = int(configNow['names'][1].split('tau-')[1]) - 1
                noFlareLC = fullLC[(fullLC['mjd'] < 56688) | (fullLC['mjd'] > 56715)]
                flareLC = fullLC[(fullLC['mjd'] > 56688) & (fullLC['mjd'] < 56715)]
                freq, power, pars = calcPSD(noFlareLC, sampling2)
                simNoFlareLC = simLC(noFlareLC, len(noFlareLC['mjd']), sampling2, pars[1], 0)
                simFullLC = getFullRangeSimLC(noFlareLC, flareLC, simNoFlareLC)
                plotLC(fullLC, simFullLC, prefix + 'simSecondaryLC.pdf')
                secondaryLC = addSimLCtoLC(fullLC, simFullLC, tau)
                plotLC(fullLC, secondaryLC, prefix + 'fullAndSecondaryLC.pdf')
                plotLC(lc1, secondaryLC, prefix + 'lowBinAndSecondaryLC.pdf')
                plotLC(lc2, secondaryLC, prefix + 'tauAndSecondaryLC.pdf')

                writeLCsForZDCF(lc1, secondaryLC, prefix + 'lc1ForSecondary.txt',
                                prefix + 'secondaryLC.txt')
                runCrossCorrZDCF(platformNow, prefix + 'ccfSecondary', 1000,
                                 prefix + 'lc1ForSecondary.txt', prefix + 'secondaryLC.txt',
                                 prefix + 'secondaryZDCF.log')

                lag, lagDn, lagUp = calcLag(platformNow, prefix + 'ccfSecondary.dcf', -200, 200,
                                            prefix + 'plikeSecondaryCCF.log')
                zdcf = np.loadtxt(prefix + 'ccfSecondary.dcf')
                title = '{} ({}, simulated {}) correlation'.format(sourceLong,
                                                                   configNow['titles'][0],
                                                                   configNow['titles'][1])
                plotZDCF(zdcf, lag, lagDn, lagUp, prefix + 'secondaryZDCF.pdf',
                         title, configNow['titles'], configNow['legend'],
                         [configNow['xInfo'], configNow['yInfo']],
                         distZDCF, lags, corrFactors)

    end = time.time()
    timeElapsed = str(datetime.timedelta(seconds=end - start))
    print('Run time: {}'.format(timeElapsed))
