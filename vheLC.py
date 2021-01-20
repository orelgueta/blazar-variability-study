#!/usr/bin/python

import sys
import numpy as np
import matplotlib as mlp
from matplotlib import pyplot as plt
from astropy.time import Time
import re
import glob
from itertools import combinations
import argparse
from astropy.stats import bayesian_blocks
from collections import defaultdict
from collections import OrderedDict
from lib import colourLogger
from lib import generalUtil as gUtil
from lib import lcUtil
from lib import dcfUtil
from lib import plotUtil as pUtil


def weightedAvgStd(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))


def weightedCov(x, y, w):
    """
    Return the weighted covariance.
    The weights are per observation, where an observation is a pair (x,y)

    x, y, w -- Numpy ndarrays with the same shape.
    """

    ave_x, std_x = weightedAvgStd(x, w)
    ave_y, std_y = weightedAvgStd(y, w)

    xVec = w*(x - ave_x)  # Easy way to include the w in the numerator
    yVec = (y - ave_y)
    numerator = xVec.T.dot(yVec)
    denominator = np.sum(w)

    return numerator/denominator


def getEnergyThresholds(fileName):

    baseFileName = fileName.split('/')[-1]
    energyThresholdRe = re.compile(r'_(\d+)')
    energyThreshold = [int(energyNow) for energyNow in
                       energyThresholdRe.findall(baseFileName)]
    energyThreshold = [200 if energyNow == 195 else energyNow for energyNow in energyThreshold]

    return energyThreshold


def getTauBin(fileName):

    baseFileName = fileName.split('/')[-1]
    tauRe = re.compile(r'_tau-(\d+)')
    tauBin = [int(tauNow) for tauNow in tauRe.findall(baseFileName)][0]

    return tauBin


def getTauBinLatex(tauBin):

    if tauBin == 1:
        return r'$\tau < 1$'
    if tauBin == 2:
        return r'$1 < \tau < 2$'
    if tauBin == 3:
        return r'$\tau$ > 2'


def getRhoLatex(tauBin1, tauBin2):

    rhoString = r'$\rho$'
    tauString = '({}, {})'.format(getTauBinLatex(tauBin1), getTauBinLatex(tauBin2))

    return rhoString + tauString


def getDCFbins(tauBin1, tauBin2):

    return 'DCF({}, {})'.format(getTauBinLatex(tauBin1), getTauBinLatex(tauBin2))


def getCorrelationFactors(veritasLCs, weighted):

    correlationFactors = dict()
    for comb in combinations([1, 2, 3], 2):
        rhoLatex = getRhoLatex(comb[0], comb[1])
        correlationFactors[rhoLatex] = calcCorr(veritasLCs[comb[0]], veritasLCs[comb[1]], weighted)

    return correlationFactors


def calcCorr(lc1, lc2, weighted):

    # The following adds explicit zero-flux points in place of missing measurements
    notIn2 = np.isin(lc1['DateMJD'], lc2['DateMJD'], invert=True)
    for dateNotIn2 in lc1[notIn2]['DateMJD']:
        arrayToInsert = [[dateNotIn2, 0, np.max(lc1['Flux Error'])]]
        lc2 = np.append(lc2, arrayToInsert, axis=0)

    notIn1 = np.isin(lc2['DateMJD'], lc1['DateMJD'], invert=True)
    for dateNotIn1 in lc2[notIn1]['DateMJD']:
        arrayToInsert = [[dateNotIn1, 0, np.max(lc2['Flux Error'])]]
        lc1 = np.append(lc1, arrayToInsert, axis=0)

    if weighted:
        combinedError = np.sqrt(np.power(lc1['Flux Error'], 2) + np.power(lc2['Flux Error'], 2))
        weights = np.min(combinedError)/combinedError

        xxCov = weightedCov(lc1['Flux'], lc1['Flux'], weights)
        yyCov = weightedCov(lc2['Flux'], lc2['Flux'], weights)
        xyCov = weightedCov(lc1['Flux'], lc2['Flux'], weights)

        corrFactor = xyCov/np.sqrt(xxCov*yyCov)

        return corrFactor
    else:
        corrFactor = np.corrcoef(np.array([lc1['Flux'], lc2['Flux']]))
        return corrFactor[0][1]


def dcfForSource(sourceNow, binning, minPeriod, maxPeriod, binWidth):

    dcf = OrderedDict()
    prefixLC = 'vheLC/dcf/lightcurves'
    fileList = glob.glob('{}/{}_{}_tau*.txt'.format(prefixLC, sourceNow, binning))
    fileList.sort()
    tauRe = re.compile(r'_tau-(\d+)')
    for comb in combinations(fileList, 2):
        tauBins = list()
        tauBins.append([int(tauNow) for tauNow in tauRe.findall(comb[0])][0])
        tauBins.append([int(tauNow) for tauNow in tauRe.findall(comb[1])][0])
        tauBins.sort()
        dcfBin = getDCFbins(tauBins[0], tauBins[1])
        dcfOutputPrefix = '{}_between_tau-{}_tau-{}_{}'.format(sourceNow, tauBins[0],
                                                               tauBins[1], binning)

        dcf[dcfBin] = dcfUtil.dcfTwoLCs('vheLC', dcfOutputPrefix,
                                        comb[0], comb[1],
                                        minPeriod, maxPeriod, binWidth)

    return dcf


def plotDCF(dcf, sourceNow, binning, binWidth, configAnal):

    sourceNowShort = list(configAnal['sources'].keys())[
                        list(configAnal['sources'].values()).index(sourceNow)]

    pu = pUtil.plotUtil()
    colors = pu.getColors('autumn')
    markers = pu.getMarkers()
    lines = pu.getLines()
    fontsize = pu.getFontsize(bigPlot=True)
    markersize = pu.getMarkersize(bigPlot=True)
    elinewidth = pu.getElinewidth(bigPlot=True)

    fig = plt.figure(figsize=(25, 15))

    for i_plt, (bins, dcfNow) in enumerate(dcf.items()):

        if i_plt == 0:
            ax1 = fig.add_subplot(311)
        else:
            plotNow = 311 + i_plt
            fig.add_subplot(plotNow, sharex=ax1)

        xTitle = 'Lag (days)'
        yTitle = 'Correlation coefficient'
        plt.errorbar(dcfNow[:, 0],
                     dcfNow[:, 1],
                     xerr=binWidth/2.,
                     yerr=dcfNow[:, 2],
                     fmt=markers[i_plt], color=colors[i_plt],
                     markersize=markersize, elinewidth=elinewidth,
                     label=bins,
                     zorder=10)

        plt.xlabel(xTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'], labelpad=0)
        plt.ylabel(yTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'])
        plt.gca().get_yaxis().get_offset_text().set_fontsize(configAnal['plot']
                                                                       ['mwlLC']
                                                                       ['axesLabelFont'])

        if i_plt == 0:
            plt.title('{}, {} bins'.format(sourceNow, binning), fontsize=30, y=1.03)

        plt.legend(loc='upper right', fontsize=configAnal['plot']['mwlLC']['legendFont'])
        plt.tick_params(axis='both', which='major',
                        labelsize=configAnal['plot']['mwlLC']['axesLabelFont'])
        plt.grid(True)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()

    fig.subplots_adjust(hspace=0)
    plt.savefig('vheLC/dcf/plots/{0}_{1}_dcf.pdf'.format(sourceNowShort, binning))

    plt.close('all')


def plotAllLC(logStdout, veritasLCs, sourceNow, binningNow,
              energyThresholds, correlationFactors, configAnal, corrTable):

    sourceNowShort = list(configAnal['sources'].keys())[
                        list(configAnal['sources'].values()).index(sourceNow)]

    pu = pUtil.plotUtil()
    colors = pu.getColors('autumn')
    markers = pu.getMarkers()
    lines = pu.getLines()
    fontsize = pu.getFontsize(bigPlot=True)
    markersize = pu.getMarkersize(bigPlot=True)
    elinewidth = pu.getElinewidth(bigPlot=True)

    fig = plt.figure(figsize=(25, 15))

    for i_plt, (tauBin, veritasFluxes) in enumerate(veritasLCs.items()):

        if i_plt == 0:
            ax1 = fig.add_subplot(311)
        else:
            plotNow = 311 + i_plt
            fig.add_subplot(plotNow, sharex=ax1)

        plotInfo = '{} TeV < E < {} TeV ({})'.format(energyThresholds[tauBin][0]/1000.,
                                                     energyThresholds[tauBin][1]/1000.,
                                                     getTauBinLatex(tauBin))

        priorFile = '{}_{}_{}_result.txt'.format('veritas', binningNow,
                                                 configAnal['nBayesSims'][sourceNowShort])
        gamma = lcUtil.readGammaFromFile(logStdout, 'vheLC', priorFile,
                                         sourceNowShort, configAnal['bayesChange'])
        logStdout.info([['p', sourceNow],
                        ['wb', '{} {} gamma -'.format('veritas', binningNow)],
                        ['b', gamma]])

        veritasBlocks = bayesian_blocks(veritasFluxes['DateMJD'], veritasFluxes['Flux'],
                                        veritasFluxes['Flux Error'], fitness='measures',
                                        gamma=gamma)

        veritasFluxMean, veritasFluxErrMean = lcUtil.calcBlockFlux(veritasFluxes['DateMJD'],
                                                                   veritasFluxes['Flux'],
                                                                   veritasFluxes['Flux Error'],
                                                                   veritasBlocks)

        # Correct for wrong start possitons of edges.
        # Implemented output of bayesian blocks is in the biddle between two observations.
        # However no statement can be made for a change between two measurments
        # so that a new block should only begin at the start of a observation.
        correctedEdges = np.array([])
        for edge in veritasBlocks:
            aboveEdge = veritasFluxes['DateMJD'] >= edge
            correctedEdges = np.append(correctedEdges,
                                       min(veritasFluxes['DateMJD'][aboveEdge] -
                                           veritasFluxes['Date Error'][aboveEdge]/2))

        veritasBlocks = correctedEdges

        xTitle = 'MJD'
        yTitle = r'$\mathcal{{F}} (\times 10^{11})$ [cm$^{-2} s^{-1}$]'
        fVeritas = 1e11

        plt.errorbar(veritasFluxes['DateMJD'],
                     fVeritas*veritasFluxes['Flux'],
                     xerr=veritasFluxes['Date Error'],
                     yerr=fVeritas*veritasFluxes['Flux Error'],
                     fmt=markers[i_plt], color=colors[i_plt], mfc='none',
                     markersize=markersize, elinewidth=elinewidth,
                     label=plotInfo,
                     zorder=10)

        # Plot Bayesian blocks
        plt.gca().step(x=veritasBlocks, y=fVeritas*veritasFluxMean, zorder=3, linewidth=3,
                       where='post', color='darkgrey', alpha=0.8)
        for i in range(len(veritasBlocks)-1):
            x1 = [veritasBlocks[i], veritasBlocks[i+1], veritasBlocks[i+1], veritasBlocks[i]]
            y1 = [fVeritas*(veritasFluxMean[i] + veritasFluxErrMean[i]),
                  fVeritas*(veritasFluxMean[i] + veritasFluxErrMean[i]),
                  fVeritas*(veritasFluxMean[i] - veritasFluxErrMean[i]),
                  fVeritas*(veritasFluxMean[i] - veritasFluxErrMean[i])]
            label = None
            if i == 0:
                label = r'Bayesian blocks'
            plt.gca().fill(x1, y1, color='darkgrey', label=label, alpha=0.4, zorder=2)

        # Get the relevant correlation factor
        corrToPrint = list()
        for comb in combinations([1, 2, 3], 2):
            if tauBin in comb:
                tauNowIndex = comb.index(tauBin)
                otherBin = 0 if tauNowIndex == 1 else 1
                # Keep two versions here, one to print on the plot
                # with the "nice" order and one with the key order
                corrToPrint.append([getRhoLatex(comb[0], comb[1]),
                                   getRhoLatex(comb[tauNowIndex], comb[otherBin])])

        correlations = '{} = {:1.2f}\n'.format(corrToPrint[0][1],
                                               correlationFactors[corrToPrint[0][0]])
        correlations += '{} = {:1.2f}'.format(corrToPrint[1][1],
                                              correlationFactors[corrToPrint[1][0]])

        plt.xlabel(xTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'], labelpad=0)
        plt.ylabel(yTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'])
        plt.gca().get_yaxis().get_offset_text().set_fontsize(configAnal['plot']
                                                                       ['mwlLC']
                                                                       ['axesLabelFont'])

        props = dict(boxstyle='round', facecolor='none')
        plt.text(configAnal['plot']['vheLC'][sourceNowShort]['xInfo'],
                 configAnal['plot']['vheLC'][sourceNowShort]['yInfo'],
                 correlations,
                 horizontalalignment='left',
                 verticalalignment='top',
                 fontsize=configAnal['plot']['mwlLC']['textFont'],
                 transform=plt.gca().transAxes,
                 bbox=props,
                 zorder=7)

        if i_plt == 0:
            plt.text(configAnal['plot']['vheLC'][sourceNowShort]['xICRC'],
                     configAnal['plot']['vheLC'][sourceNowShort]['yICRC'],
                     'ICRC 2019',
                     horizontalalignment='left',
                     verticalalignment='center',
                     fontsize=configAnal['plot']['mwlLC']['textFont'],
                     transform=plt.gca().transAxes)

        plt.legend(bbox_to_anchor=(configAnal['plot']['vheLC'][sourceNowShort]['xLegend'],
                                   configAnal['plot']['vheLC'][sourceNowShort]['yLegend']),
                   loc='center left', fontsize=configAnal['plot']['mwlLC']['legendFont'])

        if i_plt == 0:
            plt.title('{}, {} bins'.format(sourceNow, binningNow),
                      fontsize=30, y=1.18)

        ax1.set_xlim([56650, 56725])
        plt.tick_params(axis='both', which='major',
                        labelsize=configAnal['plot']['mwlLC']['axesLabelFont'])
        plt.grid(True)
        plt.gca().set_axisbelow(True)

    ax2 = ax1.twiny()
    mjds = Time(ax1.get_xticks(), format='mjd', scale='utc', out_subfmt='date')
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xticklabels(mjds.iso,
                        fontdict={'fontsize': configAnal['plot']['mwlLC']['axesLabelFont']})
    ax2.set_xlim(ax1.get_xlim())

    plt.tight_layout()
    fig.subplots_adjust(left=0.065, top=0.895, hspace=0)
    plt.savefig('vheLC/{0}/{0}_{1}_vheLCs.pdf'.format(sourceNowShort, binningNow))

    plt.close('all')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Produce VHE lightcurves and compute DCF')
    args = parser.parse_args()

    logStdout = colourLogger.initStdOutLogger()

    configAnalFile = 'configAnalysis.yaml'
    configAnal = gUtil.readYamlFile(logStdout, configAnalFile)

    # FIXME - make this local and move to configAnalysis.yaml
    corrTableFile = ('/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/'
                     'crabStability/plotLC/correctionFactors.txt')
    corrTable = lcUtil.readCorrTable(corrTableFile)

    correlationFactors = dict()
    for sourceNowShort, sourceNow in configAnal['sources'].items():
        sourceDir = 'makeLC/{}/'.format(sourceNowShort)

        logStdout.info([['wb', 'Running'],
                        ['p', sourceNow]])

        for binningNow in configAnal['binning']:
            veritasLCs = dict()
            energyThresholds = dict()
            veritasFiles = glob.glob(sourceDir + '*{}LightCurve_tau*.txt'.format(binningNow))
            veritasFiles.sort()
            for fileNow in veritasFiles:

                tauBin = getTauBin(fileNow)
                energyThresholds[tauBin] = getEnergyThresholds(fileNow)

                veritasLCs[tauBin] = lcUtil.readVeritasLC(logStdout, fileNow)
                veritasLCs[tauBin] = lcUtil.correctFluxesFromCrabLC(veritasLCs[tauBin],
                                                                    corrTable)

                # Write lightcurves for DCF
                if binningNow in configAnal['dcf'][sourceNowShort]['binning']:
                    dcfFilePrefix = '{}_{}_tau-{}'.format(sourceNowShort, binningNow, tauBin)
                    dcfUtil.saveDataForDCF('vheLC', dcfFilePrefix,
                                           veritasLCs[tauBin]['DateMJD'],
                                           veritasLCs[tauBin]['Flux'],
                                           veritasLCs[tauBin]['Flux Error'])

            correlationFactors = getCorrelationFactors(veritasLCs,
                                                       configAnal['correlation']['weighted'])
            logStdout.info([['wb', 'Plotting LCs for'],
                            ['p', sourceNow],
                            ['wb', 'with'],
                            ['bb', binningNow],
                            ['wb', 'binning']])
            plotAllLC(logStdout, veritasLCs, sourceNow, binningNow,
                      energyThresholds, correlationFactors, configAnal, corrTable)

            dcf = dcfForSource(sourceNowShort, binningNow,
                               configAnal['dcf'][sourceNowShort]['minPeriod'],
                               configAnal['dcf'][sourceNowShort]['maxPeriod'],
                               configAnal['dcf'][sourceNowShort]['binWidth'])

            logStdout.info([['wb', 'Plotting DCF for'],
                            ['p', sourceNow],
                            ['wb', 'with'],
                            ['bb', binningNow],
                            ['wb', 'binning']])
            plotDCF(dcf, sourceNow, binningNow,
                    configAnal['dcf'][sourceNowShort]['binWidth'], configAnal)

    # Merge lightcurve plots into one PDF
    fileToMerge = list()
    prefix = 'vheLC'
    for binning in configAnal['binning'].keys():
        for sourceNowShort in configAnal['sources'].keys():
            fileToMerge.append('{0}/{1}/{1}_{2}_vheLCs.pdf'.format(prefix,
                                                                   sourceNowShort,
                                                                   binning))
        gUtil.mergePDFs(logStdout, fileToMerge, '{}/{}_vheLCs.pdf'.format(prefix, binning))
        fileToMerge.clear()

    # Merge DCF plots into one PDF
    fileToMerge = defaultdict(list)
    for sourceNowShort in configAnal['sources'].keys():
        for binningNow in configAnal['dcf'][sourceNowShort]['binning']:
            fileToMerge[binningNow].append('{}/dcf/plots/{}_{}_dcf.pdf'.format(prefix,
                                                                               sourceNowShort,
                                                                               binningNow))
    for binningNow in fileToMerge.keys():
        gUtil.mergePDFs(logStdout, fileToMerge[binningNow],
                        '{}/{}_dcf.pdf'.format(prefix, binningNow))
        fileToMerge[binningNow].clear()
