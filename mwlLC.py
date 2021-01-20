#!/usr/bin/python

import os.path
import os
import sys
import glob
import subprocess
import numpy as np
import matplotlib as mlp
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import re
from astropy.stats import bayesian_blocks
from astropy.time import Time
from collections import defaultdict
from collections import OrderedDict
import argparse
from lib import colourLogger
from lib import generalUtil as gUtil
from lib import lcUtil
from lib import dcfUtil
from lib import plotUtil as pUtil


class mwl(object):
    """
    A class to study and plot MWL data
    """

    def __init__(self, *args, **kwargs):

        self.logStdout = colourLogger.initStdOutLogger()
        pu = pUtil.plotUtil()
        self.colors = pu.getColors()
        self.markers = pu.getMarkers()
        self.lines = pu.getLines()
        self.fontsize = pu.getFontsize(bigPlot=True)
        self.markersize = pu.getMarkersize(bigPlot=True)
        self.elinewidth = pu.getElinewidth(bigPlot=True)
        self.capsize = pu.getCapsize(bigPlot=True)

    def getLogger(self):
        return self.logStdout

    def compareCorrectedLC(self, veritasFluxes, corrTable, sourceNow, veritasBinning, configAnal):

        sourceNowShort = sourceNow.split('+')[0].replace(' ', '')

        origLC = veritasFluxes
        corrLC = lcUtil.correctFluxesFromCrabLC(origLC, corrTable)

        fig = plt.figure(figsize=(15, 7))

        # FIXME Why in ED the unit is /cm^2/s and here it has ergs as well?
        xTitle = 'MJD'
        yTitle = r'Flux ($\times 10^{10}$) [cm$^{-2} s^{-1}$]'
        fVeritas = 1e10
        xText = 0.285
        plt.errorbar(origLC['DateMJD'],
                     fVeritas*origLC['Flux'],
                     xerr=origLC['Date Error'],
                     yerr=fVeritas*origLC['Flux Error'],
                     markerfacecolor='none', fmt=self.markers[0], markersize=self.markersize,
                     color=self.colors[0], elinewidth=self.elinewidth, capsize=self.capsize,
                     label='Before correction')
        plt.errorbar(corrLC['DateMJD'],
                     fVeritas*corrLC['Flux'],
                     xerr=corrLC['Date Error'],
                     yerr=fVeritas*corrLC['Flux Error'],
                     markerfacecolor='none', fmt=self.markers[1], markersize=self.markersize,
                     color=self.colors[1], elinewidth=self.elinewidth, capsize=self.capsize,
                     label='After correction')
        plt.xlabel(xTitle, fontsize=20, labelpad=0)
        plt.ylabel(yTitle, fontsize=20)
        plt.gca().get_yaxis().get_offset_text().set_fontsize(15)

        energyThreshold = [1000*configAnal['spectralParameters'][sourceNowShort][i] for i in [0, 3]]
        energyThreshold = [200 if energyNow == 195 else energyNow for energyNow in energyThreshold]

        if len(energyThreshold) > 0:
            plotInfo = ('{}, VERITAS {}\n'
                        '{} TeV < E < {} TeV').format(sourceNow,
                                                      veritasBinning + ' bins',
                                                      energyThreshold[0]/1000.,
                                                      energyThreshold[1]/1000.)

        plt.text(xText, 0.91, plotInfo,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=24,
                 transform=plt.gca().transAxes)

        plt.legend(loc='upper left', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(True)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()

        fig.subplots_adjust(hspace=0)

        outputFile = ('mwlLightCurves/compareCorrectedLC/'
                      '{}_{}_compareCorrectedLC.pdf').format(veritasBinning, sourceNowShort)
        plt.savefig(outputFile)

        plt.close('all')

        return

    def dcfForSource(self, sourceNow, binning, minPeriod, maxPeriod, binWidth):

        dcf = OrderedDict()
        prefixLC = 'mwlLightCurves/dcf/lightcurves'
        veritasFile = glob.glob('{}/{}_{}_veritas*.txt'.format(prefixLC, sourceNow, binning))[0]
        fermiFile = glob.glob('{}/{}_{}_fermi*.txt'.format(prefixLC, sourceNow, binning))[0]
        swiftFile = glob.glob('{}/{}_{}_swift*.txt'.format(prefixLC, sourceNow, binning))[0]
        dcfOutputPrefix = '{}_between_{}_{}_{}'.format(sourceNow, 'veritas', 'fermi', binning)
        dcf['VERITAS--Fermi-LAT'] = dcfUtil.dcfTwoLCs('mwlLightCurves', dcfOutputPrefix,
                                                      veritasFile, fermiFile,
                                                      minPeriod, maxPeriod, binWidth)
        dcfOutputPrefix = '{}_between_{}_{}_{}'.format(sourceNow, 'veritas', 'swift', binning)
        dcf['VERITAS--Swift-XRT'] = dcfUtil.dcfTwoLCs('mwlLightCurves', dcfOutputPrefix,
                                                      veritasFile, swiftFile,
                                                      minPeriod, maxPeriod, binWidth)
        dcfOutputPrefix = '{}_between_{}_{}_{}'.format(sourceNow, 'fermi', 'swift', binning)
        dcf['Fermi-LAT--Swift-XRT'] = dcfUtil.dcfTwoLCs('mwlLightCurves', dcfOutputPrefix,
                                                        fermiFile, swiftFile,
                                                        minPeriod, maxPeriod, binWidth)

        return dcf

    def plotDCF(self, dcf, sourceNow, binning, configAnal):

        fig = plt.figure(figsize=(25, 15))

        for i_plt, (experiments, dcfNow) in enumerate(dcf.items()):

            if i_plt == 0:
                ax1 = fig.add_subplot(311)
            else:
                plotNow = 311 + i_plt
                fig.add_subplot(plotNow, sharex=ax1)

            xTitle = 'Lag (days)'
            yTitle = 'Correlation coefficient'
            plt.errorbar(dcfNow[:, 0],
                         dcfNow[:, 1],
                         xerr=configAnal['dcf'][sourceNow]['binWidth']/2.,
                         yerr=dcfNow[:, 2],
                         fmt=self.markers[i_plt], color=self.colors[i_plt],
                         markersize=self.markersize, elinewidth=self.elinewidth,
                         label=experiments,
                         zorder=10)

            plt.xlabel(xTitle, fontsize=20, labelpad=0)
            plt.ylabel(yTitle, fontsize=20)
            plt.gca().get_yaxis().get_offset_text().set_fontsize(15)

            if i_plt == 0:
                plt.title('{} ({} binning)'.format(configAnal['sources'][sourceNow], binning),
                          fontsize=30, y=1.03)

            plt.legend(loc='upper right', fontsize=20)
            plt.tick_params(axis='both', which='major', labelsize=15)
            plt.grid(True)
            plt.gca().set_axisbelow(True)
            plt.tight_layout()

        fig.subplots_adjust(hspace=0)
        plt.savefig('mwlLightCurves/dcf/plots/{}_{}_dcf.pdf'.format(sourceNow, binning))

        plt.close('all')

    def plotLC(self, veritasFluxes, fermiLC, swiftData,
               sourceNow, veritasBinning, fermiBinning, nSims, configAnal):

        sourceNowShort = list(configAnal['sources'].keys())[
                            list(configAnal['sources'].values()).index(sourceNow)]

        priorFile = '{}_{}_{}_result.txt'.format('veritas', veritasBinning, nSims)
        gamma = lcUtil.readGammaFromFile(self.logStdout, 'mwlLightCurves', priorFile,
                                         sourceNowShort, configAnal['bayesChange'])
        self.logStdout.info([['p', sourceNow],
                            ['wb', ' {} {} gamma -'.format('veritas', veritasBinning)],
                            ['b', gamma]])

        veritasBlocks = bayesian_blocks(veritasFluxes['DateMJD'], veritasFluxes['Flux'],
                                        veritasFluxes['Flux Error'], fitness='measures',
                                        gamma=gamma)

        (veritasFluxMean, veritasFluxErrMean) = lcUtil.calcBlockFlux(veritasFluxes['DateMJD'],
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

        energyThreshold = [1000*configAnal['spectralParameters'][sourceNowShort][i] for i in [0, 3]]
        energyThreshold = [200 if energyNow == 195 else energyNow for energyNow in energyThreshold]

        if len(energyThreshold) > 0:
            plotInfo = ('{}, VERITAS {}\n'
                        '{} TeV < E < {} TeV').format(sourceNow,
                                                      veritasBinning + ' bins',
                                                      energyThreshold[0]/1000.,
                                                      energyThreshold[1]/1000.)

        fermiBinningForPlotInfo = fermiBinning
        nWeeksFermiBinning = int(np.mean(np.diff(fermiLC['tmax_mjd']))/7.)
        if nWeeksFermiBinning != 4 and (fermiBinning == 'monthly' or fermiBinning == 'weekly'):
            fermiBinningForPlotInfo = '{}-week'.format(nWeeksFermiBinning)
        if nWeeksFermiBinning == 0:
            fermiBinningForPlotInfo = '{}-day'.format(int(np.mean(np.diff(fermiLC['tmax_mjd']))))

        tsUpperLimit = 0.001
        fermiUpperLimits = fermiLC[fermiLC['ts'] < tsUpperLimit]
        fermiLC = fermiLC[fermiLC['ts'] >= tsUpperLimit]

        priorFile = '{}_{}_{}_result.txt'.format('fermi', fermiBinning, nSims)
        gamma = lcUtil.readGammaFromFile(self.logStdout, 'mwlLightCurves', priorFile,
                                         sourceNowShort, configAnal['bayesChange'])
        self.logStdout.info([['p', sourceNow],
                            ['wb', ' {} {} gamma -'.format('fermi', fermiBinning)],
                            ['b', gamma]])

        fermiBlocks = bayesian_blocks(fermiLC['tmax_mjd'], fermiLC['flux'],
                                      fermiLC['flux_err'], fitness='measures', gamma=gamma)
        fermiFluxMean, fermiFluxErrMean = lcUtil.calcBlockFlux(fermiLC['tmax_mjd'],
                                                               fermiLC['flux'],
                                                               fermiLC['flux_err'],
                                                               fermiBlocks)

        rebinSwift = veritasBinning

        if rebinSwift == 'No':
            swiftDataWindow = swiftData[swiftData['mode'] == 'WT']
            swiftDataPhoton = swiftData[swiftData['mode'] == 'PC']

        priorFile = '{}_{}_{}_result.txt'.format('swift', veritasBinning, nSims)
        gamma = lcUtil.readGammaFromFile(self.logStdout, 'mwlLightCurves', priorFile,
                                         sourceNowShort, configAnal['bayesChange'])
        self.logStdout.info([['p', sourceNow],
                            ['wb', ' {} {} gamma -'.format('swift', veritasBinning)],
                            ['b', gamma]])

        swiftRateErrorAverage = (swiftData['Rate error down'] + swiftData['Rate error up'])/2.
        swiftBlocks = bayesian_blocks(swiftData['Date'], swiftData['Rate'],
                                      swiftRateErrorAverage, fitness='measures', gamma=gamma)

        swiftFluxMean, swiftFluxErrMean = lcUtil.calcBlockFlux(swiftData['Date'],
                                                               swiftData['Rate'],
                                                               swiftRateErrorAverage,
                                                               swiftBlocks)

        # Correct for wrong start possitons of edges.
        # Implemented output of bayesian blocks is in the middle between two observations.
        # However no statement can be made for a change between two measurments
        # so that a new block should only begin at the start of a observation.
        correctedEdges = np.array([])
        for edge in swiftBlocks:
            aboveEdge = swiftData['Date'] >= edge
            correctedEdges = np.append(correctedEdges,
                                       min(swiftData['Date'][aboveEdge] -
                                           (swiftData['Date error down'][aboveEdge] +
                                            swiftData['Date error up'][aboveEdge])/2))

        swiftBlocks = correctedEdges

        fig = plt.figure(figsize=(25, 15))
        ax1 = fig.add_subplot(311)

        xText = configAnal['plot']['mwlLC']['xText']
        yText = configAnal['plot']['mwlLC']['yText']
        xICRC = configAnal['plot']['mwlLC']['xICRC']
        yICRC = configAnal['plot']['mwlLC']['yICRC']

        # FIXME Why in ED the unit is /cm^2/s and here it has ergs as well?
        xTitle = 'MJD'
        fVeritas = configAnal['fluxScale']['veritas']
        yUnit = r'[cm$^{-2} s^{-1}$]'
        yFactor = ''
        if fVeritas > 0:
            yScaling = r'($\times 10^{{{0:d}}}$)'.format(int(np.log10(fVeritas)))
        yTitle = r'$\mathcal{{F}}$ {} {}'.format(yScaling, yUnit)
        plt.errorbar(veritasFluxes['DateMJD'],
                     fVeritas*veritasFluxes['Flux'],
                     xerr=veritasFluxes['Date Error'],
                     yerr=fVeritas*veritasFluxes['Flux Error'],
                     fmt=self.markers[0], color=self.colors[0], mfc='none',
                     markersize=self.markersize, elinewidth=self.elinewidth,
                     zorder=4, label=r'VERITAS')
        plt.xlabel(xTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'], labelpad=0)
        plt.ylabel(yTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'])
        plt.gca().get_yaxis().get_offset_text().set_fontsize(configAnal['plot']
                                                                       ['mwlLC']
                                                                       ['axesLabelFont'])
        plt.text(xText, yText, plotInfo,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=configAnal['plot']['mwlLC']['textFont'],
                 transform=plt.gca().transAxes)

        plt.text(xICRC, yICRC, 'ICRC 2019',
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=configAnal['plot']['mwlLC']['textFont'],
                 transform=plt.gca().transAxes)

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

        if configAnal['fVar']:
            veritasFVar = fVar(veritasFluxes['Flux'], veritasFluxes['Flux Error'])
            veritasDeltaFVar = deltaFVar(veritasFluxes['Flux'], veritasFluxes['Flux Error'])
            plt.text(0.01, 0.03,
                     r'$F_{\mathrm{var}} = %1.2f \pm %1.2f$' % (veritasFVar, veritasDeltaFVar),
                     horizontalalignment='left',
                     verticalalignment='bottom',
                     fontsize=configAnal['plot']['mwlLC']['textFont'],
                     transform=plt.gca().transAxes)

        if configAnal['constFit']:
            conF, conFcov = curve_fit(lcUtil.conFunc, veritasFluxes['DateMJD'],
                                      fVeritas*veritasFluxes['Flux'],
                                      [fVeritas*veritasFluxes['Flux'][0]],
                                      fVeritas*veritasFluxes['Flux Error'])
            conChi2 = np.sum(((fVeritas*veritasFluxes['Flux'] -
                               lcUtil.conFunc(veritasFluxes['DateMJD'], *conF))**2) /
                             ((fVeritas*veritasFluxes['Flux Error'])**2))
            plt.plot(veritasFluxes['DateMJD'], len(veritasFluxes['DateMJD'])*list(conF),
                     '--', color='black', linewidth=3,
                     label=r'Constant Fit ($p$-value = %1.2e)' %
                     (1 - stats.chi2.cdf(conChi2, (len(veritasFluxes['DateMJD']) - 1))))

        plt.legend(loc='upper left', fontsize=configAnal['plot']['mwlLC']['legendFont'])
        plt.tick_params(axis='both', which='major',
                        labelsize=configAnal['plot']['mwlLC']['axesLabelFont'])

        plt.grid(True)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()

        fig.add_subplot(312, sharex=ax1)

        plotInfo = ('{}, Fermi-LAT {}\n'
                    '{} GeV < E < {} GeV').format(sourceNow,
                                                  fermiBinningForPlotInfo + ' bins',
                                                  0.2,
                                                  300)

        xTitle = 'MJD'
        fFermi = configAnal['fluxScale']['fermi']
        yUnit = r'[cm$^{-2} s^{-1}$]'
        yFactor = ''
        if fFermi > 0:
            yScaling = r'($\times 10^{{{0:d}}}$)'.format(int(np.log10(fFermi)))
        yTitle = r'$\mathcal{{F}}$ {} {}'.format(yScaling, yUnit)
        plt.errorbar(fermiLC['tmax_mjd'],
                     fFermi*fermiLC['flux'],
                     yerr=fFermi*fermiLC['flux_err'],
                     fmt=self.markers[1], color=self.colors[1], mfc='none',
                     markersize=self.markersize, elinewidth=self.elinewidth,
                     label='Fermi-LAT', zorder=4)
        plt.errorbar(fermiUpperLimits['tmax_mjd'],
                     fFermi*fermiUpperLimits['flux_ul95'],
                     yerr=0.6*fFermi*np.mean(fermiUpperLimits['flux_ul95']),
                     fmt='_', color=self.colors[1],
                     markersize=self.markersize, elinewidth=self.elinewidth,
                     capsize=self.capsize, uplims=True, zorder=4)
        plt.text(xText, yText, plotInfo,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=configAnal['plot']['mwlLC']['textFont'],
                 transform=plt.gca().transAxes)

        # Plot Bayesian blocks
        plt.gca().step(x=fermiBlocks, y=fFermi*fermiFluxMean, linewidth=3,
                       where='post', color='darkgrey', alpha=0.8, zorder=3)
        for i in range(len(fermiBlocks)-1):
            x1 = [fermiBlocks[i], fermiBlocks[i+1], fermiBlocks[i+1], fermiBlocks[i]]
            y1 = [fFermi*(fermiFluxMean[i] + fermiFluxErrMean[i]),
                  fFermi*(fermiFluxMean[i] + fermiFluxErrMean[i]),
                  fFermi*(fermiFluxMean[i] - fermiFluxErrMean[i]),
                  fFermi*(fermiFluxMean[i] - fermiFluxErrMean[i])]
            label = None
            if i == 0:
                label = r'Bayesian blocks'
            plt.gca().fill(x1, y1, color='darkgrey', label=label, alpha=0.4, zorder=2)

        if configAnal['fVar']:
            fermiFVar = lcUtil.fVar(fermiLC['flux'], fermiLC['flux_err'])
            fermiDeltaFVar = lcUtil.deltaFVar(fermiLC['flux'], fermiLC['flux_err'])
            plt.text(0.01, 0.03,
                     r'$F_{\mathrm{var}} = %1.2f \pm %1.2f$' % (fermiFVar, fermiDeltaFVar),
                     horizontalalignment='left',
                     verticalalignment='bottom',
                     fontsize=configAnal['plot']['mwlLC']['textFont'],
                     transform=plt.gca().transAxes)

        if configAnal['constFit']:
            conF, conFcov = curve_fit(lcUtil.conFunc, fermiLC['tmax_mjd'], fFermi*fermiLC['flux'],
                                      [fFermi*fermiLC['flux'][0]], fFermi*fermiLC['flux_err'])
            conChi2 = np.sum(((fFermi*fermiLC['flux'] -
                               lcUtil.conFunc(fermiLC['tmax_mjd'], *conF))**2) /
                             ((fFermi*fermiLC['flux_err'])**2))
            plt.plot(fermiLC['tmax_mjd'], len(fermiLC['tmax_mjd'])*list(conF),
                     '--', color='forestgreen', linewidth=3,
                     label=r'Constant Fit ($p$-value = %1.2e)' %
                     (1 - stats.chi2.cdf(conChi2, (len(fermiLC['tmax_mjd']) - 1))))

        plt.xlabel(xTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'], labelpad=0)
        plt.ylabel(yTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'])
        plt.gca().get_yaxis().get_offset_text().set_fontsize(configAnal['plot']
                                                                       ['mwlLC']
                                                                       ['axesLabelFont'])
        plt.legend(loc='upper left', fontsize=configAnal['plot']['mwlLC']['legendFont'])
        plt.tick_params(axis='both', which='major',
                        labelsize=configAnal['plot']['mwlLC']['axesLabelFont'])
        plt.grid(True)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()

        fig.add_subplot(313, sharex=ax1)

        xTitle = 'MJD'
        yTitle = r'Rate [counts $s^{-1}$]'
        plotInfo = ('{}, Swift-XRT {}\n'
                    '{} keV < E < {} keV').format(sourceNow,
                                                  '{} bins'.format(rebinSwift),
                                                  0.3,
                                                  10)

        if rebinSwift == 'No':
            plt.errorbar(swiftDataWindow['Date'],
                         swiftDataWindow['Rate'],
                         xerr=[swiftDataWindow['Date error down'],
                               swiftDataWindow['Date error up']],
                         yerr=[swiftDataWindow['Rate error down'],
                               swiftDataWindow['Rate error up']],
                         fmt=self.markers[2], color=self.colors[2], mfc='none',
                         markersize=self.markersize, elinewidth=self.elinewidth,
                         label='Windowed Timing', zorder=5)
            plt.errorbar(swiftDataPhoton['Date'],
                         swiftDataPhoton['Rate'],
                         xerr=[swiftDataPhoton['Date error down'],
                               swiftDataPhoton['Date error up']],
                         yerr=[swiftDataPhoton['Rate error down'],
                               swiftDataPhoton['Rate error up']],
                         fmt=self.markers[3], color=self.colors[3], mfc='none',
                         markersize=self.markersize, elinewidth=self.elinewidth,
                         label='Photon Counting', zorder=4)  # tomato color is also OK
        else:
            plt.errorbar(swiftData['Date'],
                         swiftData['Rate'],
                         xerr=[swiftData['Date error down'], swiftData['Date error up']],
                         yerr=[swiftData['Rate error down'], swiftData['Rate error up']],
                         fmt=self.markers[2], color=self.colors[2], mfc='none',
                         markersize=self.markersize, elinewidth=self.elinewidth,
                         label='Swift-XRT', zorder=5)

        # Plot Bayesian blocks
        plt.gca().step(x=swiftBlocks, y=swiftFluxMean, zorder=3, linewidth=3,
                       where='post', color='darkgrey', alpha=0.8)

        for i in range(len(swiftBlocks)-1):
            x1 = [swiftBlocks[i], swiftBlocks[i+1], swiftBlocks[i+1], swiftBlocks[i]]
            y1 = [(swiftFluxMean[i] + swiftFluxErrMean[i]),
                  (swiftFluxMean[i] + swiftFluxErrMean[i]),
                  (swiftFluxMean[i] - swiftFluxErrMean[i]),
                  (swiftFluxMean[i] - swiftFluxErrMean[i])]
            label = None
            if i == 0:
                label = r'Bayesian blocks'
            plt.gca().fill(x1, y1, color='darkgrey', label=label, alpha=0.4, zorder=2)

        plt.xlabel(xTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'], labelpad=0)
        plt.ylabel(yTitle, fontsize=configAnal['plot']['mwlLC']['axesTitleFont'])
        plt.gca().get_yaxis().get_offset_text().set_fontsize(configAnal['plot']
                                                                       ['mwlLC']
                                                                       ['axesLabelFont'])
        plt.text(xText, yText, plotInfo,
                 horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=configAnal['plot']['mwlLC']['textFont'],
                 transform=plt.gca().transAxes)

        if configAnal['fVar']:
            swiftFVar = lcUtil.fVar(swiftData['Rate'], swiftData['Rate error down'])
            swiftDeltaFVar = lcUtil.deltaFVar(swiftData['Rate'], swiftData['Rate error down'])
            plt.text(0.01, 0.03,
                     r'$F_{\mathrm{var}} =%1.2f \pm %1.2f$' % (swiftFVar, swiftDeltaFVar),
                     horizontalalignment='left',
                     verticalalignment='bottom',
                     fontsize=configAnal['plot']['mwlLC']['textFont'],
                     transform=plt.gca().transAxes)

        if configAnal['constFit']:
            conF, conFcov = curve_fit(lcUtil.conFunc, swiftData['Date'], swiftData['Rate'],
                                      [swiftData['Rate'][0]], swiftData['Rate error down'])
            conChi2 = np.sum(((swiftData['Rate'] - lcUtil.conFunc(swiftData['Date'], *conF))**2) /
                             ((swiftData['Rate error down'])**2))
            plt.plot(swiftData['Date'], len(swiftData['Date'])*list(conF),
                     '--', color='dodgerblue', linewidth=3,
                     label=r'Constant Fit ($p$-value = %1.2e)' %
                     (1 - stats.chi2.cdf(conChi2, (len(swiftData['Date']) - 1))))

        plt.legend(loc='upper left', fontsize=configAnal['plot']['mwlLC']['legendFont'])
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
        fig.subplots_adjust(hspace=0)

        outputFile = 'mwlLightCurves/{}_{}_lightcurve.pdf'.format(veritasBinning, sourceNowShort)
        plt.savefig(outputFile, bbox_inches='tight')

        plt.close('all')

        return


if __name__ == '__main__':

    # FIXME implement also the determine prior code as part of the full package.

    mwl = mwl()
    logStdout = mwl.getLogger()

    parser = argparse.ArgumentParser(description=('Compare the old and new spectra.'))
    parser.add_argument('-c', '--compareCorrectedLC', action='store_true', default=False,
                        help='Compare corrected lightcurves')

    args = parser.parse_args()

    configAnalFile = 'configAnalysis.yaml'
    configAnal = gUtil.readYamlFile(logStdout, configAnalFile)

    # FIXME - make this local and move to configAnalysis.yaml
    corrTableFile = ('/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/'
                     'crabStability/plotLC/correctionFactors.txt')
    corrTable = lcUtil.readCorrTable(corrTableFile)

    for sourceNowShort, sourceNow in configAnal['sources'].items():

        for veritasBinning, fermiBinning in configAnal['binning'].items():

            logStdout.info([['wb', 'Running'],
                            ['p', ' {}'.format(sourceNow)],
                            ['wb', ' {}'.format(veritasBinning)]])

            veritasDirectory = os.path.join('./makeLC/', sourceNowShort)
            veritasFile = glob.glob(os.path.join(veritasDirectory,
                                    '{}*fullEnergyRange*.txt'.format(veritasBinning)))[0]

            veritasLC = lcUtil.readVeritasLC(logStdout, veritasFile)

            if args.compareCorrectedLC:
                mwl.compareCorrectedLC(veritasLC, corrTable,
                                       sourceNow, veritasBinning, configAnal)

            veritasLC = lcUtil.correctFluxesFromCrabLC(veritasLC, corrTable)

            # FIXME
            # For now put the maximum dates hardcoded here
            if (configAnal['dates']['veritas'][0] != 54375
                    or configAnal['dates']['veritas'][1] != 58465):
                veritasLC = lcUtil.lcBetweenDates(logStdout, veritasLC,
                                                  configAnal['dates']['veritas'], 'veritas')
            veritasObsFile = os.path.join(os.path.join('./makeSpectrum/', sourceNowShort),
                                          'fluxPerRun.txt')
            veritasObs = lcUtil.readVeritasObservations(logStdout, veritasObsFile)

            fermiBinningNow = fermiBinning
            if fermiBinning == 'yearly':
                fermiBinningNow = 'monthly'
            fermiDirNow = os.path.join(configAnal['fermi']['baseDir'],
                                       '{}LightCurves'.format(fermiBinningNow),
                                       sourceNowShort, sourceNowShort)
            fermiFileName = gUtil.getSourceNameFermi(logStdout, sourceNow) + '_lightcurve.npy'
            fermiFile = os.path.join(fermiDirNow, fermiFileName)
            fermiLC = lcUtil.readFermiLC(logStdout, fermiFile)
            if fermiBinning == 'yearly':
                fermiLC = lcUtil.rebinFermi(fermiLC, veritasObs['date'],
                                            configAnal['dates']['fermi'])
            # FIXME
            # For now put the maximum dates hardcoded here
            if (configAnal['dates']['fermi'][0] != 53423
                    or configAnal['dates']['fermi'][1] != 58465):
                fermiLC = lcUtil.lcBetweenDates(logStdout, fermiLC,
                                                configAnal['dates']['fermi'], 'fermi')

            if sourceNowShort in ['PG1553', '1ES1011']:
                # Remove crazy points in finely binned LCs
                # (affects one point in PG1553 and one in 1ES1011)
                fermiLC = fermiLC[np.isfinite(fermiLC['flux'])]
                fermiLC = fermiLC[np.where(fermiLC['flux_err'] != np.max(fermiLC['flux_err']))]
                fermiLC = fermiLC[[fermiLC['flux_err']/fermiLC['flux'] < 1]]
                tsUpperLimit = 0.001
                fermiLC = fermiLC[fermiLC['ts'] >= tsUpperLimit]

            swiftFile = os.path.join(configAnal['swift']['baseDir'], sourceNowShort,
                                     'dailyBins', '{}_lightcurve.qdp'.format(sourceNowShort))

            swiftData = lcUtil.readSwiftLC(logStdout, swiftFile)
            swiftData = lcUtil.rebinSwiftLC(logStdout, swiftData, veritasBinning,
                                            veritasObs['date'], configAnal['dates']['swift'])

            # FIXME
            # For now put the maximum dates hardcoded here
            if (configAnal['dates']['swift'][0] != 53423
                    or configAnal['dates']['swift'][1] != 58465):
                swiftData = lcUtil.lcBetweenDates(logStdout, swiftData,
                                                  configAnal['dates']['swift'], 'swift')

            if veritasBinning in configAnal['dcf'][sourceNowShort]['binning']:
                dcfFilePrefix = '{}_{}_{}'.format(sourceNowShort, veritasBinning, 'veritas')
                dcfUtil.saveDataForDCF('mwlLightCurves', dcfFilePrefix,
                                       veritasLC['DateMJD'], veritasLC['Flux'],
                                       veritasLC['Flux Error'])
                dcfFilePrefix = '{}_{}_{}'.format(sourceNowShort, veritasBinning, 'fermi')
                tsUpperLimit = 0.001
                fermiLC_forDCF = fermiLC[fermiLC['ts'] >= tsUpperLimit]
                dcfUtil.saveDataForDCF('mwlLightCurves', dcfFilePrefix,
                                       fermiLC_forDCF['tmax_mjd'], fermiLC_forDCF['flux'],
                                       fermiLC_forDCF['flux_err'])
                dcfFilePrefix = '{}_{}_{}'.format(sourceNowShort, veritasBinning, 'swift')
                dcfUtil.saveDataForDCF('mwlLightCurves', dcfFilePrefix,
                                       swiftData['Date'], swiftData['Rate'],
                                       swiftData['Rate error average'])

            if not args.compareCorrectedLC:
                mwl.plotLC(veritasLC, fermiLC, swiftData,
                           sourceNow, veritasBinning, fermiBinning,
                           configAnal['nBayesSims'][sourceNowShort], configAnal)

        for binningNow in configAnal['dcf'][sourceNowShort]['binning']:
            dcfNow = mwl.dcfForSource(sourceNowShort, binningNow,
                                      configAnal['dcf'][sourceNowShort]['minPeriod'],
                                      configAnal['dcf'][sourceNowShort]['maxPeriod'],
                                      configAnal['dcf'][sourceNowShort]['binWidth'])
            mwl.plotDCF(dcfNow, sourceNowShort, binningNow, configAnal)

    if args.compareCorrectedLC:
        # Merge lightcurve plots into one PDF
        fileToMerge = list()
        prefix = 'mwlLightCurves/compareCorrectedLC'
        for binning in configAnal['binning'].keys():
            for sourceNowShort in configAnal['sources'].keys():
                fileToMerge.append('{}/{}_{}_compareCorrectedLC.pdf'.format(prefix,
                                                                            binning,
                                                                            sourceNowShort))
            gUtil.mergePDFs(logStdout, fileToMerge, '{}/{}CompareLCs.pdf'.format(prefix, binning))
            fileToMerge.clear()
    else:
        # Merge lightcurve plots into one PDF
        fileToMerge = list()
        prefix = 'mwlLightCurves'
        for binning in configAnal['binning'].keys():
            for sourceNowShort in configAnal['sources'].keys():
                fileToMerge.append('{}/{}_{}_lightcurve.pdf'.format(prefix,
                                                                    binning,
                                                                    sourceNowShort))
            gUtil.mergePDFs(logStdout, fileToMerge, '{}/{}LCs.pdf'.format(prefix, binning))
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
                            '{}/mwl_{}_dcf.pdf'.format(prefix, binningNow))
            fileToMerge[binningNow].clear()
