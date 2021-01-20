#!/usr/bin/python

import os.path
import numpy as np
import ROOT
import argparse
import sys
from collections import defaultdict
sys.path.append("..")
from lib import colourLogger
from lib import generalUtil as gUtil


def calcSpectrum(source, configAnal):

    anasumFile = configAnal['anasumFiles'][source]
    dirNow = os.path.join(os.getcwd(), source)
    sourceTitle = configAnal['sources'][source]

    parameters = dict()
    parameters['minEnergy'] = configAnal['spectralParameters'][source][0]
    parameters['normalizationEnergy'] = configAnal['spectralParameters'][source][1]
    parameters['spectralIndex'] = configAnal['spectralParameters'][source][2]
    parameters['maxEnergy'] = configAnal['spectralParameters'][source][3]

    vpah = ROOT.VPlotAnasumHistograms(anasumFile)
    ca1 = vpah.plot_radec(0)  # plot type = 0 for significance map

    # We can also add more things to the plot, like the locations of stars
    # the 'Hipparcos...' is a catalogue kept in $VERITAS_EVNDISP_AUX_DIR/AstroData/Catalogues/
    # You can pick other catalogues from there, but the Hipparchos MAG9 is the preferred
    vpah.plot_catalogue(ca1, 'Hipparcos_MAG9_1997.dat')
    vpah.plot_catalogue(ca1, 'tevcat.dat')

    # When the 'background' or OFF events are calculated, certain regions are excluded
    # this will plot those excluded regions
    vpah.plot_excludedRegions(ca1)

    w = 500
    h = 500
    # This will set the size of the plot, in pixels
    # when printed to a file, this is how big the image file will be
    ca1.SetCanvasSize(w, h)
    # This will set the size of the plot's window, for when it displays
    ca1.SetWindowSize(w + (w - ca1.GetWw()), h + (h - ca1.GetWh()))

    titleText = ROOT.TLatex(0.3, 0.92, 'Significance map - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    ca1.Modified()
    ca1.Update()

    ca1.Print(os.path.join(dirNow, 'significanceMap.pdf'))

    ca2 = vpah.plot_radec(1)  # plot type = 1 for excess map

    # When the 'background' or OFF events are calculated, certain regions are excluded
    # this will plot those excluded regions
    vpah.plot_excludedRegions(ca2)

    ca2.SetCanvasSize(w, h)
    ca2.SetRightMargin(0.15)
    # pl2.GetZaxis()->SetTitleOffset(1.3)
    ca2.SetWindowSize(w + (w - ca2.GetWw()), h + (h - ca2.GetWh()))

    titleText = ROOT.TLatex(0.3, 0.92, 'Excess map - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    ca2.Modified()
    ca2.Update()
    ca2.Print(os.path.join(dirNow, 'excessMap.pdf'))

    sigCanvas = vpah.plot_significanceDistributions()
    sigCanvas.SetCanvasSize(w, h)
    sigCanvas.SetWindowSize(w + (w - sigCanvas.GetWw()), h + (h - sigCanvas.GetWh()))
    sigCanvas.cd()

    titleText = ROOT.TLatex(0.05, 0.92, 'Significance distribution - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    sigCanvas.Modified()
    sigCanvas.Update()
    sigCanvas.Print(os.path.join(dirNow, 'sigDist.pdf'))

    myRunSummary = ROOT.VPlotRunSummary(anasumFile)
    cumulativeCanvas = myRunSummary.plot_cumSignificance()

    cumulativeCanvas.SetCanvasSize(w, h)
    cumulativeCanvas.SetWindowSize(w + (w - cumulativeCanvas.GetWw()),
                                   h + (h - cumulativeCanvas.GetWh()))

    titleText = ROOT.TLatex(0.3, 0.92, 'Cumalative significance - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    cumulativeCanvas.Modified()
    cumulativeCanvas.Update()
    cumulativeCanvas.Print(os.path.join(dirNow, 'cumulativeSig.pdf'))

    ################################################################################################
    # Plot energy spectrum
    ################################################################################################

    logStdout.info([['p', source],
                    ['wb', 'spectrum,'],
                    ['g', '{} < E < {} TeV'.format(parameters['minEnergy'],
                                                   parameters['maxEnergy'])]])

    vEnergy = ROOT.VEnergySpectrum(anasumFile)
    # vEnergy.setDebug(3)
    vEnergy.setEnergyRangeLinear(parameters['minEnergy'], 30.)
    vEnergy.setPlottingEnergyRangeLinear(0.10, 4.)
    vEnergy.setPlottingYaxis(1.e-15, 1.e-9)
    vEnergy.setEnergyBinning(0.2)
    vEnergy.setSpectralFitFunction(0)

    # Setting desired confidence level (needs to be set before applying fit)
    vEnergy.setConfidenceLevel(0.68)

    vEnergy.setEnergyThresholdDefinition(2, -1., 0.1)

    vEnergy.setSpectralFitFluxNormalisationEnergy(parameters['normalizationEnergy'])
    vEnergy.setSpectralFitRangeLin(parameters['minEnergy'], 2.5)
    # vEnergy.setSignificanceParameters(1.999, 4.999, 0.99, 17, 0)
    vEnergy.setSignificanceParameters(-9, -9, 0.99, 17, 0)
    vEnergy.setPlottingMultiplierIndex(0.0)

    c_Spectrum = vEnergy.plot()

    fEnergy = vEnergy.fitEnergySpectrum()
    # Getting TGraphError Confidence Interval (needs to be called after fit)
    gConInt = vEnergy.getEnergySpectrumConfidenceInterval()

    fEnergy.SetLineColor(ROOT.kRed)
    # Making it semi Transparent
    gConInt.SetFillColorAlpha(ROOT.kRed, 0.5)
    gConInt.SetFillStyle(1001)
    # Fills area between error bars
    gConInt.Draw('3')
    vEnergy.plotFitValues()

    titleText = ROOT.TLatex(0.3, 0.92, 'Energy spectrum - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    vEnergy.printDifferentialFluxes()

    c_Spectrum.Print(os.path.join(dirNow, 'energySpectrum.pdf'))

    # --------------------------------------------------------------------------------------
    # Plot also with E^2
    # --------------------------------------------------------------------------------------

    vEnergy.setPlottingMultiplierIndex(2.0)

    c_Spectrum = vEnergy.plot()

    fEnergy = vEnergy.fitEnergySpectrum()
    # Getting TGraphError Confidence Interval (needs to be called after fit)
    gConInt = vEnergy.getEnergySpectrumConfidenceInterval()

    fEnergy.SetLineColor(ROOT.kRed)
    # Making it semi Transparent
    gConInt.SetFillColorAlpha(ROOT.kRed, 0.5)
    gConInt.SetFillStyle(1001)
    # Fills area between error bars
    gConInt.Draw('3')
    vEnergy.plotFitValues()

    titleText = ROOT.TLatex(0.3, 0.92, 'Energy spectrum - {}'.format(sourceTitle))
    titleText.SetNDC()
    titleText.SetTextSize(0.031)
    titleText.SetTextColor(ROOT.kBlack)
    titleText.Draw()

    c_Spectrum.Print(os.path.join(dirNow, 'energySpectrumSquared.pdf'))

    ##############################################################################################

    # Write out the info of each run to the log file (weird way of doing it...)
    vEnergy = ROOT.VEnergySpectrum(anasumFile)
    vEnergy.setDebug(1)
    vEnergy.setEnergyRangeLinear(parameters['minEnergy'], 30.)
    vEnergy.setEnergyThresholdDefinition(2, -1., 0.1)
    vEnergy.combineRuns()

    ##############################################################################################

    fluxCalc = ROOT.VFluxCalculation(anasumFile)
    # set spectral parameters (start of flux calculation [TeV],
    # normalization at energy E0 [TeV], spectral index, end of flux calculation [TeV])
    fluxCalc.setSpectralParameters(parameters['minEnergy'], parameters['normalizationEnergy'],
                                   parameters['spectralIndex'], parameters['maxEnergy'])
    fluxCalc.setSignificanceParameters(2)
    fluxCalc.setTimeBinnedAnalysis(True)
    fluxCalc.calculateIntegralFlux(parameters['minEnergy'])
    # fluxCalc.plotFluxesVSMJDDaily();
    fluxCalc.writeResults(os.path.join(dirNow, 'fluxPerRun.root'))

    fluxPerRunFile = ROOT.TFile.Open(os.path.join(dirNow, 'fluxPerRun.root'), 'r')

    observations = defaultdict(list)
    for event in fluxPerRunFile.fluxes:
        if event.Run > 0:
            observations['run'].append(event.Run)
            observations['date'].append(event.MJD)
            observations['flux'].append(event.Flux)
            observations['fluxError'].append(event.FluxE)
            observations['significance'].append(event.Signi)
            observations['ze'].append(event.Ze)

    obsData = np.c_[observations['run'],
                    observations['date'],
                    observations['flux'],
                    observations['fluxError'],
                    observations['significance'],
                    observations['ze']]

    headersType = {'names': ('run', 'date', 'flux', 'fluxError',
                             'significance', 'ze'),
                   'formats': ('f8', 'f8', 'f8', 'f8',
                               'f8', 'f8')}

    obsData = np.core.records.fromarrays(obsData.transpose(), dtype=headersType)

    np.savetxt(os.path.join(dirNow, 'fluxPerRun.txt'), obsData, delimiter='    ',
               fmt=['%d', '%.1f', '%.4e', '%.4e', '%.3f', '%.1f'],
               header='     '.join(obsData.dtype.names))

    del vpah

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Produce skymaps and spectra.')
    parser.add_argument('source', action='store',
                        help='Source to produce the spectrum for')
    parser.add_argument('configAnalFile', action='store',
                        help='YAML file with the analysis configuration')

    args = parser.parse_args()

    logStdout = colourLogger.initStdOutLogger()

    configAnal = gUtil.readYamlFile(logStdout, args.configAnalFile)

    ROOT.gROOT.SetBatch(True)
    ROOT.TFile.Open._creates = True

    # load shared library
    ROOT.gSystem.Load("$EVNDISPSYS/lib/libVAnaSum.so")

    calcSpectrum(args.source, configAnal)
