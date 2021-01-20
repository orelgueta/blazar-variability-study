#!/usr/bin/python

import os.path
import os
import sys
import glob
import subprocess
import numpy as np
import matplotlib as mlp
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.legend_handler import HandlerTuple
from matplotlib.legend_handler import HandlerLine2D
import scipy.constants as con
import argparse
import filecmp
from scipy import interpolate
from astropy.io import ascii
from astropy import units as u
from astropy.table import QTable
from lib import colourLogger
from lib import generalUtil as gUtil
from lib import plotUtil as pUtil
from lib import sedUtil as sUtil


class sed(object):
    """
    A class to study and plot MWL data
    """

    def __init__(self, *args, **kwargs):

        self.logStdout = colourLogger.initStdOutLogger()
        pu = pUtil.plotUtil()
        self.colors = pu.getColors()
        self.markers = pu.getMarkers()
        self.lines = pu.getLines()
        self.fontsize = pu.getFontsize()
        self.markersize = pu.getMarkersize()
        self.elinewidth = pu.getElinewidth()
        self.capsize = pu.getCapsize()

        self.deabsorbFactorVeritas = list()
        self.deabsorbFactorVeritasDn = list()
        self.deabsorbFactorVeritasUp = list()
        self.deabsorbFactorFermi = list()
        self.deabsorbFactorFermiDn = list()
        self.deabsorbFactorFermiUp = list()

    def getLogger(self):
        return self.logStdout

    def saveSEDs(self, vFluxes, sedFermi, sourceNowShort):

        veritasSed = np.c_[vFluxes['E'],
                           vFluxes['E_min'],
                           vFluxes['E_max'],
                           vFluxes['dE'],
                           vFluxes['dN/dE'],
                           vFluxes['dN/dE error down'],
                           vFluxes['dN/dE error up']]

        headersType = {'names': ('E', 'E_min', 'E_max', 'dE',
                                 'dN/dE', 'dN/dE error down', 'dN/dE error up'),
                       'formats': ('f8', 'f8', 'f8', 'f8',
                                   'f8', 'f8', 'f8')}

        veritasSed = np.core.records.fromarrays(veritasSed.transpose(), dtype=headersType)

        np.savetxt('sed/veritasSed/{}_veritas_sed.dat'.format(sourceNowShort),
                   veritasSed, fmt=['%e']*len(veritasSed.dtype.names),
                   header='     '.join(veritasSed.dtype.names))

        fermiSed = np.c_[sedFermi['e_ctr'],
                         sedFermi['dnde'],
                         sedFermi['dnde_err_lo'],
                         sedFermi['dnde_err_hi']]

        headersType = {'names': ('E',
                                 'dN/dE', 'dN/dE error down', 'dN/dE error up'),
                       'formats': ('f8', 'f8', 'f8', 'f8')}

        fermiSed = np.core.records.fromarrays(fermiSed.transpose(), dtype=headersType)

        np.savetxt('sed/fermiSed/{}_fermi_sed.dat'.format(sourceNowShort),
                   fermiSed, fmt=['%e']*len(fermiSed.dtype.names),
                   header='     '.join(fermiSed.dtype.names))

        return

    def writeEcsvTable(self, ecsvFile, columns):
        """
        This function recieves a dict of columns with the keys
        as the titles of the columns. For each key (column), a list of values
        is given, where it is assumed that the unit is already included in it.
        """

        tableToWrite = QTable()
        for title, values in columns.items():
            tableToWrite[title] = values

        ascii.write(tableToWrite, ecsvFile, format='ecsv', overwrite=True)
        del tableToWrite

        return

    def saveEcsvSED(self, sourceNowShort, vheData, heData, xrayData, optData=list()):

        ############################################################################################
        # VERITAS
        ############################################################################################

        sed = dict()
        sed['energy'] = (vheData['E']*u.TeV).to('eV')
        sed['energy_error'] = (vheData['dE']*u.TeV).to('eV')
        sed['flux'] = (vheData['dN/dE']*(1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_lo'] = (vheData['dN/dE error down'] *
                                (1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_hi'] = (vheData['dN/dE error up'] *
                                (1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['ul'] = 0*vheData['E'].astype(np.int)

        self.writeEcsvTable('sed/forNaima/{0}/veritas_{0}.ecsv'.format(sourceNowShort), sed)

        ############################################################################################
        # Deabsorbed VERITAS
        ############################################################################################

        # FIXME notice that the uncertainty on the deabsorption is not taken into account here.
        sed = dict()
        sed['energy'] = (vheData['E']*u.TeV).to('eV')
        sed['energy_error'] = (vheData['dE']*u.TeV).to('eV')
        sed['flux'] = (self.deabsorbFactorVeritas *
                       vheData['dN/dE']*(1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_lo'] = (self.deabsorbFactorVeritas *
                                vheData['dN/dE error down'] *
                                (1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_hi'] = (self.deabsorbFactorVeritas *
                                vheData['dN/dE error up'] *
                                (1./(u.TeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['ul'] = 0*vheData['E'].astype(np.int)

        self.writeEcsvTable('sed/forNaima/{0}/veritas_deabsorbed_{0}.ecsv'.format(sourceNowShort),
                            sed)

        ############################################################################################
        # Fermi
        ############################################################################################

        sed = dict()
        sed['energy'] = (heData['e_ctr']*u.MeV).to('eV')
        sed['energy_error'] = (0.*heData['e_ctr']*u.MeV).to('eV')
        sed['flux'] = (heData['dnde']*(1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_lo'] = (0.5*heData['dnde_err'] *
                                (1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_hi'] = (0.5*heData['dnde_err'] *
                                (1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['ul'] = 0*heData['e_ctr'].astype(np.int)

        self.writeEcsvTable('sed/forNaima/{0}/fermi_{0}.ecsv'.format(sourceNowShort), sed)

        ############################################################################################
        # Fermi deabsorbed
        ############################################################################################

        sed = dict()
        sed['energy'] = (heData['e_ctr']*u.MeV).to('eV')
        sed['energy_error'] = (0.*heData['e_ctr']*u.MeV).to('eV')
        sed['flux'] = (self.deabsorbFactorFermi * heData['dnde'] *
                       (1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_lo'] = (0.5 * self.deabsorbFactorFermi * heData['dnde_err'] *
                                (1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['flux_error_hi'] = (0.5 * self.deabsorbFactorFermi * heData['dnde_err'] *
                                (1./(u.MeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
        sed['ul'] = 0*heData['e_ctr'].astype(np.int)

        self.writeEcsvTable('sed/forNaima/{0}/fermi_deabsorbed_{0}.ecsv'.format(sourceNowShort),
                            sed)

        ############################################################################################
        # Swift-XRT
        ############################################################################################

        if len(xrayData) > 0:
            sed = dict()
            sed['energy'] = (xrayData['E']*u.GeV).to('eV')
            sed['energy_error'] = (((xrayData['dE down'] + xrayData['dE up'])/2.) *
                                   u.GeV).to('eV')
            sed['flux'] = ((xrayData['E^2 dN/dE']/np.power(xrayData['E'], 2)) *
                           (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['flux_error_lo'] = ((xrayData['E^2 dN/dE error down']/np.power(xrayData['E'], 2)) *
                                    (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['flux_error_hi'] = ((xrayData['E^2 dN/dE error up']/np.power(xrayData['E'], 2)) *
                                    (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['ul'] = 0*xrayData['E'].astype(np.int)

            self.writeEcsvTable('sed/forNaima/{0}/xrt_{0}.ecsv'.format(sourceNowShort), sed)

        ############################################################################################
        # Swift-UVOT
        ############################################################################################

        if len(optData) > 0:
            sed = dict()
            sed['energy'] = (optData['E']*u.GeV).to('eV')
            sed['energy_error'] = (((optData['dE down'] + optData['dE up'])/2.) *
                                   u.GeV).to('eV')
            sed['flux'] = ((optData['E^2 dN/dE']/np.power(optData['E'], 2)) *
                           (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['flux_error_lo'] = ((optData['E^2 dN/dE error down']/np.power(optData['E'], 2)) *
                                    (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['flux_error_hi'] = ((optData['E^2 dN/dE error up']/np.power(optData['E'], 2)) *
                                    (1./(u.GeV*u.cm*u.cm*u.s))).to(1./(u.eV*u.cm*u.cm*u.s))
            sed['ul'] = 0*optData['E'].astype(np.int)

            self.writeEcsvTable('sed/forNaima/{0}/uvot_{0}.ecsv'.format(sourceNowShort), sed)

        return

    def plotSED(self, sedVeritas, sedFermi, sedXrt, sedUvot, eblDict,
                sourceNow, configAnal):

        sourceNowShort = sourceNow.split('+')[0].replace(' ', '')
        redshiftNow = configAnal['redshift'][sourceNowShort]
        if isinstance(redshiftNow, list):
            redshiftErrDn = redshiftNow[1]
            redshiftErrUp = redshiftNow[2]
            redshiftNow = redshiftNow[0]
        else:
            redshiftErrDn = 0
            redshiftErrUp = 0

        vFluxes = sedVeritas[sedVeritas['dN/dE error'] > 0]
        vUpperLimits = sedVeritas[sedVeritas['dN/dE error'] == 0]

        xTitle = 'E [eV]'
        titleNow = '{}, z={}'.format(sourceNow, redshiftNow)
        titleOffset = 1.2
        if '0033' in sourceNow and len(sedXrt) < 1 and len(sedUvot) < 1:
            xTitle = 'E [GeV]'
            titleNow = sourceNow
            titleOffset = 1.02
        # yTitle = '{}'.format(xTitle.split('[')[1].split(']')[0])
        # yTitle = r'$E^{2}$dN/d$E$ [{xTitle} cm$^{-2}$ s$^{-1}$]'.format(
        #                                                             xTitle.split('[')[1].split(']')[0])
        yTitle = r'$E^{2}$dN/d$E$ ['
        yTitle += xTitle.split('[')[1].split(']')[0]
        yTitle += r' cm$^{-2}$ s$^{-1}$]'

        leg = dict()
        leg['veritas'] = list()
        leg['fermi'] = list()
        leg['paiano'] = list()
        ax = host_subplot(111)

        tevToEv, gevToEv, mevToEv, kevToEv = 1e12, 1e9, 1e6, 1e3
        mevToTeV = 1e-6

        fVeritas = tevToEv
        fFermi = mevToEv
        if '0033' in sourceNow and len(sedXrt) < 1 and len(sedUvot) < 1:
            fVeritas = 1e3
            fFermi = 1e-3

        self.logStdout.info([['wb', 'Plotting SED for '],
                             ['p', '{}'.format(sourceNow)],
                             ['bb', '(z = {})'.format(redshiftNow)]])

        lH = ax.errorbar(fVeritas*vFluxes['E'],
                         fVeritas*np.power(vFluxes['E'], 2)*vFluxes['dN/dE'],
                         xerr=fVeritas*vFluxes['dE'],
                         yerr=[fVeritas*np.power(vFluxes['E'], 2) *
                               vFluxes['dN/dE error down'].astype(np.float),
                               fVeritas*np.power(vFluxes['E'], 2) *
                               vFluxes['dN/dE error up'].astype(np.float)],
                         markerfacecolor=self.colors[0], fmt=self.markers[0],
                         markersize=self.markersize, color=self.colors[0],
                         elinewidth=self.elinewidth, capsize=self.capsize,
                         label='VERITAS')
        ax.errorbar(fVeritas*vUpperLimits['E'],
                    fVeritas*np.power(vUpperLimits['E'], 2)*vUpperLimits['dN/dE'],
                    yerr=0.3*fVeritas*np.power(vUpperLimits['E'], 2)*vUpperLimits['dN/dE'],
                    fmt='_', color=self.colors[0], markersize=self.markersize,
                    elinewidth=self.elinewidth, markeredgewidth=2, capsize=3, uplims=True)

        leg['veritas'].append(lH)

        # Deabsorb VERITAS
        eblAbsorb = interpolate.interp1d(eblDict[redshiftNow]['E [TeV]'],
                                         eblDict[redshiftNow]['e^tau'], kind='linear',
                                         bounds_error=False, fill_value=1.)
        # Deabsorb with Paiano redshift
        eblAbsorbPaiano = interpolate.interp1d(eblDict[0.467]['E [TeV]'],
                                               eblDict[0.467]['e^tau'], kind='linear',
                                               bounds_error=False, fill_value=1.)
        if '0033' in sourceNow and len(sedXrt) < 1 and len(sedUvot) < 1:
            # The uncertainties include both statistical and model uncertaintiy.
            # The round is necessary to avoid floating point issues.
            zUp = round(redshiftNow + redshiftErrUp, 3)
            zDn = round(redshiftNow - redshiftErrDn, 3)
            eblAbsorbUncertUp = interpolate.interp1d(eblDict[zUp]['E [TeV]'],
                                                     eblDict[zUp]['e^tau'],
                                                     kind='linear', bounds_error=False,
                                                     fill_value=1.)
            eblAbsorbUncertDn = interpolate.interp1d(eblDict[zDn]['E [TeV]'],
                                                     eblDict[zDn]['e^tau'],
                                                     kind='linear', bounds_error=False,
                                                     fill_value=1.)

        absorbFactor = eblAbsorb(vFluxes['E'])
        # Remove points with crazy fluxes due to deabsorption
        maskNow = (fVeritas*np.power(vFluxes['E'], 2) *
                   absorbFactor*vFluxes['dN/dE'] < 100)
        vFluxes = vFluxes[maskNow]
        self.deabsorbFactorVeritas = absorbFactor[maskNow]

        absorbFactorPaiano = eblAbsorbPaiano(vFluxes['E'])

        lH = ax.errorbar(fVeritas*vFluxes['E'],
                         fVeritas*np.power(vFluxes['E'], 2) * self.deabsorbFactorVeritas *
                         vFluxes['dN/dE'],
                         xerr=fVeritas*vFluxes['dE'],
                         yerr=[fVeritas*np.power(vFluxes['E'], 2) * self.deabsorbFactorVeritas *
                               vFluxes['dN/dE error down'].astype(np.float),
                               fVeritas*np.power(vFluxes['E'], 2) * self.deabsorbFactorVeritas *
                               vFluxes['dN/dE error up'].astype(np.float)],
                         markerfacecolor='none', fmt=self.markers[0], markersize=self.markersize,
                         color=self.colors[0], elinewidth=self.elinewidth, capsize=self.capsize,
                         label='VERITAS (deabsorbed)')

        leg['veritas'].append(lH)

        # Get the full array back (look later how to reset the mask)
        vFluxes = sedVeritas[sedVeritas['dN/dE error'] > 0]

        if '0033' in sourceNow and len(sedXrt) < 1 and len(sedUvot) < 1:
            absorbFactorDn = eblAbsorbUncertDn(vFluxes['E'])
            bracketZdownFluxes = vFluxes[maskNow]
            self.deabsorbFactorVeritasDn = absorbFactorDn[maskNow]

            absorbFactorUp = eblAbsorbUncertUp(vFluxes['E'])
            bracketZupFluxes = vFluxes[maskNow]
            self.deabsorbFactorVeritasUp = absorbFactorUp[maskNow]

            lH = plt.fill_between(fVeritas*bracketZdownFluxes['E'],
                                  fVeritas*np.power(bracketZdownFluxes['E'], 2) *
                                  self.deabsorbFactorVeritasDn*bracketZdownFluxes['dN/dE'],
                                  fVeritas*np.power(bracketZupFluxes['E'], 2) *
                                  self.deabsorbFactorVeritasUp*bracketZupFluxes['dN/dE'],
                                  alpha=0.3, edgecolor=self.colors[0],
                                  facecolor=self.colors[0], zorder=1,
                                  label='Deabsorbtion uncertainty')

            leg['deabsorb'] = lH

            lH = ax.errorbar(fVeritas*vFluxes['E'],
                             fVeritas*np.power(vFluxes['E'], 2) * absorbFactorPaiano *
                             vFluxes['dN/dE'],
                             xerr=fVeritas*vFluxes['dE'],
                             yerr=[fVeritas*np.power(vFluxes['E'], 2) * absorbFactorPaiano *
                                   vFluxes['dN/dE error down'].astype(np.float),
                                   fVeritas*np.power(vFluxes['E'], 2) * absorbFactorPaiano *
                                   vFluxes['dN/dE error up'].astype(np.float)],
                             markerfacecolor='none', fmt=self.markers[0],
                             markersize=self.markersize, color='grey', alpha=0.5,
                             elinewidth=self.elinewidth, capsize=self.capsize,
                             label='VERITAS (deabsorbed z=0.467)', zorder=0)

            leg['paiano'].append(lH)

        # Plot Fermi data
        lH = ax.errorbar(fFermi*sedFermi['e_ctr'],
                         fFermi*sedFermi['e2dnde'],
                         yerr=[fFermi*sedFermi['e2dnde_err_lo'],
                               fFermi*sedFermi['e2dnde_err_hi']],
                         markerfacecolor=self.colors[1],
                         fmt=self.markers[1], markersize=self.markersize,
                         color=self.colors[1], elinewidth=self.elinewidth, capsize=self.capsize,
                         label='Fermi/LAT')

        leg['fermi'].append(lH)

        # Deabsorb Fermi/LAT
        self.deabsorbFactorFermi = eblAbsorb(mevToTeV*sedFermi['e_ctr'])
        # don't plot deabsorbed points where EBL effect is negligible
        maskNow = (self.deabsorbFactorFermi > 1.1)
        sedFermiMasked = sedFermi[maskNow]
        absorbFactor = self.deabsorbFactorFermi[maskNow]
        absorbFactorFermiPaiano = eblAbsorbPaiano(mevToTeV*sedFermiMasked['e_ctr'])

        if '0033' in sourceNow and len(sedXrt) < 1 and len(sedUvot) < 1:
            self.deabsorbFactorFermiDn = eblAbsorbUncertDn(mevToTeV*sedFermiMasked['e_ctr'])
            bracketZdownFluxes = self.deabsorbFactorFermiDn*sedFermiMasked['e2dnde']

            self.deabsorbFactorFermiUp = eblAbsorbUncertUp(mevToTeV*sedFermiMasked['e_ctr'])
            bracketZupFluxes = self.deabsorbFactorFermiUp*sedFermiMasked['e2dnde']

            plt.fill_between(fFermi*sedFermiMasked['e_ctr'],
                             fFermi*self.deabsorbFactorFermiDn*sedFermiMasked['e2dnde'],
                             fFermi*self.deabsorbFactorFermiUp*sedFermiMasked['e2dnde'],
                             alpha=0.3, edgecolor=self.colors[1],
                             facecolor=self.colors[1], zorder=2)

            lH = ax.errorbar(fFermi*sedFermiMasked['e_ctr'],
                             fFermi*absorbFactor*sedFermiMasked['e2dnde'],
                             yerr=[fFermi*absorbFactor*sedFermiMasked['e2dnde_err_lo'],
                                   fFermi*absorbFactor*sedFermiMasked['e2dnde_err_hi']],
                             markerfacecolor='none', fmt=self.markers[1],
                             markersize=self.markersize, color=self.colors[1],
                             elinewidth=self.elinewidth, capsize=self.capsize,
                             label='Fermi/LAT (deabsorbed)')

            leg['fermi'].append(lH)

            lH = ax.errorbar(fFermi*sedFermiMasked['e_ctr'],
                             fFermi*absorbFactorFermiPaiano*sedFermiMasked['e2dnde'],
                             yerr=[fFermi*absorbFactorFermiPaiano*sedFermiMasked['e2dnde_err_lo'],
                                   fFermi*absorbFactorFermiPaiano*sedFermiMasked['e2dnde_err_hi']],
                             markerfacecolor='none', fmt=self.markers[1],
                             markersize=self.markersize, color='grey', alpha=0.5,
                             elinewidth=self.elinewidth, capsize=self.capsize,
                             label='Fermi/LAT (deabsorbed z=0.467)')

            leg['paiano'].append(lH)

        # Plot Swift data
        if len(sedXrt) > 1:
            lH = ax.errorbar(gevToEv*sedXrt['E'],
                             gevToEv*sedXrt['E^2 dN/dE'],
                             xerr=[gevToEv*sedXrt['dE down'], gevToEv*sedXrt['dE up']],
                             yerr=[gevToEv*sedXrt['E^2 dN/dE error down'],
                                   gevToEv*sedXrt['E^2 dN/dE error up']],
                             markerfacecolor='none', fmt=self.markers[2],
                             markersize=self.markersize, color=self.colors[2],
                             elinewidth=self.elinewidth, capsize=self.capsize,
                             label='Swift/XRT')

            leg['xrt'] = lH

        if len(sedUvot) > 1:
            lH = ax.errorbar(gevToEv*sedUvot['E'],
                             gevToEv*sedUvot['E^2 dN/dE'],
                             xerr=[gevToEv*sedUvot['dE down'], gevToEv*sedUvot['dE up']],
                             yerr=[gevToEv*sedUvot['E^2 dN/dE error down'],
                                   gevToEv*sedUvot['E^2 dN/dE error up']],
                             markerfacecolor='none', fmt=self.markers[3],
                             markersize=self.markersize, color=self.colors[3],
                             elinewidth=self.elinewidth, capsize=self.capsize,
                             label='Swift/UVOT')

            leg['uvot'] = lH

            evTicks = [1e1, 1e3, 1e5, 1e7, 1e9, 1e11]
        else:
            evTicks = [1e3, 1e5, 1e7, 1e9, 1e11]

        plt.xlabel(xTitle, fontsize=self.fontsize)
        plt.ylabel(yTitle, fontsize=self.fontsize)
        plt.tick_params(axis='both', which='major', labelsize=self.fontsize)
        plt.title(titleNow, fontsize=15, y=titleOffset)

        if len(sedXrt) > 1 or len(sedUvot) > 1:
            freqTicks = [str(np.round(np.log10(v*(con.e/con.h)))) for v in evTicks]
            freqAx = ax.twin()  # freqAx is responsible for top axis
            freqAx.set_xticks(evTicks)
            freqAx.set_xticklabels(freqTicks, fontsize=self.fontsize)
            freqAx.set_xlabel(r'log($\nu$) [Hz]', fontsize=self.fontsize, labelpad=10)

            freqAx.axis['right'].major_ticklabels.set_visible(False)
            freqAx.axis['right'].major_ticks.set_visible(False)
            freqAx.axis['top'].major_ticklabels.set_visible(True)

        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

        legHandlers = [tuple(leg['veritas']), tuple(leg['fermi'])]
        legTitles = ['VERITAS (deabsorbed, z=0.33)', 'Fermi/LAT (deabsorbed, z=0.33)']
        if 'deabsorb' in leg:
            legHandlers.append(leg['deabsorb'])
            legTitles.append('Deabs. uncertainty')
        if len(leg['paiano']) > 0:
            legHandlers.append(tuple(leg['paiano']))
            legTitles.append('Deabsorbed with z=0.467')
        if 'xrt' in leg:
            legHandlers.append(leg['xrt'])
            legTitles.append('Swift/XRT')
        if 'uvot' in leg:
            legHandlers.append(leg['uvot'])
            legTitles.append('Swift/UVOT')
        legend = plt.legend(legHandlers, legTitles,
                            handler_map={tuple: HandlerTuple(ndivide=None)},
                            fontsize='large')

        plt.grid(True)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()
        plt.savefig('sed/' + sourceNowShort + '_sed.pdf')

        plt.clf()

        # Save tables of the SED in an easy to read format.
        self.saveSEDs(vFluxes, sedFermi, sourceNowShort)
        # Save tables of the SED in ECSV format for SED modelling.
        self.saveEcsvSED(sourceNowShort, vFluxes, sedFermi, sedXrt, sedUvot)

        return


if __name__ == '__main__':

    sed = sed()
    logStdout = sed.getLogger()

    parser = argparse.ArgumentParser(description=('Plot SED.'))
    parser.add_argument('-c', '--copy', action='store_true', default=False,
                        help='Copy new log files containing the VERITAS SED or plot the spectra')
    parser.add_argument('-e', '--ebl', action='store_true', default=False,
                        help='Plot only HE and VHE with the EBL uncertainty (only for 1ES 0033)')
    args = parser.parse_args()

    configAnalFile = 'configAnalysis.yaml'
    configAnal = gUtil.readYamlFile(logStdout, configAnalFile)

    eblDict = sUtil.readFranceschini(logStdout, configAnal['EBL']['franceschini2008']['file'])

    prefixLog = '/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/spectra/'
    suffixLog = '/log/runSpectrum.sh.o'

    for sourceNowShort, sourceNow in configAnal['sources'].items():

        logStdout.info([['wb', 'Running'],
                        ['p', ' {}'.format(sourceNow)]])

        if args.copy:
            listOfFiles = glob.glob('makeSpectrum/{}/log/runSpectrum.sh.o*'.format(sourceNowShort))
            latestFile = max(listOfFiles, key=os.path.getctime)
            logStdout.info([['wb', 'Copying'],
                            ['y', latestFile],
                            ['wb', 'to'],
                            ['y', 'sed/veritasLogs/{}.log'.format(sourceNowShort)]])
            subprocess.call('/bin/cp {0} sed/veritasLogs/{1}.log'.format(latestFile,
                                                                         sourceNowShort),
                            shell=True)
        else:
            # Check whether a newer VERITAS spectrum (log file) exists
            # and issue a warning if it does
            veritasFileName = 'sed/veritasLogs/{}.log'.format(sourceNowShort)
            listOfFiles = glob.glob('makeSpectrum/{}/log/runSpectrum.sh.o*'.format(sourceNowShort))
            latestFile = max(listOfFiles, key=os.path.getctime)
            if not filecmp.cmp(latestFile, veritasFileName):
                logStdout.warn([['wr', 'NOTICE!'],
                                ['wb', 'There is a newer spectrum file:'],
                                ['y', '\n{}'.format(latestFile)],
                                ['wb', '\nare you sure you do not want to run'],
                                ['g', 'python sed.py -c'],
                                ['wb', 'first?']])

            veritasData = sUtil.readVeritasData(logStdout, veritasFileName, fromLog=True)

            fermiSourceName = gUtil.getSourceNameFermi(logStdout, sourceNow)
            fermiFileName = os.path.join(configAnal['fermi']['baseDir'],
                                         'yearlyLightCurves',
                                         sourceNowShort,
                                         sourceNowShort,
                                         fermiSourceName + '_sed.npy')

            fermiData = sUtil.readFermiData(logStdout, fermiFileName)

            if configAnal['swift']['swiftSource'][sourceNowShort] == 'mireia':
                swiftFile = 'swift/swiftFromMireia/{}.yml'.format(sourceNowShort)
            else:
                swiftFile = 'swift/onlineTool/{}/spectrum/xrtSed{}.dat'.format(
                            sourceNowShort,
                            configAnal['swift']['xrtMethod'][sourceNowShort])

            (xrtData,
             uvotData) = sUtil.readSwiftData(logStdout,
                                             swiftFile,
                                             configAnal['swift']['swiftSource'][sourceNowShort])

            if args.ebl and '0033' in sourceNow:
                sed.plotSED(veritasData, fermiData, list(), list(), eblDict,
                            sourceNow, configAnal)
            else:
                sed.plotSED(veritasData, fermiData, xrtData, uvotData, eblDict,
                            sourceNow, configAnal)

    if not args.copy:

        # Merge DCF plots into one PDF
        fileToMerge = list()
        for sourceNowShort in configAnal['sources'].keys():
            fileToMerge.append('sed/{}_sed.pdf'.format(sourceNowShort))
        gUtil.mergePDFs(logStdout, fileToMerge, 'sed/seds.pdf')
        fileToMerge.clear()
