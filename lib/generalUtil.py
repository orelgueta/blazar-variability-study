#!/usr/bin/python

import os
import sys
import inspect
import yaml
import numpy as np
from astropy.io import fits
import astropy.units as u
from lib import colourLogger
from lib import PyPDF2


def checkFileExists(logger, fileName, whoCalledMe, exitIfNot=True):

    if not os.path.isfile(fileName):
        if exitIfNot:
            logger.critical([['c', 'From {}:'.format(whoCalledMe)],
                             ['rb', fileName + ' is not a file!']])
            sys.exit(1)
        else:
            logger.warn([['c', 'From {}:'.format(whoCalledMe)],
                         ['r', fileName + ' is not a file, continuing anyway']])
            return False
    else:
        return True


def readYamlFile(logger, fileName):

    checkFileExists(logger, fileName, inspect.stack()[0][3])
    logger.info([['wb', 'Reading'],
                 ['y', ' {}'.format(fileName)]])

    with open(fileName, 'r') as stream:
        try:
            outDict = yaml.load(stream, Loader=yaml.SafeLoader)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)

    return outDict


def getSourceNameFermi(logger, source):

    hdulist = fits.open('fermiCatalog/gll_psc_8year_v5.fit')
    sourceCatalog = hdulist[1].data

    sourceNameFermi = ''

    # Find the source name in the Fermi catalog.
    # This should have been done with a simple "index()" command,
    # but it doesn't work for some reason
    for i_src, sourceFermi in enumerate(sourceCatalog['ASSOC_TEV']):
        if source in sourceFermi:
            sourceNameFermi = sourceCatalog['Source_Name'][i_src].replace(' ', '_').lower()

    if sourceNameFermi == '':
        logger.critical([['rb', 'Could not find '],
                         ['p', source],
                         ['rb', ' in the Fermi catalog']])
        sys.exit(1)
    else:
        logger.info([['wb', 'Found '],
                     ['p', source],
                     ['wb', ' in the Fermi catalog, '],
                     ['p', sourceNameFermi]])

    return sourceNameFermi


def mergePDFs(logger, pdfs, outputName):

    logger.info([['wb', 'Will attempt to merge ['],
                 ['y', ', '.join(pdfs)],
                 ['wb', '] into'],
                 ['y', outputName]])

    merger = PyPDF2.PdfFileMerger()

    for pdfNow in pdfs:
        checkFileExists(logger, pdfNow, inspect.stack()[0][3])
        merger.append(pdfNow)

    merger.write(outputName)


def ecpl(e, pars):
    """
    One dimensional power law with an exponential cutof model function.
    The energy to evaluate is e and pars is a dictionary containing the parameters:
    'reference' - the reference/normalisation energy
    'amplitude' - the flux at the reference energy
    'index' - the power law index
    'e_c' - the cutoff energy
    """

    xx = e / pars['reference']
    return pars['amplitude'] * xx**(-pars['index']) * np.exp(-(e / pars['e_c']))


def crabFluxAbove(logger, e):

    preCalc = readYamlFile(logger, '../lib/crabFlux.yml')

    if e in preCalc['fluxAboveEnergy']:
        return preCalc['fluxAboveEnergy'][e]
    else:
        hess_ecpl = {'amplitude': 3.76e-11 * u.Unit('1 / (cm2 s TeV)'),
                     'index': 2.39,
                     'e_c': (14.3 * u.TeV),
                     'reference': 1 * u.TeV}

        dE = 0.001
        energies = np.arange(e, 30, dE)
        integral = 0
        for eNow in energies:
            integral += ecpl(eNow*u.TeV, hess_ecpl)*dE*u.TeV

        return integral


def getTableTex(table, captionTex, **additionalPars):
    """Produce a Latex table from the structured array in "table"
    and the caption text in "captionTex". It is assumed that the headers
    of the table are given in the names of the structured array.
    To avoid long tables, print only up to 30 lines.
    """

    columnsTex = ' & '.join(table.dtype.names) + ' \\\\'
    boxSize = ''
    if len(columnsTex) > 90:
        boxSize = r'\resizebox{\linewidth}{!}'

    if additionalPars['firstColumnLeftAligned']:
        nCols = 'l ' + 'c '*(len(table.dtype.names) - 1)
    else:
        nCols = 'c '*len(table.dtype.names)
    tableTex = """%
                  \\begin{{table}}[!htbp]
                      \\begin{{center}}
                      {boxSize:s}{{
                          \\begin{{tabular}}{{{nCols:s}}}
                              \\toprule
                              {colHead:s}
                              \\midrule""".format(nCols=nCols,
                                                  colHead=columnsTex, boxSize=boxSize)

    if 'nLines' in additionalPars:
        nLines = additionalPars['nLines']
    else:
        nLines = 30 if table.shape[0] > 30 else table.shape[0]
    for i_line in range(nLines):
        for i_entry, entry in enumerate(table[i_line]):
            if isinstance(entry, float):
                if 'scientific' in additionalPars:
                    if 'precision' in additionalPars:
                        table[i_line][i_entry] = ('{0:.{1:d}e}'.
                                                  format(entry,
                                                         additionalPars['precision']))
                    else:
                        table[i_line][i_entry] = '{0:e}'.format(entry)
                else:
                    if 'precision' in additionalPars:
                        table[i_line][i_entry] = ('{0:.{1:d}f}'.
                                                  format(entry,
                                                         additionalPars['precision']))
                    else:
                        if abs(entry) < 1e-10 and entry != 0:
                            table[i_line][i_entry] = 0
        datTex = ' & '.join(str(line) for line in table[i_line]) + ' \\\\'
        tableTex += """
        {datTex:s}""".format(datTex=datTex)

    tableTex += """
                \\bottomrule
                        \\end{{tabular}}
                    }}
                        \\caption{{{captionTex:s}}}
                    \\end{{center}}
                \\end{{table}}
                %""".format(captionTex=captionTex)

    return tableTex
