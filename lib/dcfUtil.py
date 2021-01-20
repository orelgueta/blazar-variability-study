#!/usr/bin/python

import os
import numpy as np
import glob
import inspect
from lib import generalUtil as gUtil


def saveDataForDCF(directory, filePrefix, time, flux, fluxError):

    reducedData = np.column_stack((time,
                                   flux,
                                   fluxError))
    header = 'Date (MJD)   Flux   Flux error'

    np.savetxt('{}/dcf/lightcurves/{}.txt'.format(directory, filePrefix),
               reducedData, delimiter=',   ', header=header)

    return


def dcfTwoLCs(directory, outputPrefix, lcFile1, lcFile2,
              minPeriod, maxPeriod, binWidth):

    outputLog = '{}/dcf/log/{}.txt'.format(directory, outputPrefix)
    dcfCmd = ('python pydcf/dcf.py -o -np {0} {1} {2} {3} {4} >& {5}'.format
              (lcFile1, lcFile2, minPeriod, maxPeriod, binWidth, outputLog))

    os.system(dcfCmd)

    dcfOutputFile = 'dcf_output.csv'
    targetLoc = '{}/dcf/results/{}.csv'.format(directory, outputPrefix)
    os.system('mv dcf_output.csv {}'.format(targetLoc))

    dcfOutput = np.loadtxt(targetLoc, comments='#', delimiter=',')

    return dcfOutput
