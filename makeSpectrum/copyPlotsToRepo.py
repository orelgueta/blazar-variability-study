#!/usr/bin/python

import subprocess
from lib import colourLogger
from lib import generalUtil as gUtil

if __name__ == '__main__':

    logStdout = colourLogger.initStdOutLogger()

    configAnal = gUtil.readYamlFile(logStdout, '../configAnalysis.yaml')

    for sourceNow in configAnal['sources'].keys():
        dest = ('/afs/ifh.de/group/cta/scratch/ogueta/vts/variabilityStudy/'
                'blazarvariability/plots/spectra/{0}/').format(sourceNow)
        logStdout.info([['wb', 'Copying spectrum plots of'],
                        ['p', sourceNow],
                        ['wb', 'to'],
                        ['y', dest]])
        subprocess.call('/bin/cp {0}/*.pdf {1}'.format(sourceNow, dest), shell=True)
