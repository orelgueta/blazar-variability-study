#!/usr/bin/python

import subprocess
import sys
import os
import argparse
from lib import colourLogger
from lib import generalUtil as gUtil


def submitJobs(logStdout, sourceNow, mode, configAnalFile, configAnal):

    if mode.lower() == 'lc':
        suffix = 'LC'
    elif mode.lower() == 'spectrum' or mode.lower() == 'spectra':
        suffix = 'Spectrum'
    else:  # Shouldn't get here anyway
        logStdout.critical([['r', 'Do not know the running mode:'],
                            ['bb', mode]])
        sys.exit(1)

    workDir = '{}/make{}'.format(os.getcwd(), suffix)
    logDir = os.path.join(workDir, sourceNow, 'log')
    try:
        subprocess.check_call(['mkdir', '-p', logDir])
    except subprocess.CalledProcessError as e:
        logStdout.critical([['r', 'Could not create log directory'],
                            ['y', logDir]])
        sys.exit(1)

    outLog = '{}/out.log'.format(logDir)
    errLog = '{}/err.log'.format(logDir)
    if configAnal['batchFarm']:
        try:
            logStdout.info([['wb', 'Sending'],
                            ['p', sourceNow],
                            ['wb', 'to the qsub batch farm to produce a'],
                            ['g', suffix]])
            with open(outLog, 'w') as out, open(errLog, 'w') as err:
                subprocess.check_call(['qsub', '-js', '9', '-l', 'h_cpu=00:08:00',
                                       '-l', 'tmpdir_size=40G', '-l', 'h_rss=4G',
                                       '-V', '-o', logDir, '-e', logDir,
                                       '{}/run{}.sh'.format(workDir, suffix), workDir, sourceNow,
                                       '{}/{}'.format(os.getcwd(), configAnalFile)],
                                      cwd=workDir, stdout=out,
                                      stderr=err)
        except subprocess.CalledProcessError as e:
            logStdout.critical([['r', 'Could not submit batch job for'],
                                ['p', sourceNow],
                                ['wb', 'Reported error is\n'],
                                ['rb', e]])
            sys.exit(1)

    else:
        logStdout.info([['wb', 'Running'],
                        ['p', sourceNow],
                        ['wb', 'locally to produce a'],
                        ['g', suffix]])
        logStdout.info([['wb', 'Patience, output is printed to\n'],
                        ['y', outLog]])
        logStdout.info([['wb', 'Error is printed to\n'],
                        ['y', errLog]])
        try:
            with open(outLog, 'w') as out, open(errLog, 'w') as err:
                subprocess.check_call(['python', '{}/calc{}.py'.format(workDir, suffix),
                                       sourceNow, '{}/{}'.format(os.getcwd(), configAnalFile)],
                                      cwd=workDir, stdout=out,
                                      stderr=err)
        except subprocess.CalledProcessError as e:
            logStdout.critical([['r', 'Could not run'],
                                ['p', sourceNow],
                                ['wb', 'Reported error is\n'],
                                ['rb', e]])
            sys.exit(1)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Which VERITAS output to generate?')
    parser.add_argument('mode', action='store',
                        choices=['lc', 'LC', 'spectrum', 'spectra'],
                        help='Produce lightcurve or spectra')

    args = parser.parse_args()

    logStdout = colourLogger.initStdOutLogger()

    configAnalFile = 'configAnalysis.yaml'
    configAnal = gUtil.readYamlFile(logStdout, configAnalFile)

    for sourceNowShort in configAnal['sources'].keys():
        submitJobs(logStdout, sourceNowShort, args.mode, configAnalFile, configAnal)
