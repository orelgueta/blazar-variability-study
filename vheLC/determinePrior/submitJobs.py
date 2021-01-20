#!/usr/bin/python

import subprocess
import sys
import os

if __name__ == '__main__':

    sources = {'1ES0033': '1ES 0033+595',
               '1ES0502': '1ES 0502+675',
               '1ES1011': '1ES 1011+496',
               '1ES1218': '1ES 1218+304',
               '1ES0229': '1ES 0229+200',
               'RGBJ0710': 'RGB J0710+591',
               'PG1553':  'PG 1553+113',
               'PKS1424': 'PKS 1424+240'
               }

    bins = [
            'yearly',
            'monthly',
            'weekly'
            ]

    workDir = os.getcwd()

    for dirNow, sourceNow in sources.items():

        for binning in bins:

            logDir = os.path.join(workDir, dirNow, 'log/')
            try:
                subprocess.check_call(['mkdir', '-p', logDir])
            except subprocess.CalledProcessError as e:
                print('Could not create log directory')
                sys.exit(1)

            try:
                with open('{}/out.log'.format(logDir), 'w') as out, open('{}/err.log'.format(logDir), 'w') as err:
                    subprocess.check_call(['qsub', '-js', '9', '-l', 'h_cpu=24:00:00',
                                           '-l', 'tmpdir_size=40G', '-l', 'h_rss=4G',
                                           '-V', '-o', logDir, '-e', logDir,
                                           '{}/runPriors.sh'.format(workDir), workDir,
                                           dirNow, binning],
                                          cwd=workDir, stdout=out,
                                          stderr=err)
            except subprocess.CalledProcessError as e:
                print('Could not submit jobs')
                sys.exit(1)
