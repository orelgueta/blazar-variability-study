#!/bin/bash

targetDir='/Users/ogueta/work/veritas/variabilityStudy/ICRC2019/proceeding/figures'
sourceDir=('1ES1011_VERITAS_Fermi_nightly' 
           # '1ES1011_VERITAS_Swift_nightly'
           # '1ES0033_VERITAS_Swift_weekly'
           # '1ES0229_VERITAS_Swift_nightly'
           '1ES1218_VERITAS_Fermi_weekly'
           '1ES1218_VERITAS_Swift_nightly'
           # '1ES1218_swiftFermi_weekly'
           'PG1553_VERITAS_Fermi_nightly'
           'PG1553_VERITAS_Swift_nightly'
           # 'PG1553_tau-1_tau-2_nightly'
           # 'PG1553_tau-1_tau-3_nightly'
           '1ES1011_tau-1_tau-2_nightly'
           '1ES1011_tau-1_tau-3_nightly'
          )

for source in "${sourceDir[@]}"; do
    echo 'copying' ${source}/zdcf.pdf 'to' ${targetDir}/${source}_zdcf.pdf
    /bin/cp ${source}/zdcf.pdf ${targetDir}/${source}_zdcf.pdf
done
