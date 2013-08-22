#!/bin/bash

#Real job
jobsub -g --X509_USER_PROXY /scratch/zpli/grid/zpli.lbne.proxy -N 1 -q text.sh `whoami` `pwd`

#Test job
#jobsub -g --X509_USER_PROXY /scratch/zpli/grid/zpli.lbne.proxy  -N 1000 -q optical_test.sh `whoami` `pwd`
