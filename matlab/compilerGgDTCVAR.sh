#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < compilerGgDTCVAR.m

exit 0