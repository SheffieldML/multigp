#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < compilerGgwhiteDTCVAR.m

exit 0