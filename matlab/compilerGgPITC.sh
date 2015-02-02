#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < compilerGgPITC.m

exit 0