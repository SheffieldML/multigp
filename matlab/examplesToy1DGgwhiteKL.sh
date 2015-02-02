#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < demSpmgpNoiseToy1KL.m
matlab -nodisplay -nojvm < demSpmgpNoiseToy2KL.m

exit 0


