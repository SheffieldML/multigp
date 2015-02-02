#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < demSpmgpGgToy3KL.m
matlab -nodisplay -nojvm < demSpmgpGgToy4KL.m 
matlab -nodisplay -nojvm < demSpmgpGgToy5KL.m 
exit 0
