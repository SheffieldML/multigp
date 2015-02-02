#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd1KL.m
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd2KL.m
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd3KL.m
exit 0

