#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < scriptBatchJuraGgFullCd.m
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd1.m
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd2.m
matlab -nodisplay -nojvm < scriptBatchJuraGgSpmgpCd3.m
exit 0
