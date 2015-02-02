#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < scriptBatchJuraGgwhiteFullCd.m
matlab -nodisplay -nojvm < scriptBatchJuraGgwhiteSpgmpCd1.m
matlab -nodisplay -nojvm < scriptBatchJuraGgwhiteSpgmpCd2.m
exit 0
