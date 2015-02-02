#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
matlab -nodisplay -nojvm < scriptBatchJuraGgwhiteSpgmpCd1KL.m
matlab -nodisplay -nojvm < scriptBatchJuraGgwhiteSpgmpCd2KL.m
exit 0




