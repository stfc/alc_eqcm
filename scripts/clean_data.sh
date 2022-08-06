#!/usr/bin/env bash

[ -d 'ANALYSIS_EQCM' ]    && rm -R ANALYSIS_EQCM
[ -d 'ATOMISTIC_MODELS' ] && rm -R ATOMISTIC_MODELS 
[ -d 'RESTART' ]          && rm -R RESTART
[ -f 'OUT_EQCM' ]         && rm -R OUT_EQCM
