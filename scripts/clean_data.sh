#!/usr/bin/env bash

[ -d 'ANALYSIS_EQCM' ]    && rm -R ANALYSIS_EQCM
[ -d 'ATOMISTIC_MODELS' ] && rm -R ATOMISTIC_MODELS 
[ -d 'RESTART' ]          && rm -R RESTART
[ -f 'OUT_EQCM' ]         && rm OUT_EQCM
[ -f 'SET_SIMULATION' ]   && rm SET_SIMULATION
[ -f 'SET_KPOINTS'    ]   && rm SET_KPOINTS
[ -f 'POTCAR'         ]   && rm POTCAR
[ -f 'RECORD_MODELS'  ]   && rm RECORD_MODELS
