#!/usr/bin/env bash

folder="test-ifort"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON 
make
ctest --output-on-failure
