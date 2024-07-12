#!/usr/bin/env bash

folder="test-intel"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON 
make
ctest --output-on-failure
