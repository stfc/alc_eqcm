#!/usr/bin/env bash

folder="build-intel-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
