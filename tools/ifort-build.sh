#!/usr/bin/env bash

folder="build-ifort"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=Off
make
