#!/usr/bin/env bash

folder="build-intel"
rm -rf $folder && mkdir $folder && cd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=Off
make
