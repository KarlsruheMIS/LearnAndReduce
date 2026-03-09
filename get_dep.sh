#!/bin/bash

git clone https://github.com/KennethLangedal/CHILS.git
cd CHILS 

git checkout b3f3dbeb38efe44ca7bedbd736df56a94d5e3a36 

make libCHILS.a 

mkdir -p ../extern/CHILS
mv -f libCHILS.a ../extern/CHILS
cp -f include/chils.h ../extern/CHILS

cd ..
rm -rf CHILS 