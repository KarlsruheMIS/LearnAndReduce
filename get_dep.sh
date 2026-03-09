#!/bin/bash

git clone https://github.com/KarlsruheMIS/CHILS.git
cd CHILS 

git checkout e2391453739b5977c46fc6915109345b6b5eacd4 

make libCHILS.a 

mkdir -p ../extern/CHILS
mv -f libCHILS.a ../extern/CHILS
cp -f include/chils.h ../extern/CHILS

cd ..
rm -rf CHILS 