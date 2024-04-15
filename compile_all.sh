#!/bin/bash
if (( $# != 1 )); then
    >&2 echo "Usage: $0 <buildtype:Release/Debug/Profile> "
    buildtype=Release
fi

buildtype=$1 # Release or Debug 

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

# compile mmwis and hils

rm -rf deploy
# rm -rf build
mkdir -p build
cd build

if [[ "$buildtype" == "Profile" ]]; then
  cmake ../ -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
else
    cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_BUILD_TYPE=${buildtype}
fi
make -j $NCORES
cd ..

mkdir -p deploy
cp -f ./build/solver/branch_reduce/branch_reduce_convergence    deploy/branch_reduce
cp -f ./build/solver/branch_reduce/kernelization                deploy/kernelization


