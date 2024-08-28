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
    if [[ "$buildtype" == "Release" ]]; then
        # cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_BUILD_TYPE=${buildtype} 
        cmake ../ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
    else
        cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_BUILD_TYPE=${buildtype} -DREDUCTION_INFO=ON
    fi
fi
make -j $NCORES
cd ..

mkdir -p deploy
cp -f ./build/solver/branch_reduce/branch_reduce_convergence    deploy/branch_reduce
cp -f ./build/solver/branch_reduce/kernelization                deploy/kernelization
cp -f ./build/solver/branch_reduce/generate_full_graph_training_data       deploy/generate_training_data
cp -f ./build/metis_to_cosmo                                    deploy/metis_to_cosmo
cp -f ./build/graphchecker                                      deploy/graphchecker


