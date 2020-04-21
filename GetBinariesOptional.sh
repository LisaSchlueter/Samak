#!/bin/bash
# get optinal Samak binaries: Fit results, twins, etc.
# input argument: 1) MPP login 
GetBinariesOptional(){
#---------------------------
# create directories if necessary
mkdir -p ./tritium-data/mat/TwinKnm1/
mkdir -p ./tritium-data/mat/TwinKnm2/
mkdir -p ./tritium-data/fit/Knm1/Uniform/
mkdir -p ./tritium-data/fit/Knm2/Uniform/
mkdir -p ./tritium-data/fit/Knm1/SingleRingFit/	
mkdir -p ./tritium-data/fit/Knm2/SingleRingFit/	
mkdir -p ./tritium-data/sensitivity/Knm1/
mkdir -p ./tritium-data/sensitivity/Knm2/
mkdir -p ./inputs/CovMat/RF/CM/
mkdir -p ./inputs/CovMat/RF/CMPartVar/
mkdir -p ./inputs/CovMat/LongPlasma/CM/
mkdir -p ./inputs/CovMat/TASR/CM/
mkdir -p ./inputs/CovMat/FSD/CM/

#---------------- start synchronization: 

# fit results
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/fit/ ./tritium-data/fit/

# sensitivites
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/sensitivity/ ./tritium-data/sensitivity/

# twins
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm1/*.mat ./tritium-data/mat/TwinKnm1/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm2/*.mat ./tritium-data/mat/TwinKnm2/

# covariance matrices
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/CovMat/RF/CM/*.mat ./inputs/CovMat/RF/CM/ 
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/CovMat/RF/CMPartVar/*.mat ./inputs/CovMat/RF/CMPartVar/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/CovMat/LongPlasma/CM/*.mat ./inputs/CovMat/LongPlasma/CM/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/CovMat/FSD/CM/*.mat ./inputs/CovMat/FSD/CM/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/CovMat/TASR/CM/*.mat ./inputs/CovMat/TASR/CM/

}
GetBinariesOptional $1
