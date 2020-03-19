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

#---------------- start synchronization: 

# fit results
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/fit/ ./tritium-data/fit/

# sensitivites
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/sensitivity/ ./tritium-data/sensitivity/

# twins
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm1/*.mat ./tritium-data/mat/TwinKnm1/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm2/*.mat ./tritium-data/mat/TwinKnm2/

}
GetBinariesOptional $1
