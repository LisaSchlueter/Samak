#!/bin/bash
# copy hdf5 files from christians homefolder on the tristan server to my local h5 folder
# input argument: 1) MPP login 
GetBinaries(){
#---------------------------
# create directories if necessary
mkdir -p ./tritium-data/hdf5/Knm1/
mkdir -p ./tritium-data/hdf5/Knm2/
mkdir -p ./tritium-data/mat/TwinKnm1/
mkdir -p ./tritium-data/mat/TwinKnm2/
mkdir -p ./tritium-data/fit/Knm1/Uniform/
mkdir -p ./tritium-data/fit/Knm2/Uniform/
mkdir -p ./tritium-data/fit/Knm1/SingleRingFit/	
mkdir -p ./tritium-data/fit/Knm2/SingleRingFit/	
mkdir -p ./tritium-data/sensitivity/Knm1/
mkdir -p ./tritium-data/sensitivity/Knm2/
mkdir -p ./tritium-data/NPfactor/
mkdir -p ./inputs/FSD/
mkdir -p ./inputs/ELossFunction/
mkdir -p ./inputs/WGTSMACE/
mkdir -p ./fitting/fminuit/

#---------------- start synchronization: 

# synchronize hdf5 data folder with webtrium data folder
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/home/iwsatlas1/karlch/.local/share/webtrium/Knm1/RunFiles/H5/ ./tritium-data/hdf5/Knm1/
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/home/iwsatlas1/karlch/.local/share/webtrium/Knm2/RunFiles/H5/ ./tritium-data/hdf5/Knm2/

# inputs: FSD, energy loss function, wgts density profile
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/FSD/ ./inputs/FSD/
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/ELossFunction/ ./inputs/ELossFunction/  
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/WGTSMACE/wgts_density_profile_for_thierry.dat  ./inputs/WGTSMACE/wgts_density_profile_for_thierry.dat 

# fit results, standard sensitivities, non-poiss factor
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/fit/ ./tritium-data/fit/
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/sensitivity/ ./tritium-data/sensitivity/
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/NPfactor/*.mat ./tritium-data/NPfactor/

# twins
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm1/*.mat ./tritium-data/mat/TwinKnm1/
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/mat/TwinKnm2/*.mat ./tritium-data/mat/TwinKnm2/

# fminuit
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/fitting/fminuit/ ./fitting/fminuit/

# background ROI info
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/studies/ROICoverage_S_B/ ./studies/ROICoverage_S_B/

# FPD viewer frame
rsync -i -a -v $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/FPDFrame.fig ./inputs/

}
GetBinaries $1
