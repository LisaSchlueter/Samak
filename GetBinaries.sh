#!/bin/bash
# copy hdf5 files from christians homefolder on the tristan server to my local h5 folder
# input argument: 1) MPP login 
GetBinaries(){
#---------------------------
# create directories if necessary
mkdir -p ./tritium-data/hdf5/Knm1/
mkdir -p ./tritium-data/hdf5/Knm2/
mkdir -p ./tritium-data/NPfactor/
mkdir -p ./inputs/FSD/
mkdir -p ./inputs/ELossFunction/
mkdir -p ./inputs/WGTSMACE/
mkdir -p ./fitting/fminuit/
mkdir -p ./inputs/MonitorSpec/
mkdir -p ./inputs/WGTSMACE/WGTS_ISProb/
mkdir -p ./inputs/RunLists/KNM2/

#---------------- start synchronization: 

# synchronize hdf5 data folder with webtrium data folder
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/home/iwsatlas1/karlch/.local/share/webtrium/Knm1/RunFiles/H5/ ./tritium-data/hdf5/Knm1/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/home/iwsatlas1/karlch/.local/share/webtrium/Knm2/RunFiles/H5/ ./tritium-data/hdf5/Knm2/

# inputs: FSD, energy loss function, wgts density profile
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/FSD/ ./inputs/FSD/
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/ELossFunction/ ./inputs/ELossFunction/  
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/WGTSMACE/wgts_density_profile_for_thierry.dat  ./inputs/WGTSMACE/wgts_density_profile_for_thierry.dat 

# non-poiss factor
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/tritium-data/NPfactor/*.mat ./tritium-data/NPfactor/

# fminuit
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/fitting/fminuit/ ./fitting/fminuit/

# background ROI info
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/studies/ROICoverage_S_B/ ./studies/ROICoverage_S_B/

# FPD viewer frame
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/FPDFrame.fig ./inputs/

# Main spectrometer qU drift from Monitor spectrometer (KNM2)
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/MonitorSpec/* ./inputs/MonitorSpec/

# energy-dependent inel. scattering cross section: pre-calculated mesh grid
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/WGTSMACE/WGTS_ISProb/InitISProbMeshGrid* ./inputs/WGTSMACE/WGTS_ISProb/

# golden run list
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/inputs/RunLists/KNM2/ ./inputs/RunLists/KNM2/

}
GetBinaries $1
