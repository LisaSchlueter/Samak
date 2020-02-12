#!/bin/bash
# copy hdf5 files from christians homefolder on the tristan server to my local h5 folder
# input argument: 1) login 2) folder which you want to copy
GetH5FilesServerFT(){

FTdir = FirstTritium.katrin
# synchronize hdf5 data folder with FT data folder                                                                                                                                             
rsync -i -a --ignore-existing -v $1@pcltr-01.mpp.mpg.de:/home/iwsatlas1/karlch/RunFiles/$FTdir/H5Files/ ./hdf5/$FTdir/

#scp -r $1@pcltr-01.mpp.mpg.de:/home/iwsatlas1/karlch/RunFiles/$2/H5Files/  ./hdf5/$2/
#mv ./hdf5/H5Files ./hdf5/$2
}
GetH5FilesServerFT $1

