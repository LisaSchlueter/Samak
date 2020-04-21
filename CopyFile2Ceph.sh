#!/bin/bash 

# copy folder from local machine to ceph storage server
# to be run from your local mashine                                                                             
# input argument: 1) MPP login 2) folder which you want to copy                                                      
# start path after Samak3.0 e.g. tritium-data/hdf5/Knm1/

CopyFile2Server(){                                                                                          
ssh $1@csltr-02.mpp.mpg.de "mkdir -p /remote/ceph/user/s/schluete/Samak/$2/"
scp -r ./$2/* $1@csltr-02.mpp.mpg.de:/remote/ceph/user/s/schluete/Samak/$2/
}
CopyFile2Server $1 $2



