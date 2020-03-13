#!/bin/bash 
# to be run from your local mashine
# input argument: 1) login 2) folder which you want to copy                                                                                                                                          
# start path after Samak3.0 e.g. tritium-data/hdf5/Knm1/
CopyFile2Server(){
scp -r ./$2/* $1@csltr-02.mpp.mpg.de:./Samak3.0/$2
}
CopyFile2Server $1 $2



