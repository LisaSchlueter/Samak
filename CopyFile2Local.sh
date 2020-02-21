#!/bin/bash 
# to be run from your local mashine                                                                                                                                                                                                                                                                                      
# input argument: 1) login 2) folder which you want to copy from server to local machine: start after Samak2.0 e.g. tritium-data/hdf5/Knm1/                                                                                                                                          
CopyFile2Local(){
scp $1@pcltr-01.mpp.mpg.de:./Samak3.0/$2/* ./$2/
}
CopyFile2Local $1 $2



