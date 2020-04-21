#!/bin/bash 
# to be run from your local mashine
# input argument: 1) login 
# will copy automatically all necessary Samak files to server                                                                                                                  
CopyFile2Server(){
path1=inputs/ELossFunction
path2=inputs/FSD
#path3 = inputs/WGTSMACE
path4=tritium-data/mat/Knm1
scp -r ./$path1/* $1@pcltr-01.mpp.mpg.de:./Samak2.0/$path1/
scp -r ./$path2/* $1@pcltr-01.mpp.mpg.de:./Samak2.0/$path2/
scp -r ./$path4/* $1@pcltr-01.mpp.mpg.de:./Samak2.0/$path4/
}
./CopyAllNecessary2Server.sh $1



