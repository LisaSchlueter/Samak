#!/bin/bash
# input: 1) MPP login 2) Samak directory
SyncSamakDir(){
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:Samak3.0/$2/ ./$2/
}
SyncSamakDir $1 $2
