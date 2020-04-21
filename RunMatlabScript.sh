#!/bin/bash
#input arguments: 1. script name 2. name of output file (e.g. MatlabOutput.txt)
function RunMatlabScript(){
/remote/ceph/group/katrin/software/MATLAB/./matlab_start.sh $1 $2
}
RunMatlabScript $1 $2






