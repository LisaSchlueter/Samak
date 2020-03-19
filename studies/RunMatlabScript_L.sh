#!/bin/bash

#Simple Script to Run Matlab Script in Batchmode on server
#Execute with: ./RunMatlabScript.sh >> outputfile.txt &
#Matlab output & Error messages are directed to outputfile.txt   

path='/home/iwsatlas1/schluete/samak/studies/tf'
scriptname='fite0_RFCovMat.m'
output='MatlabOutput_fite0_RFCM.txt'

nice nohup /home/iwsatlas1/fuchsdom/matlab2017b/R2017b/bin/matlab -nosplash -nodesktop -noFigureWindows -r "try; cd('$path'); tic; run('$scriptname'); catch ME;  fprintf('MATLAB error: %s \n',ME.message); end; toc; quit" >> $output &
