#!/bin/bash

#Simple Script to Run Matlab Script in Batchmode on server
#Make script executable with chmod 755 scriptname.sh 
#Execute with: ./RunMatlabScript_Lisa.sh >> outputfile.txt &   
#Matlab output & Error messages are directed to outputfile.txt

path='/home/lisa/Dokumente/Studium/Masterstudium/Masterarbeit/samak/studies/tf'
scriptname='How2Use_CovarianceMatrixClass.m'
output='MatlabOutput.txt'

nohup ~/Programme/Matlab/bin/matlab -nosplash -nodesktop -noFigureWindows -r "try; cd('$path'); tic; run('$scriptname'); catch ME;  fprintf('MATLAB error: %s \n Aborting!\n',ME.message); end; toc; quit" >> $output &
