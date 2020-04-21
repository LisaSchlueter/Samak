#!/bin/bash
#input arguments: script name
function RunMatlabScriptServer(){
output='MatlabOutput.txt'
nice nohup /home/iwsatlas1/fuchsdom/matlab2017b/R2017b/bin/matlab -nosplash -nodesktop -noFigureWindows -r "try; tic; run('$1'); catch ME;  fprintf('MATLAB error: \
%s \n',ME.message); end; toc; quit" >> $output &
}
RunMatlabScriptServer $1
