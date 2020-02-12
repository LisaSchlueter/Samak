#!/bin/bash
# copy hdf5 files from christians homefolder on the tristan server to my local h5 folder
# input argument: 1) login 2) folder which you want to copy
GetH5FilesFromServer(){
		
# synchronize hdf5 data folder with webtrium data folder
rsync -i -a -v $1@pcltr-01.mpp.mpg.de:/home/iwsatlas1/karlch/.local/share/webtrium/$2/RunFiles/H5/ ./hdf5/$2/

# convert to mat file DOESNT WORK LIKE THIS
#matlab -nosplash -nodesktop -noFigureWindows -r "try; run('HDF5readallruns($2)'); catch ME;  fprintf('MATLAB error: %s \n',ME.message); end;  quit" 

}
GetH5FilesFromServer $1 $2

# Cron Job Command /10 min - Use VI as crontab -e editor
#*/10 * * * * /Users/thierrylasserre1/Samak2.0/tritium-data/GetH5FilesWebtrium.sh lasserre Knm1 >/dev/null 2>&1

