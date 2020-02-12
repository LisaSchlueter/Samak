### Simulation and Analysis in Matlab for Katrin (SAMAK; Bulg.: castle) ###

Samak is a framework that is devoted to:
- Monte Carlo simulation for Tritium and Krypton spectra 
- Systematics investigations
- Fit KATRIN Tritium and Krypton data
- Construct confidence intervals and belts

Contact: thierry.lasserre@cea.fr and lisa.schlueter@mpp.mpg.de

Structure:
----------
- .git/             	- git repository program
- doc			- samak documentation
- inputs            	- experiment specific data (FSD, Response functions, etc.) and covariance matrices
- fitting 	        - fit class and fit functions 
- krypton-data          - KATRIN Krypton data
- RunAnalysis 	        - Classes helpful for KATRIN tritium data analysis
- simulation        	- KATRIN specific simulation code
- studies           	- user-specific analysis code
- knm1ana               - KNM1 analysis and studies
- knm2ana               - KNM2 analysis and studies
- tools             	- miscellaneous functions
- tritium data  		- KATRIN tritium data and KATRIN MC data (for example twins)

BINARY Files
----------
- crucial binary files are stored with git lfs:
- if not installed on your machine: download and install git lfs
- run 'git lfs pull' to retreive all binary files

