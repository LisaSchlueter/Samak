% Samak Analysis

% Define a Run Analysis
R=RunAnalysis('RunNr',51870,'DataType','Twin');

% Extract Response Function
[teBins qUbins] = size(R.ModelObj.RF); 
fprintf('Samak RF: teBins=%.0f qUbins=%.0f\n',teBins,qUbins);
R.ModelObj.DisplayWGTSMACEInfo

% Plot Response Function
R.ModelObj.PlotKTF