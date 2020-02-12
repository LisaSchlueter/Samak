function IEsuppressionFactor = GetIEsuppressionFactor(IEvoltage)
%------------------------------------------------------------------------------------------
% function to compute the inner electrodes (IE) backgrond suppression
% factor as a function of IE voltage
% Parameterization from: Nikolaus Trost (private communication + phD Thesis for reference)
% ------------------------------------------------------------------------------------------
% Lisa Schlüter (September 2018)
%-------------------------------------------------------------------------------------------

%work in progress....
 IEPar = [0.85429773, -0.13817633,  8.07469579]; %baked spetrometer
 %IEPar = [1.06770558, -0.10807942,  3.27214832]; %unbaked spectrometer
 
 GetRate   = @(x)arrayfun(@(x) (IEPar(1)*(x+IEPar(3))^IEPar(2)),x); % Compute Background Rate, input: IE voltage
 GetFactor = @(x)arrayfun(@(x) GetRate(x)./GetRate(0),x);           % Compute supression factor
 IEsuppressionFactor = GetFactor(IEvoltage);
end
 