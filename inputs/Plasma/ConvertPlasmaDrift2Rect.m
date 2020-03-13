function RectWidth = ConvertPlasmaDrift2Rect(Slope,TimeSec)

% Convert slope (eV/day) to width of rectangular function of FSDs (eV)
TimeDay = TimeSec./(60*60*24);
RectWidth = Slope.*TimeDay;

end