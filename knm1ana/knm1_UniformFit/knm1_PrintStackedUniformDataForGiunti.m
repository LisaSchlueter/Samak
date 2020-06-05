format long;
HV              = round(R.RunData.qU(R.exclDataStart:end),3);
CountPerHV      = round(R.RunData.TBDIS(R.exclDataStart:end));
TimePerHVs      = round(R.RunData.TimeSec .* R.RunData.qUfrac(R.exclDataStart:end));
Rate            = round(R.RunData.TBDIS(R.exclDataStart:end)./ (R.RunData.TimeSec .* R.RunData.qUfrac(R.exclDataStart:end)),4);

Tfull = table(HV,CountPerHV,TimePerHVs,Rate,'VariableNames',{'Retarding Energy (eV)';'Counts';'Time (Seconds)';'Rate (s-1)'});
writetable(Tfull,'KATRIN_KNM1full_27HV.dat','Delimiter',' ');

Tshort = table(HV,Rate,'VariableNames',{'Retarding Energy (eV)';'Rate (s-1)'});
writetable(Tshort,'KATRIN_KNM1_27HV.dat','Delimiter',' ');

Tshort
