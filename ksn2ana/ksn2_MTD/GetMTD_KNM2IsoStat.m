% construct new KNM2 IsoStat MTD


% load knm2 object
fname = [getenv('SamakPath'),'local_Syst_ThesisPlots/results/KNM2Obj.mat'];
d = importdata(fname);
M = d.A.ModelObj; % TBD object

RunData = d.A.RunData;

qU = M.qU; % use original qU
qUfrac_reg = M.qUfrac;
TimeSec = M.TimeSec;

SignalIdx = qU-18574<0;
qUfracNormSignal = sum(qUfrac_reg(SignalIdx)); % <1, because of deleted 300 eV point
Rate= M.TBDIS(SignalIdx)./(qUfrac_reg(SignalIdx).*TimeSec);%-M.BKG_RateSec;
qUfrac     = zeros(size(qU));
qUfrac(SignalIdx)   = 1./Rate;
qUfrac(SignalIdx)  = qUfrac(SignalIdx).*(qUfracNormSignal./sum(qUfrac(SignalIdx))); % normalize to original 
qUfrac(~SignalIdx)  = qUfrac_reg(~SignalIdx);


%sum(qUfrac)-sum(qUfrac_reg)
GetFigure;
plot(qU-18574,qUfrac_reg,'-o')
hold on;
plot(qU-18574,qUfrac,'-.o')
xlim([-40 20])
sname = [getenv('SamakPath'),'ksn2ana/ksn2_MTD/ref_files/KSN2IsoStatMTD_new.mat'];
save(sname,'qU','qUfrac','TimeSec','RunData');

