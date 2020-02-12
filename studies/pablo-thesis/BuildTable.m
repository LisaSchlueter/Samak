function savethis = BuildTable(varargin)

if nargin == 0
doExcel = true;
doSaveWorkspace = false;
else
    doSaveWorkspace = true;
    doExcel = false;
end

addpath(genpath('../../../Samak2.0'));


% if ~exist('SPAll','var')
    data_version = '2';
    SPAll = load(['data/StackedPixelsResultsThesis',data_version]);
    SPSR = load(['data/StackedPixelsStackedRunsResultsThesis',data_version]);
    R = load(['data/StackedRingsStackedRunsResultsThesis',data_version]);
    SPix = load(['data/SinglePixelResultsThesis',data_version]);
    MPix = load(['data/MultiPixelResultsThesis',data_version]);
% end

E0_i = 18573.7;

%% Statistics
% Mean StackP All runs

w_SPAll_s = 1./(SPAll.TableShortStat(:,6)).^2;
E_SPAll_s = wmean(SPAll.TableShortStat(:,2),w_SPAll_s)+ E0_i;
s_SPAll_s = 1./sqrt(sum(w_SPAll_s));
chi2_SPAll_s = 'See fig. \ref{f:Chi2E0Stat}';

w_SPAll_m = 1./(SPAll.TableMedStat(:,6)).^2;
E_SPAll_m = wmean(SPAll.TableMedStat(:,2),w_SPAll_m)+ E0_i;
s_SPAll_m = 1./sqrt(sum(w_SPAll_m));
chi2_SPAll_m = 'See fig. \ref{f:Chi2E0Stat}';


w_SPAll_l = 1./(SPAll.TableLongStat(:,6)).^2;
E_SPAll_l = wmean(SPAll.TableLongStat(:,2),w_SPAll_l)+ E0_i;
s_SPAll_l = 1./sqrt(sum(w_SPAll_l));
chi2_SPAll_l = 'See fig. \ref{f:Chi2E0Stat}';

% StackP StackRuns

E_SPSR_s = SPSR.ResultShortStat(1,2)+ E0_i;
s_SPSR_s = SPSR.ResultShortStat(1,6);
chi2_SPSR_s = [num2str(round(SPSR.ResultShortStat(1,9),2)),'/',num2str(SPSR.ResultShortStat(1,10))];

E_SPSR_m = SPSR.ResultMedStat(1,2)+ E0_i;
s_SPSR_m = SPSR.ResultMedStat(1,6);
chi2_SPSR_m = [num2str(round(SPSR.ResultMedStat(1,9),2)),'/',num2str(SPSR.ResultMedStat(1,10))];

E_SPSR_l = SPSR.ResultLongStat(1,2)+ E0_i;
s_SPSR_l = SPSR.ResultLongStat(1,6);
chi2_SPSR_l = [num2str(round(SPSR.ResultLongStat(1,9),2)),'/',num2str(SPSR.ResultLongStat(1,10))];

% Stacked Ring

w_R_s = 1./(R.TableShortStat(:,8)).^2;
E_R_s = wmean(R.TableShortStat(:,2),w_R_s)+ E0_i;
s_R_s = 1./sqrt(sum(w_R_s));
chi2_R_s = 'See fig. \ref{f:RingE0Stat}';

w_R_m = 1./(R.TableMedStat(:,8)).^2;
E_R_m = wmean(R.TableMedStat(:,2),w_R_m)+ E0_i;
s_R_m = 1./sqrt(sum(w_R_m));
chi2_R_m = 'See fig. \ref{f:RingE0Stat}';

% Single Pixel Stacked Runs

w_SPix_s = 1./(SPix.TableShortStat(:,8)).^2;
E_SPix_s = wmean(SPix.TableShortStat(:,2),w_SPix_s)+ E0_i;
s_SPix_s = 1./sqrt(sum(w_SPix_s));
chi2_SPix_s = 'See fig. \ref{f:PixelE0Stat}';

w_SPix_m = 1./(SPix.TableMedStat(:,8)).^2;
E_SPix_m = wmean(SPix.TableMedStat(:,2),w_SPix_m)+ E0_i;
s_SPix_m = 1./sqrt(sum(w_SPix_m));
chi2_SPix_m = 'See fig. \ref{f:PixelE0Stat}';

% Multipixel Stacked Runs

E_MPix_s = MPix.ResultShortStat(:,2) + E0_i;
s_MPix_s = MPix.ResultShortStat(:,246);
chi2_MPix_s = [num2str(round(MPix.ResultShortStat(1,489),2)),'/',num2str(MPix.ResultShortStat(1,490))];

E_MPix_m = MPix.ResultMedStat(:,2) + E0_i;
s_MPix_m = MPix.ResultMedStat(:,246);
chi2_MPix_m = [num2str(round(MPix.ResultMedStat(1,489),2)),'/',num2str(MPix.ResultMedStat(1,490))];

%% Systematics
savethis = zeros(1,12);
% Mean StackP All runs
chi2_sys_SPAll_s = 'See fig. \ref{f:Chi2E0Sys}';
[E_sys_SPAll_s,s_sys_SPAll_s,chi2,dof] = SysUncertainty(SPAll.TableShortSys(:,6),SPAll.TableShortStat(:,6),SPAll.TableShortSys(:,2));
E_sys_SPAll_s = E_sys_SPAll_s + E0_i;
if chi2pvalue(chi2,dof) > 0.05
    savethis(1) = 1;
end

chi2_sys_SPAll_m = 'See fig. \ref{f:Chi2E0Sys}';
[E_sys_SPAll_m,s_sys_SPAll_m,chi2,dof] = SysUncertainty(SPAll.TableMedSys(:,6),SPAll.TableMedStat(:,6),SPAll.TableMedSys(:,2));
E_sys_SPAll_m = E_sys_SPAll_m + E0_i;
if chi2pvalue(chi2,dof) > 0.05
    savethis(2) = 1;
end

chi2_sys_SPAll_l = 'See fig. \ref{f:Chi2E0Sys}';
[E_sys_SPAll_l,s_sys_SPAll_l,chi2,dof] = SysUncertainty(SPAll.TableLongSys(:,6),SPAll.TableLongStat(:,6),SPAll.TableLongSys(:,2));
E_sys_SPAll_l = E_sys_SPAll_l + E0_i;
if chi2pvalue(chi2,dof) > 0.05
    savethis(3) = 1;
end

% StackP StackRuns

E_sys_SPSR_s = SPSR.ResultShortSys(1,2) + E0_i;
s_sys_SPSR_s = SPSR.ResultShortSys(1,6);
chi2_sys_SPSR_s = [num2str(round(SPSR.ResultShortSys(1,9),2)),'/',num2str(SPSR.ResultShortSys(1,10))];
if chi2pvalue(SPSR.ResultShortSys(1,9),SPSR.ResultShortSys(1,10)) > 0.05
    savethis(4) = 1;
end

E_sys_SPSR_m = SPSR.ResultMedSys(1,2) + E0_i;
s_sys_SPSR_m = SPSR.ResultMedSys(1,6);
chi2_sys_SPSR_m = [num2str(round(SPSR.ResultMedSys(1,9),2)),'/',num2str(SPSR.ResultMedSys(1,10))];
if chi2pvalue(SPSR.ResultMedSys(1,9),SPSR.ResultMedSys(1,10)) > 0.05
    savethis(5) = 1;
end

E_sys_SPSR_l = SPSR.ResultLongSys(1,2) + E0_i;
s_sys_SPSR_l = SPSR.ResultLongSys(1,6);
chi2_sys_SPSR_l = [num2str(round(SPSR.ResultLongSys(1,9),2)),'/',num2str(SPSR.ResultLongSys(1,10))];
if chi2pvalue(SPSR.ResultLongSys(1,9),SPSR.ResultLongSys(1,10)) > 0.05
    savethis(6) = 1;
end

% Stacked Ring

[E_sys_R_s,s_sys_R_s,chi2,dof] = SysUncertainty(R.TableShortSys(:,8),R.TableShortStat(:,8),R.TableShortSys(:,2));
E_sys_R_s = E_sys_R_s + E0_i;
chi2_sys_R_s = 'See fig. \ref{f:RingE0Sys}';
if chi2pvalue(chi2,dof) > 0.05
    savethis(7) = 1;
end

[E_sys_R_m,s_sys_R_m,chi2,dof] = SysUncertainty(R.TableMedSys(:,8),R.TableMedStat(:,8),R.TableMedSys(:,2));
E_sys_R_m = E_sys_R_m + E0_i;
chi2_sys_R_m = 'See fig. \ref{f:RingE0Sys}';
if chi2pvalue(chi2,dof) > 0.05
    savethis(8) = 1;
end

% Single Pixel Stacked Runs

[E_sys_SPix_s,s_sys_SPix_s,chi2,dof] = SysUncertainty(SPix.TableShortSys(:,8),SPix.TableShortStat(:,8),SPix.TableShortSys(:,2));
E_sys_SPix_s = E_sys_SPix_s + E0_i;
chi2_sys_SPix_s = 'See fig. \ref{f:PixelE0Sys}';
if chi2pvalue(chi2,dof) > 0.05
    savethis(9) = 1;
end

[E_sys_SPix_m,s_sys_SPix_m,chi2,dof] = SysUncertainty(SPix.TableMedSys(:,8),SPix.TableMedStat(:,8),SPix.TableMedSys(:,2));
E_sys_SPix_m = E_sys_SPix_m + E0_i;
chi2_sys_SPix_m = 'See fig. \ref{f:PixelE0Sys}';
if chi2pvalue(chi2,dof) > 0.05
    savethis(10) = 1;
end

% Multipixel Stacked Runs

E_sys_MPix_s = MPix.ResultShortSys(:,2) + E0_i;
s_sys_MPix_s = MPix.ResultShortSys(:,246);
chi2_sys_MPix_s = [num2str(round(MPix.ResultShortSys(1,489),2)),'/',num2str(MPix.ResultShortSys(1,490))];
if chi2pvalue(MPix.ResultShortSys(1,489),MPix.ResultShortSys(1,490)) > 0.05
    savethis(11) = 1;
end

E_sys_MPix_m = MPix.ResultMedSys(:,2) + E0_i;
s_sys_MPix_m = MPix.ResultMedSys(:,246);
chi2_sys_MPix_m = [num2str(round(MPix.ResultMedSys(1,489),2)),'/',num2str(MPix.ResultMedSys(1,490))];
if chi2pvalue(MPix.ResultMedSys(1,489),MPix.ResultMedSys(1,490)) > 0.05
    savethis(12) = 1;
end


%% Build Table for excel

HeaderRow = {'Segm.','(eV)','Stat.','$\chi^2$','Stat.+Sys.','$\chi^2$'};
FirstColumn = {'NO','','','','','','','','Ring','','','Single Pix.','','','Multi Pix.','',''}';
AvRuns = {'Runs Av.','','','',''};
Stacked = {'Stacked Runs','','','',''};

round_SPAll = 2; round_SPAll_s = 2;
SPAll_t = {'Sh.',[num2str(round(E_SPAll_s,round_SPAll-1)),'  $\pm$  ',num2str(round(s_SPAll_s,round_SPAll_s))],chi2_SPAll_s,...
    [num2str(round(E_sys_SPAll_s,round_SPAll-1)),'  $\pm$  ',num2str(round(s_sys_SPAll_s,round_SPAll_s))],chi2_sys_SPAll_s;...
    'Med.',[num2str(round(E_SPAll_m,round_SPAll)),'  $\pm$  ',num2str(round(s_SPAll_m,round_SPAll_s))],chi2_SPAll_m,...
    [num2str(round(E_sys_SPAll_m,round_SPAll-1)),'  $\pm$  ',num2str(round(s_sys_SPAll_m,round_SPAll_s))],chi2_sys_SPAll_m;...
    'L.',[num2str(round(E_SPAll_l,round_SPAll)),'  $\pm$  ',num2str(round(s_SPAll_l,round_SPAll_s))],chi2_SPAll_l,...
    [num2str(round(E_sys_SPAll_l,round_SPAll-1)),'  $\pm$  ',num2str(round(s_sys_SPAll_l,round_SPAll_s))],chi2_sys_SPAll_l};

round_SPSR = 2; round_SPSR_s = 2;
SPSR_t = {'Sh.',[num2str(round(E_SPSR_s,round_SPSR-1)),'  $\pm$  ',num2str(round(s_SPSR_s,round_SPSR_s))],chi2_SPSR_s,...
    [num2str(round(E_sys_SPSR_s,round_SPSR-1)),'  $\pm$  ',num2str(round(s_sys_SPSR_s,round_SPSR_s))],chi2_sys_SPSR_s;...
    'Med.',[num2str(round(E_SPSR_m,round_SPSR)),'  $\pm$  ',num2str(round(s_SPSR_m,round_SPSR_s))],chi2_SPSR_m,...
    [num2str(round(E_sys_SPSR_m,round_SPSR-1)),'  $\pm$  ',num2str(round(s_sys_SPSR_m,round_SPSR_s))],chi2_sys_SPSR_m;...
    'L.',[num2str(round(E_SPSR_l,round_SPSR)),'  $\pm$  ',num2str(round(s_SPSR_l,round_SPSR_s))],chi2_SPSR_l,...
    [num2str(round(E_sys_SPSR_l,round_SPSR-1)),'  $\pm$  ',num2str(round(s_sys_SPSR_l,round_SPSR_s))],chi2_sys_SPSR_l};

round_R = 2; round_R_s = 2;
R_t = {'Sh.',[num2str(round(E_R_s,round_R-1)),'  $\pm$  ',num2str(round(s_R_s,round_R_s))],chi2_R_s,...
    [num2str(round(E_sys_R_s,round_R-1)),'  $\pm$  ',num2str(round(s_sys_R_s,round_R_s))],chi2_sys_R_s;...
    'Med.',[num2str(round(E_R_m,round_R)),'  $\pm$  ',num2str(round(s_R_m,round_R_s))],chi2_R_m,...
    [num2str(round(E_sys_R_m,round_R-1)),'  $\pm$  ',num2str(round(s_sys_R_m,round_R_s))],chi2_sys_R_m};

round_SPix = 2; round_SPix_s = 2;
SPix_t = {'Sh.',[num2str(round(E_SPix_s,round_SPix-1)),'  $\pm$  ',num2str(round(s_SPix_s,round_SPix_s))],chi2_SPix_s,...
    [num2str(round(E_sys_SPix_s,round_SPix-1)),'  $\pm$  ',num2str(round(s_sys_SPix_s,round_SPix_s))],chi2_sys_SPix_s;...
    'Med.',[num2str(round(E_SPix_m,round_SPix)),'  $\pm$  ',num2str(round(s_SPix_m,round_SPix_s))],chi2_SPix_m,...
    [num2str(round(E_sys_SPix_m,round_SPix-1)),'  $\pm$  ',num2str(round(s_sys_SPix_m,round_SPix_s))],chi2_sys_SPix_m};

round_MPix = 2; round_MPix_s = 2;
MPix_t = {'Sh.',[num2str(round(E_MPix_s,round_MPix-1)),'  $\pm$  ',num2str(round(s_MPix_s,round_MPix_s))],chi2_MPix_s,...
    [num2str(round(E_sys_MPix_s,round_MPix-1)),'  $\pm$  ',num2str(round(s_sys_MPix_s,round_MPix_s))],chi2_sys_MPix_s;...
    'Med.',[num2str(round(E_MPix_m,round_MPix)),'  $\pm$  ',num2str(round(s_MPix_m,round_MPix_s))],chi2_MPix_m,...
    [num2str(round(E_sys_MPix_m,round_MPix-1)),'  $\pm$  ',num2str(round(s_sys_MPix_m,round_MPix_s))],chi2_sys_MPix_m};


ContentE0 = [AvRuns;SPAll_t;Stacked;SPSR_t;Stacked;R_t;Stacked;SPix_t;Stacked;MPix_t];
ContentE0FC = [FirstColumn,ContentE0];

w_totalStat = (1./([s_SPAll_s,s_SPAll_m,s_SPAll_l,...
    s_SPSR_s,s_SPSR_m,s_SPSR_l,s_R_s,s_R_m,s_SPix_s,s_SPix_m,...
    s_MPix_s,s_MPix_m]).^2)';
E_totalStat = [E_SPAll_s,E_SPAll_m,E_SPAll_l,E_SPSR_s,E_SPSR_m,E_SPSR_l,...
    E_R_s,E_R_m,E_SPix_s,E_SPix_m,E_MPix_s,E_MPix_m]';
totalStat = wmean(E_totalStat,w_totalStat);
s_totalStat = 1./sqrt(sum(w_totalStat));

TotUncVec = [s_sys_SPAll_s,s_sys_SPAll_m,s_sys_SPAll_l,...
    s_sys_SPSR_s,s_sys_SPSR_m,s_sys_SPSR_l,s_sys_R_s,s_sys_R_m,s_sys_SPix_s,s_sys_SPix_m,...
    s_sys_MPix_s,s_sys_MPix_m]';
E_totalSys = [E_sys_SPAll_s,E_sys_SPAll_m,E_sys_SPAll_l,E_sys_SPSR_s,E_sys_SPSR_m,E_sys_SPSR_l,...
    E_sys_R_s,E_sys_R_m,E_sys_SPix_s,E_sys_SPix_m,E_sys_MPix_s,E_sys_MPix_m]';
StaUncVec = [s_SPAll_s,s_SPAll_m,s_SPAll_l,...
    s_SPSR_s,s_SPSR_m,s_SPSR_l,s_R_s,s_R_m,s_SPix_s,s_SPix_m,...
    s_MPix_s,s_MPix_m]';
[E_mean,Unc_mean,chi2,dof] = SysUncertainty(TotUncVec,StaUncVec,E_totalSys);


TableE0 = [HeaderRow;ContentE0FC];
if doExcel
    xlswrite('TableE0.xlsx',TableE0,1)
end

if doSaveWorkspace
    save(['BuildTableWorkspace',data_version],'-v7.3');
end



end



