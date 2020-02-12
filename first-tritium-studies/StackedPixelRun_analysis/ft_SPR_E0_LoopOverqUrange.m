function  ft_SPR_E0_LoopOverqUrange(varargin)
% Stacked Pixel Range
% Fit parameters as a function of qU range
% Stat only & Stat+Sys
%
% Input:
%
%
% Output: Excel Table
%
%
addpath(genpath('../../../Samak2.0'));
index =  -[1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
ELossFlag = 'Abdurashitov'; % 'Aseev','Abdurashitov','CW_GLT'

p = inputParser;
p.addParameter('RunList','StackCD100all'); %StackCD100all %StackCD100_3hours
p.addParameter('mode','Data',@(x)ismember(x,{'Sim','Data'}));
p.addParameter('sysFlag','OFF',@(x)ismember(x,{'ON','OFF'})); %
p.addParameter('fixPar','1 5 6',@(x)ischar(x));
p.parse(varargin{:});
RunList       = p.Results.RunList;
mode          = p.Results.mode;
sysFlag       = p.Results.sysFlag;
fixPar        = p.Results.fixPar;

ranges = 1:18; 
Nranges = length(ranges);

% Choose fixed parameters (optional)
FixParLabel = fixPar(~isspace(fixPar));

% If simulation
StatFluct  = 'ON';
SysFluct   = 'ON';
nSamples   = 1;

% Choose fitter
fitter = 'minuit';

i=0;
for range = 1:Nranges
    i=i+1;
    % Choose data range to analyze
    dataStart = ranges(range);
    
    
    switch mode
        case 'Data'
            
            MRA_stat(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                'chi2','chi2Stat','RunList',RunList,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,...
                'DataEffcorr','RunSummary',...
                'ELossFlag',ELossFlag);

            switch sysFlag
                case 'ON'
                    MRA_sys(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                        'chi2','chi2CMShape','RunList',RunList,...
                        'fitter',fitter,'exclDataStart',dataStart,...
                        'fixPar',fixPar,...
                        'DataEffcorr','RunSummary',...
                        'ELossFlag',ELossFlag);
            end
            
        case 'Sim'
            
            MRA_stat(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                'chi2','chi2Stat','RunList',RunList,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,...
                'DataEffcorr','OFF',...
                'ELossFlag',ELossFlag);
            
            switch sysFlag
                case 'ON'
                    MRA_sys(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                        'chi2','chi2CMShape','RunList',RunList,...
                        'fitter',fitter,'exclDataStart',dataStart,...
                        'fixPar',fixPar,...
                        'DataEffcorr','OFF',...
                        'ELossFlag',ELossFlag);
            end
            
            %Compute MC Data set
            switch StatFluct
                case 'ON'
                    if strcmp(SysFluct,'OFF')
                        TBDIS_Sim_stat = mvnrnd(MRA_stat(i).ModelObj.TBDIS,MRA_stat(i).ModelObj.TBDIS',nSamples)'; % nSamples simulated integral spectra
                        TBDIS_Sim_sys = mvnrnd(MRA_sys(i).ModelObj.TBDIS,MRA_sys(i).ModelObj.TBDIS',nSamples)'; % nSamples simulated integral spectra
                    elseif strcmp(SysFluct,'ON')
                        TBDIS_Sim_stat = mvnrnd(MRA_stat(i).ModelObj.TBDIS,MRA_stat(i).ModelObj.TBDIS',nSamples)'; % nSamples simulated integral spectra
                        TBDIS_Sim_sys = mvnrnd(MRA_sys(i).ModelObj.TBDIS,MRA_sys(i).FitCM+diag(MRA_sys(i).ModelObj.TBDIS),nSamples)';
                    end
                case 'OFF'
                    if strcmp(SysFluct,'OFF')
                        TBDIS_Sim_stat = repmat(MRA_stat(i).ModelObj.TBDIS,1,nSamples);
                        TBDIS_Sim_sys = repmat(MRA_sys(i).ModelObj.TBDIS,1,nSamples);
                    elseif  strcmp(SysFluct,'ON')
                        TBDIS_Sim_stat = mvnrnd(MRA_stat(i).ModelObj.TBDIS,MRA_stat(i).FitCM,nSamples)';
                        TBDIS_Sim_sys = mvnrnd(MRA_sys(i).ModelObj.TBDIS,MRA_sys(i).FitCM,nSamples)';
                    end
            end
            
            MRA_stat(i).RunData.TBDIS   = TBDIS_Sim_stat;
            MRA_stat(i).RunData.TBDISE  = sqrt(TBDIS_Sim_stat);
            
            switch sysFlag
                case 'ON'
                    MRA_sys(i).RunData.TBDIS   = TBDIS_Sim_sys;
                    MRA_sys(i).RunData.TBDISE  = sqrt(TBDIS_Sim_sys);
            end
            
    end

    MRA_stat(i).Fit();
    
    switch sysFlag
        case 'ON'
            MRA_sys(i).Fit();
    end
    
end

switch sysFlag
        case 'OFF'
            MRA_sys = MRA_stat;
end

%% Gather Results
range = [];
e0_stat    = []; e0err_stat  = [];
m_stat     = []; merr_stat   = [];
b_stat     = []; berr_stat   = [];
n_stat     = []; nerr_stat   = [];
pgs_stat   = []; pgserr_stat = [];
chi2_stat  = []; dof_stat    = [];
e0_sys     = []; e0err_sys   = [];
m_sys      = []; merr_sys    = [];
b_sys      = []; berr_sys    = [];
n_sys      = []; nerr_sys    = [];
pgs_sys    = []; pgserr_sys  = [];
chi2_sys   = []; dof_sys     = [];

for i=1:numel(ranges)
    
    % m 
    m_stat         = [m_stat       MRA_stat(i).FitResult.par(1)];   
    merr_stat      = [merr_stat    MRA_stat(i).FitResult.err(1)];   
    m_sys          = [m_sys        MRA_sys(i).FitResult.par(1)];   
    merr_sys       = [merr_sys     MRA_sys(i).FitResult.err(1)];    
    % e0 
    e0_stat        = [e0_stat      MRA_stat(i).FitResult.par(2)+MRA_stat(i).ModelObj.Q_i];   
    e0err_stat     = [e0err_stat   MRA_stat(i).FitResult.err(2)];   
    e0_sys         = [e0_sys       MRA_sys(i).FitResult.par(2)+MRA_sys(i).ModelObj.Q_i];   
    e0err_sys      = [e0err_sys    MRA_sys(i).FitResult.err(2)];  
    % b 
    b_stat         = [b_stat       MRA_stat(i).FitResult.par(3)+MRA_stat(i).ModelObj.BKG_RateSec_i];   
    berr_stat      = [berr_stat    MRA_stat(i).FitResult.err(3)];   
    b_sys          = [b_sys        MRA_sys(i).FitResult.par(3)+MRA_sys(i).ModelObj.BKG_RateSec_i];   
    berr_sys       = [berr_sys     MRA_sys(i).FitResult.err(3)]; 
    % n
    n_stat         = [n_stat       MRA_stat(i).FitResult.par(4)];   
    nerr_stat      = [nerr_stat    MRA_stat(i).FitResult.err(4)];   
    n_sys          = [n_sys        MRA_sys(i).FitResult.par(4)];   
    nerr_sys       = [nerr_sys     MRA_sys(i).FitResult.err(4)]; 
    % gs
    pgs_stat       = [pgs_stat     MRA_stat(i).FitResult.par(5)+MRA_stat(i).ModelObj.DTNormGS_i];   
    pgserr_stat    = [pgserr_stat  MRA_stat(i).FitResult.err(5)];   
    pgs_sys        = [pgs_sys      MRA_sys(i).FitResult.par(5)+MRA_sys(i).ModelObj.DTNormGS_i];   
    pgserr_sys     = [pgserr_sys   MRA_sys(i).FitResult.err(5)];
    %chi2
    chi2_stat       = [chi2_stat   MRA_stat(i).FitResult.chi2min];
    chi2_sys        = [chi2_sys    MRA_sys(i).FitResult.chi2min];
    %dof
    dof_stat        = [dof_stat   MRA_stat(i).FitResult.dof];
    dof_sys         = [dof_sys    MRA_sys(i).FitResult.dof];
end

% Latex Table
myindex = index(ranges);
    system('touch ft_SPR_E0_LoopOverqUrange.xls; rm -rf ft_SPR_E0_LoopOverqUrange.xls');
    TexTable = table(myindex',m_stat',merr_stat',m_sys',merr_sys',e0_stat',e0err_stat',e0_sys',e0err_sys',b_stat',berr_stat',b_sys',berr_sys',n_stat',nerr_stat',n_sys',nerr_sys',pgs_stat',pgserr_stat',pgs_sys',pgserr_sys',chi2_stat',chi2_sys',dof_stat',dof_sys',...
        'VariableNames',{'ranges','m_stat','merr_stat','m_sys','merr_sys','e0_stat','e0err_stat','e0_sys','e0err_sys','b_stat','berr_stat','b_sys','berr_sys','n_stat','nerr_stat','n_sys','nerr_sys','pgs_stat','pgserr_stat','pgs_sys','pgserr_sys','chi2_stat','chi2_sys','dof_stat','dof_sys'});
    writetable(TexTable,'ft_SPR_E0_LoopOverqUrange.xls');
    FinalTableName = sprintf('ft_SPR_E0_LoopOverqUrange_fixPar%s_%s.xls',FixParLabel,mode);
    system(sprintf('mv ft_SPR_E0_LoopOverqUrange.xls %s',FinalTableName));
    
%% Endpoint Verus Range
myMainTitle=[sprintf('KATRIN First Tritium Samak Fit - %s - qU Scan - FixPar = %s',mode,fixPar)];                 
maintitle=myMainTitle;
savefile=sprintf('plots/E0_%s_qUScan_fixpar%s_%s.pdf',RunList,FixParLabel,mode);
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
subplot(2,1,1)
hstat=errorbar(myindex+2,e0_stat,e0err_stat,'o','Color',rgb('IndianRed'),'MarkerSize',7,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2)
hold on
hsys=errorbar(myindex-2,e0_sys,e0err_sys,'s','Color',rgb('DarkBlue'),'MarkerSize',7,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2)
hold off
ylabel('Effective E_0 (eV)')
xlabel('Lower qU (V)'); 
sp1 = sprintf('stat - <(E_0)_{eff}>_{uncorrelated}=%.2f eV - std = %.2f eV',wmean(e0_stat,1./e0err_stat.^2),std(e0_stat)); 
sp2 = sprintf('stat+sys - <(E_0)_{eff}>_{uncorrelated}=%.2f eV - std = %.2f eV',wmean(e0_sys,1./e0err_sys.^2),std(e0_sys)); 
leg=legend([hstat,hsys],sp1,sp2,'Location','NorthWest');                   
leg.Color = 'none'; leg.FontSize = 18; legend boxoff;
grid on
axis([-1700 -50 (wmean(e0_sys,1./e0err_sys.^2)-5*std(e0_sys)) (wmean(e0_sys,1./e0err_sys.^2)+5*std(e0_sys))]);
PrettyFigureFormat
set(gca,'FontSize',24);
set(gca, 'XScale', 'log')
subplot(2,1,2)
hstat=plot(myindex+3,chi2pvalue(chi2_stat,dof_stat),'o','Color',rgb('IndianRed'),'MarkerSize',7,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2)
hold on
hsys=plot(myindex-3,chi2pvalue(chi2_sys,dof_sys),'s','Color',rgb('DarkBlue'),'MarkerSize',7,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2)
hold off
ylabel('p-value')
xlabel('Lower qU (V)'); 
sp1 = sprintf('stat'); 
sp2 = sprintf('stat+sys'); 
leg=legend([hstat,hsys],sp1,sp2,'Location','SouthEast');                   
leg.Color = 'none'; legend boxoff;
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
grid on
xlim([-1700 -50]);
PrettyFigureFormat
set(gca,'FontSize',24);
publish_figurePDF(gcf,savefile);

%% Pgs Verus Range
myMainTitle=[sprintf('KATRIN First Tritium Samak Fit - %s - qU Scan - FixPar = %s',mode,fixPar)];                 
maintitle=myMainTitle;
savefile=sprintf('plots/Pgs_%s_qUScan_fixpar%s_%s.pdf',RunList,FixParLabel,mode);
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
hstat=errorbar(myindex+2,pgs_stat,pgserr_stat,'o','Color',rgb('IndianRed'),'MarkerSize',7,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2)
% hold on
% hsys=errorbar(myindex-2,pgssys,pgssyserr,'s','Color',rgb('DarkBlue'),'MarkerSize',7,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2)
% hold off
ylabel('P_{gs} ')
xlabel('Lower qU (V)'); 
sp1 = sprintf('stat - <P_{gs}>_{uncorrelated}=%.2f eV - std = %.2f eV',wmean(pgs_stat,1./pgserr_stat.^2),std(pgs_stat)); 
sp2 = sprintf('stat+sys - <(E_0)_{eff}>_{uncorrelated}=%.2f eV - std = %.2f eV',wmean(pgs_sys,1./pgserr_sys.^2),std(pgs_sys)); 
%leg=legend([hstat,hsys],sp1,sp2,'Location','NorthWest');                   
%leg.Color = 'none'; leg.FontSize = 18; legend boxoff;
grid on
axis([-1700 -50 0.45 0.7]);
PrettyFigureFormat
set(gca,'FontSize',24);
set(gca, 'XScale', 'log')
publish_figurePDF(gcf,savefile);
                    