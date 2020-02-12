function  ft_SPR_RhoD_LoopOverqUrange(varargin)
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

p = inputParser;
p.addParameter('RunList','StackCD100all');
p.addParameter('mode','Data',@(x)ismember(x,{'Sim','Data'}));
p.addParameter('fixPar','1 5 6',@(x)ischar(x));
p.parse(varargin{:});
RunList       = p.Results.RunList;
mode          = p.Results.mode;
fixPar        = p.Results.fixPar;

ranges = 1:14;
Nranges = length(ranges);
index =  -[1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
myindex = index(ranges);

% Choose fixed parameters (optional)
FixParLabel = fixPar(~isspace(fixPar));


% Choose fitter
fitter = 'minuit';
    
MRA_init = MultiRunAnalysis('AnaFlag','StackPixel',...
    'chi2','chi2Stat','RunList',RunList,...
    'fitter',fitter,'exclDataStart',1,...
    'fixPar',fixPar,...
    'DataEffcorr','RunSummary');

% Check if file is already there
FinalTableName = sprintf('ft_SPR_RhoD_LoopOverqUrange_fixPar%s_%s.xls',FixParLabel,mode);
if~exist(FinalTableName,'file')
    
    % If simulation
    StatFluct  = 'ON';
    SysFluct   = 'ON';
    nSamples   = 1;
    
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
                    'DataEffcorr','RunSummary');
                
                MRA_sys(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                    'chi2','chi2CMShape','RunList',RunList,...
                    'fitter',fitter,'exclDataStart',dataStart,...
                    'fixPar',fixPar,...
                    'DataEffcorr','RunSummary');
                
                MRA_sys(i).ComputeCM('DataDriven','OFF','WGTS_TASR_RelErr',0.005);

            case 'Sim'
                
                MRA_stat(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                    'chi2','chi2Stat','RunList',RunList,...
                    'fitter',fitter,'exclDataStart',dataStart,...
                    'fixPar',fixPar,...
                    'DataEffcorr','OFF');
                
                MRA_sys(i) = MultiRunAnalysis('AnaFlag','StackPixel',...
                    'chi2','chi2CMShape','RunList',RunList,...
                    'fitter',fitter,'exclDataStart',dataStart,...
                    'fixPar',fixPar,...
                    'DataEffcorr','OFF');
                
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
                
                MRA_sys(i).RunData.TBDIS   = TBDIS_Sim_sys;
                MRA_sys(i).RunData.TBDISE  = sqrt(TBDIS_Sim_sys);
        end
        
        MRA_stat(i).RhoDScan;
        MRA_sys(i).RhoDScan;
        
    end
    %% Gather Results
    cd_stat    = []; cd1smin_stat  = []; cd1smax_stat  = [];
    cd_sys     = []; cd1smin_sys   = []; cd1smax_sys   = [];
    
    for i=1:numel(ranges)
        
        cd_stat         = [cd_stat       MRA_stat(i).CDmin];
        cd_sys          = [cd_sys        MRA_sys(i).CDmin];
        
        cd1smin_stat    = [cd1smin_stat  MRA_stat(i).rhoDlowerUnc-MRA_stat(i).CDmin];
        cd1smin_sys     = [cd1smin_sys   MRA_sys(i).rhoDlowerUnc-MRA_sys(i).CDmin];
        
        cd1smax_stat    = [cd1smax_stat  MRA_stat(i).rhoDupperUnc-MRA_stat(i).CDmin];
        cd1smax_sys     = [cd1smax_sys   MRA_sys(i).rhoDupperUnc-MRA_sys(i).CDmin];
        
    end
    
    % Latex Table
    system('touch ft_SPR_RhoD_LoopOverqUrange.xls; rm -rf ft_SPR_RhoD_LoopOverqUrange.xls');
    TexTable = table(myindex',cd_stat',cd1smin_stat',cd1smax_stat',cd_sys',cd1smin_sys',cd1smax_sys',...
        'VariableNames',{'ranges','cdstat','cd1sminstat','cd1smaxstat','cdsys','cd1sminsys','cd1smaxsys'});
    writetable(TexTable,'ft_SPR_RhoD_LoopOverqUrange.xls');
    system(sprintf('mv ft_SPR_RhoD_LoopOverqUrange.xls %s',FinalTableName));
    
else
    
    [ndata, text, alldata] = xlsread(FinalTableName);
    
    cd_stat         = cell2mat(alldata(:,2));
    cd_sys          = cell2mat(alldata(:,5));
    
    cd1smin_stat    = abs(cell2mat(alldata(:,3)));
    cd1smin_sys     = abs(cell2mat(alldata(:,6)));
    
    cd1smax_stat    = abs(cell2mat(alldata(:,4)));
    cd1smax_sys     = abs(cell2mat(alldata(:,7)));    
end

    
%% RhoD Verus Range
myMainTitle=[sprintf('KATRIN First Tritium Samak RhoD Fit - %s - qU Scan - FixPar = %s',mode,fixPar)];                 
maintitle=myMainTitle;
savefile=sprintf('plots/RhoD_%s_qUScan_fixpar%s_%s.pdf',RunList,FixParLabel,mode);
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
hstat=errorbar(myindex+2,cd_stat,cd1smin_stat,+cd1smax_stat,'o','Color',rgb('IndianRed'),'MarkerSize',10,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2)
hold on
hsys=errorbar(myindex-2,cd_sys,cd1smin_sys,cd1smax_sys,'s','Color',rgb('DarkBlue'),'MarkerSize',10,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2)
line([-2000 0],[MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2],'Color','red','LineStyle','--');
[ll5 la5]= boundedline([-2000 0],...
                [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2],...
                [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2]*0.05,...
                'alpha','cmap',rgb('green'));
            ll5.LineStyle='--'; ll5.LineWidth=3; 
hold off
ylabel('Rho x D molecules/cm^2')
xlabel('Lower qU (V)'); 
sp1 = sprintf('stat - <RhoD>_{uncorrelated}=%.3g mol/cm2 - std = %.2g mol/cm2',wmean(cd_stat,1./((cd1smin_stat+cd1smax_stat)/2).^2),std(cd_stat)); 
sp2 = sprintf('stat+sys - <RhoD>_{uncorrelated}=%.3g mol/cm2 - std = %.2g mol/cm2',wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2),std(cd_sys)); 
leg=legend([hstat,hsys,la5],sp1,sp2,'RhoD 5% uncertainty','Location','NorthWest');                   
leg.Color = 'none'; leg.FontSize = 18; legend boxoff;
grid on
axis([-1700 -50 wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2)-5*std(cd_sys) wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2)+5*std(cd_sys)]);
PrettyFigureFormat
set(gca,'FontSize',24);
set(gca, 'XScale', 'lin')
publish_figurePDF(gcf,savefile);

%% RhoD Verus Range
myMainTitle=[sprintf('KATRIN First Tritium Column Density Fit - %s - qU Scan - FixPar = %s (Samak)',mode,fixPar)];                 
maintitle=myMainTitle;
savefile=sprintf('plots/RhoD_%s_qUScan_fixpar%s_%s_2.pdf',RunList,FixParLabel,mode);
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
hsys=errorbar(myindex,cd_sys,cd1smin_sys,cd1smax_sys,'s','Color',rgb('DarkBlue'),'MarkerSize',13,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('PowderBlue'),'LineWidth',2)
hold on
line([-2000 0],[MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2],'Color','red','LineStyle','--');
[ll1 la1]= boundedline([-2000 0],...
    [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2],...
    [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2]*0.01,...
    'alpha','cmap',rgb('Orange'));
ll5.LineStyle='--'; ll5.LineWidth=3;
[ll5 la5]= boundedline([-2000 0],...
    [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2],...
    [MRA_init.RunData.WGTS_CD_MolPerCm2 MRA_init.RunData.WGTS_CD_MolPerCm2]*0.05,...
    'alpha','cmap',rgb('goldenRod'));
ll5.LineStyle='--'; ll5.LineWidth=3;
hold off
ylabel('Column Density (mol/cm^2)')
xlabel('Lower qU (V)'); 
sp2 = sprintf('Column Density Fit (stat+sys), Uncorrelated Average = %.3g mol/cm^2',wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2)); 
leg=legend([hsys,ll5,la1,la5],sp2,sprintf('Column Density Expected Value : %.3g mol/cm^2',MRA_init.RunData.WGTS_CD_MolPerCm2),'Column Density 1% uncertainty band','Column Density 5% uncertainty band','Location','NorthWest');                   
leg.Color = 'none'; leg.FontSize = 18; legend boxoff;
grid on
%axis([-1700 -50 wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2)-6*std(cd_sys) wmean(cd_sys,1./((abs(cd1smin_sys)+abs(cd1smax_sys))/2).^2)+6*std(cd_sys)]);
axis([-1700 -50 4.1e17 4.95e17]);
PrettyFigureFormat
set(gca,'FontSize',26);
set(gca, 'XScale', 'lin')
publish_figurePDF(gcf,savefile);