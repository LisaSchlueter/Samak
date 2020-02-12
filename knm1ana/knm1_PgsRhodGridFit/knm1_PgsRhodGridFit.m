% Perform Final KNM1 Neutrino Mass Fit
% Scan around best rhoD: -2% -1% nominal +1% +2%
% Scan around best Pgs TT: -1% -0.5% nominal +0.5% +1%
% Plot Results on Grid for
% - Neutrino mass squared best fit
% - Neutrino mass squared error 
% - chi2
% 
% T. Lasserre
% Last Modified, 10/07/2019

%% Initialize Model/Data
RunList = 'KNM1';
clear TwinFSDrhoD;
chi2 = 'chi2CMShape';
if ~exist('TwinFSDrhoD','var')
   TwinFSDrhoD = MultiRunAnalysis(...
       'RunList',RunList,...
       'chi2',chi2,...
       'DataType','Real',...
       'exclDataStart',14,...
       'NonPoissonScaleFactor',1.064);
end
BkgCM                 = 'OFF';
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');

%% Loop on rhoD & Pgs Modeling and Perform Fits
 rhod_local     = TwinFSDrhoD.ModelObj.WGTS_CD_MolPerCm2.*[0.98 0.99 1 1.01 1.02];
 pgs_local      = TwinFSDrhoD.ModelObj.TTNormGS_i.*[0.98 0.99 1 1.01 1.02];
for i = 1:numel(rhod_local)
    for j = 1:numel(pgs_local)
        TwinFSDrhoD = MultiRunAnalysis(...
            'RunList',RunList,...
            'chi2',chi2,...
            'DataType','Real',...
            'exclDataStart',14,...
            'NonPoissonScaleFactor',1.064,...
            'i_TTGS',pgs_local(j)-TwinFSDrhoD.ModelObj.TTNormGS_i,...
            'i_TTES',-(pgs_local(j)-TwinFSDrhoD.ModelObj.TTNormGS_i));
        TwinFSDrhoD.ModelObj.WGTS_CD_MolPerCm2 = rhod_local(i);
        TwinFSDrhoD.ModelObj.AdjustRF; 
        TwinFSDrhoD.ModelObj.ComputeTBDDS();
        TwinFSDrhoD.ModelObj.ComputeTBDIS;
        if ~strcmp(TwinFSDrhoD.chi2,'chi2Stat') %Initialize CovMat with statistics
            TwinFSDrhoD.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM,'FSDNorm_RelErr',0);
        end
        TwinFSDrhoD.fixPar = '5 6 7 8 9 10 11'; % fix only qUOffset!
        % Fit
        TwinFSDrhoD.Fit();
        fprintf(2,'Fit with RhoD = %.3g and Pgs = %.3f (%.3f)...\n',rhod_local(i),TwinFSDrhoD.ModelObj.TTNormGS,pgs_local(j)-TwinFSDrhoD.ModelObj.TTNormGS_i);
        rhod_pgs_struc(i,j) = TwinFSDrhoD.FitResult;
    end
end

%% Display chi2 matrix
for i = 1:numel(rhod_local)
    for j = 1:numel(pgs_local)
        chi2_m(i,j) = rhod_pgs_struc(i,j).chi2min;
    end
end
range    = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart)-TwinFSDrhoD.ModelObj.Q_i);
qUmin = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart));
qUmax = round(TwinFSDrhoD.ModelObj.qU(end));
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Fit Grid Scan - %d Runs - [%.0f - %0.f] eV - Stat+Sys - \chi^2',...
    numel(TwinFSDrhoD.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_UniformStackedGridScan_%d_%s_%.0feVbelowE0-chi2.png',...
    numel(TwinFSDrhoD.RunList),TwinFSDrhoD.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit Grid Scan','NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

cm=0:1/1000:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]); imagesc(chi2_m); c=colorbar; c.Label.String = '\chi ^2';
set(gca,'XTick',[1:numel(rhod_local)],'XTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},...
    'YTick',[1:numel(pgs_local)],'YTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},'YDir','normal');
[x,y] = meshgrid(1:numel(rhod_local),1:numel(pgs_local));
text(x(:),y(:),num2str(chi2_m(:)),'HorizontalAlignment','center');
ylabel('Column Density'); 
xlabel('FSD T-T Ground State Probability'); 
PrettyFigureFormat
export_fig(gcf,savefile,'-m3');

%% Display m2 matrix
for i = 1:numel(rhod_local)
    for j = 1:numel(pgs_local)
        m2_m(i,j) = rhod_pgs_struc(i,j).par(1);
    end
end
range    = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart)-TwinFSDrhoD.ModelObj.Q_i);
qUmin = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart));
qUmax = round(TwinFSDrhoD.ModelObj.qU(end));
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Fit Grid Scan - %d Runs - [%.0f - %0.f] eV - Stat+Sys - m^2 (eV^2)',...
    numel(TwinFSDrhoD.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_UniformStackedGridScan_%d_%s_%.0feVbelowE0-m2.png',...
    numel(TwinFSDrhoD.RunList),TwinFSDrhoD.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit Grid Scan','NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

cm=0:1/1000:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
imagesc(m2_m); c=colorbar;  c.Label.String = 'eV^2';
set(gca,'XTick',[1:numel(rhod_local)],'XTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},...
    'YTick',[1:numel(pgs_local)],'YTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},'YDir','normal');
[x,y] = meshgrid(1:numel(rhod_local),1:numel(pgs_local));
text(x(:),y(:),num2str(m2_m(:)),'HorizontalAlignment','center');
ylabel('Column Density'); 
xlabel('FSD T-T Ground State Probability'); 
PrettyFigureFormat
export_fig(gcf,savefile,'-m3');

%% Display m2 error matrix
for i = 1:numel(rhod_local)
    for j = 1:numel(pgs_local)
        m2err_m(i,j) = rhod_pgs_struc(i,j).err(1);
    end
end
range    = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart)-TwinFSDrhoD.ModelObj.Q_i);
qUmin = round(TwinFSDrhoD.ModelObj.qU(TwinFSDrhoD.exclDataStart));
qUmax = round(TwinFSDrhoD.ModelObj.qU(end));
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Fit Grid Scan - %d Runs - [%.0f - %0.f] eV - Stat+Sys - \\sigma_{m^2} (eV^2)',...
    numel(TwinFSDrhoD.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_UniformStackedGridScan_%d_%s_%.0feVbelowE0-m2error.png',...
    numel(TwinFSDrhoD.RunList),TwinFSDrhoD.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit Grid Scan','NumberTitle','off','rend','painters','pos',[10 10 1000 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

cm=0:1/1000:1; my_cm_blue=[1-cm;1-cm;ones(size(cm))]';my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);imagesc(m2err_m); c=colorbar;  c.Label.String = 'eV^2';
set(gca,'XTick',[1:numel(rhod_local)],'XTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},...
    'YTick',[1:numel(pgs_local)],'YTickLabel',{'-2%';'-1%';'nominal';'+1%';'+2%'},'YDir','normal');
[x,y] = meshgrid(1:numel(rhod_local),1:numel(pgs_local));
text(x(:),y(:),num2str(m2err_m(:)),'HorizontalAlignment','center');
ylabel('Column Density'); 
xlabel('FSD T-T Ground State Probability'); 
PrettyFigureFormat
export_fig(gcf,savefile,'-m3');
