% position of profile best fits
% based on ksn2_NuMassSensitivity.m
% impact of 3+1 model extionsion on neutrino mass sensitivity
% plus nu-mass isoline 
chi2 = 'chi2CMShape';
DataType ='Twin';% 'Real';%
SavePlt = 'ON';

CombiStereo = 'OFF';

if strcmp(DataType,'Real')
mNuSq = -1:0.01:2.5;
else
mNuSq = -1:0.01:2.5;
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
if strcmp(CombiStereo,'OFF')
    CombiStr = '';
else
    CombiStr = 'CombiStereo';
end
savefile = sprintf('%sksn2_NuMassSensitivity%s_%s_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.mat',...
    savedir,CombiStr,DataType,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));%CombiStereo

if exist(savefile,'file') 
    load(savefile);
     fprintf('load file  %s \n',savefile)
else
     fprintf('file does not exist: %s \n',savefile)
    return
end


GetFigure;


mNuSq_leg = [-1:0.25:0,0.3,0.5,1,2];
Mycolors= [autumn(sum(mNuSq_leg<=0));winter(sum(mNuSq_leg>0))];

pHandles = cell(numel(mNuSq_leg),1);
legStr   = cell(numel(mNuSq_leg),1);


for i=1:numel(mNuSq)
    pchi2min =  plot(sin2T4_bf(i),mNu4Sq_bf(i),'.','MarkerSize',12,'Color',rgb('Orange'));
    hold on;
 LineWidthIso = 1;
    if max(diff(mNu4Sq_iso{i}))>100
        mNu4Sq_iso_tmp = mNu4Sq_iso{i};
        sin2T4_iso_tmp = sin2T4_iso{i};
        nContourIdx = find(abs(diff(mNu4Sq_iso{i}))>=100);
        if numel(nContourIdx)==1
            plot(sin2T4_iso_tmp(1:nContourIdx),mNu4Sq_iso_tmp(1:nContourIdx),'-','Color',rgb('Silver'),'LineWidth',LineWidthIso);
            piso =  plot(sin2T4_iso_tmp(nContourIdx+1:end),mNu4Sq_iso_tmp(nContourIdx+1:end),'-','Color',rgb('Silver'),'LineWidth',LineWidthIso);
        else
             StartIdx = 1;
             StopIdx = nContourIdx(1);
 
            for j=1:numel(nContourIdx)
                plot(sin2T4_iso_tmp(StartIdx:StopIdx-1),mNu4Sq_iso_tmp(StartIdx:StopIdx-1),'-','Color',rgb('Silver'),'LineWidth',LineWidthIso);
                StartIdx = nContourIdx(j)+1;
                if j==numel(nContourIdx)
                    StopIdx = numel(sin2T4_iso_tmp);   
                else
                    StopIdx = nContourIdx(j+1);
                end
            end 

        end
    else
        piso =  plot(sin2T4_iso{i},mNu4Sq_iso{i},'-','Color',rgb('Silver'),'LineWidth',LineWidthIso);
    end

end
%
plot(sin2T4_contour_bf,mNu4Sq_contour_bf,'x','LineWidth',2,'Color',rgb('ForestGreen'));
pfree = plot(sin2T4_contour,mNu4Sq_contour,'-','LineWidth',2,'Color',rgb('ForestGreen'));
PrettyFigureFormat('FontSize',22);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_{4}^2 (eV^{ 2})'));
%
leg = legend([piso,pfree,pchi2min],sprintf('{\\itm}_\\nu^2 best fit isoline'),...
  sprintf('Exclusion curve: {\\itm}_\\nu^2 unconstrained'),...
  sprintf('Smallest \\chi^2 on {\\itm}_\\nu^2-isoline'),...
  'Location','southwest');
%leg.Title.String = sprintf('Local best fit for {\\itm}_\\nu^2 (eV^2)');
%leg.Title.FontSize = get(gca,'FontSize'); leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);
xlim([1e-03 0.5]);
ylim([0.1 40^2]);

%
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/plots/'];
MakeDir(pltdir)
pltfile = sprintf('%sksn2_NuMassSensitivity%s_%s_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.pdf',...
    pltdir,CombiStr,DataType,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));%CombiStereo
export_fig(pltfile)
%print(pltfile,'-dpng','-r350');
