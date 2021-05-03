% plot significant contours (delta chi^2>5.99)
Hypothesis = 'H0';
InterpMode = 'lin';
SavePlt = 'ON';
switch Hypothesis
    case 'H0'
        randMC =[1001:1250,1358:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1' 
        randMC = 1:1e3;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples.mat',savedir,InterpMode,numel(randMC));
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC));
end

if exist(savefile,'file')
    fprintf('load file %s \n',savefile);
    d = importdata(savefile);
end
%
IdxSig = find(d.ClosedLog95); % find contours with significance
fprintf('%.0f out of %.0f  (%.1f%%) with delta chi^2 > 5.99 \n',numel(IdxSig),numel(randMC),100*numel(IdxSig)/numel(randMC));
%%
tmpPltNo = 1;
f1 = figure('Units','normalized','Position',[0.1,0.1,1,1]);
for i=1:numel(IdxSig)
    if mod(i,10)==0
        if strcmp(SavePlt,'ON')
            MakeDir('./tmp/');
            export_fig(gcf,sprintf('tmp/TmpPlt%.0f.pdf',tmpPltNo));
        end
        f1 = figure('Units','normalized','Position',[0.1,0.1,1,1]);
        tmpPltNo = tmpPltNo+1;
        SubPltIdx = 1;
    elseif mod(i,9)==0
        SubPltIdx = 9;
    else
        SubPltIdx = mod(i,9);
    end
    
    subplot(3,3,SubPltIdx);
    pA = plot(d.sin2T4_contour_Asimov,d.mNu4Sq_contour_Asimov,'k-','LineWidth',2);
    hold on;
    pR = plot(d.sin2T4_contour{IdxSig(i)},d.mNu4Sq_contour{IdxSig(i)},'Color',rgb('Orange'),'LineWidth',2);
    pRbf=plot(d.sin2T4_bf(IdxSig(i)),d.mNu4Sq_bf(IdxSig(i)),'x','LineWidth',2,'Color',rgb('Orange'));
    leg = legend(pR(1),sprintf('RandMC %.0f (|{\\itU}_{e4}|^2=%.3g, {\\itm}_4^2 = %.1f eV^2)',randMC(IdxSig(i)),d.sin2T4_bf(IdxSig(i)),d.mNu4Sq_bf(IdxSig(i))),...
        'Location','southwest'); 
    PrettyLegendFormat(leg);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
    PrettyFigureFormat('FontSize',15);
    hold off;
    xlim([1e-03 0.5]);
    ylim([0.1,1600]);
    if i==numel(IdxSig) && strcmp(SavePlt,'ON')
        export_fig(gcf,sprintf('tmp/TmpPlt%.0f.pdf',tmpPltNo));
    end
end

if strcmp(SavePlt,'ON')
    pltname = strrep(strrep(savefile,'results','plots'),'.mat','_SigContours.pdf');
    system(sprintf('"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o %s tmp/*.pdf',pltname));
    system('rm -r ./tmp/');
    close all
end