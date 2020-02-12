% fit  with different energy loss functions
% plot impact on neutrino mass

% config
DataType ='Twin';
RunList = 'KNM1';
exclDataStart = 14;
chi2 = 'chi2Stat';
ReFit = 'OFF';

%label results and load if possible
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('NuMassShift_ELoss_%s_%s_%.0f.mat',DataType,RunList,exclDataStart)];
if exist(savename,'file') && strcmp(ReFit,'OFF')
    load(savename)
else 
    % init models with different eloss functions
    CommongArg = {'RunList',RunList,'chi2','chi2Stat','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
        'exclDataStart',exclDataStart,'chi2',chi2,'minuitOpt','min;minos'};
    M_T2 = MultiRunAnalysis(CommongArg{:},'ELossFlag','KatrinT2');
    M_D2 = MultiRunAnalysis(CommongArg{:},'ELossFlag','KatrinD2');
    M_Abd = MultiRunAnalysis(CommongArg{:},'ELossFlag','Abdurashitov');
   % M_Ase = MultiRunAnalysis(CommongArg{:},'ELossFlag','Aseev');
    
    %  statfit
    M_T2.Fit;
    M_D2.Fit;
    M_Abd.Fit;
    mNuSq = zeros(3,1); mNuSqErr = zeros(3,1);
    mNuSq(1) = M_T2.FitResult.par(1);  mNuSqErr(1) = M_T2.FitResult.err(1);
    mNuSq(2) = M_D2.FitResult.par(1);  mNuSqErr(2) = M_D2.FitResult.err(1);
    mNuSq(3) = M_Abd.FitResult.par(1); mNuSqErr(3) = M_Abd.FitResult.err(1);
   
    % CM fit
    M_T2.chi2 = 'chi2CMShape';
    M_T2.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF');
    M_T2.Fit;
    
    M_D2.chi2 = 'chi2CMShape';
    M_D2.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF');
    M_D2.Fit;
    
    M_Abd.chi2 = 'chi2CMShape';
    M_Abd.ComputeCM('SysEffects',struct('RF_EL','ON'),'BkgCM','OFF');
    M_Abd.Fit;

    mNuSqCM = zeros(3,1); mNuSqCMErr = zeros(3,1);
    mNuSqCM(1) = M_T2.FitResult.par(1);  mNuSqCMErr(1) = M_T2.FitResult.err(1);
    mNuSqCM(2) = M_D2.FitResult.par(1);  mNuSqCMErr(2) = M_D2.FitResult.err(1);
    mNuSqCM(3) = M_Abd.FitResult.par(1); mNuSqCMErr(3) = M_Abd.FitResult.err(1);
  
    range = abs(round(M_T2.RunData.qU(exclDataStart)-18575));
    eloss_str = {sprintf('Katrin T_2'),sprintf('Katrin D_2'),'Abdurashitov'};%,'Aseev'};
    save(savename,'mNuSq','mNuSqErr','mNuSqCM','mNuSqCMErr','eloss_str','range');
end

%% plot neutrino mass shift as a function of energy loss function
f22 = figure('Renderer','opengl');
set(f22, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

x= 1:numel(mNuSq);
%plot(x,mNuSq,'--o','LineWidth',3,'MarkerSize',8,...%mNuSqErr,
 %   'MarkerFaceColor',rgb('CadetBlue'),'Color',rgb('SteelBlue'));
%hold on;
plot(x,mNuSqCM,'--o','LineWidth',3,'MarkerSize',8,...%mNuSqCMErr,
    'MarkerFaceColor',rgb('IndianRed'),'Color',rgb('FireBrick'));

%plot(x,wmean(mNuSq,mNuSqErr).*ones(numel(x),1));
xticks(x); xticklabels(eloss_str);
ylabel(sprintf('m^2_\\nu (eV^2)'));
PrettyFigureFormat;
grid on;
plotname = strrep(strrep(savename,'results','plots'),'.mat','.png');
leg = legend('stat + sys'); legend boxoff;%,'stat + sys (EL)'); legend boxoff;
leg.Title.String =  sprintf('%.0feV range',range);
leg.Location ='northwest';leg.NumColumns=2;
set(gca,'FontSize',22);
xlim([min(x) max(x)]);%-0.1,+0.2
%ylim([-0.9 1.2]);
print(f22,plotname,'-dpng','-r450');

