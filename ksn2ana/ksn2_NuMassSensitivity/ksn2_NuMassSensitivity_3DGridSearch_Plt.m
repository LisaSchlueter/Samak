% ksn2 calculate chi2-prifle for nu-mass
% perform grid searches for differed nu-masses

%% settings that might change
DataType = 'Real';
switch DataType
    case 'Twin'
        FixmNuSq_all = -1:0.1:1;
    case 'Real'
        FixmNuSq_all = sort(round([-1.0:0.1:2.5],2));
end

RecomputeFlag = 'ON';
Maxm4Sq    = 36^2; % interpolation
chi2min = zeros(numel(FixmNuSq_all),1);
sin2T4_bf =  zeros(numel(FixmNuSq_all),1);
mNu4Sq_bf =  zeros(numel(FixmNuSq_all),1);
CombiSTEREO = 'ON';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefileCombi = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_CombiSTEREO-%s.mat',...
    savedir,DataType,CombiSTEREO);

if exist(savefileCombi,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefileCombi);
else
    %% configure RunAnalysis object
    
    range = 40;
    chi2 = 'chi2CMShape';
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %%
    switch DataType
        case 'Real'
            LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON','CheckSmallerN','ON'};
            nGridSteps = 40;
        case 'Twin'
            LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
            nGridSteps = 30;
    end
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',LoadGridArg};
    S = SterileAnalysis(SterileArg{:});
    %%
    if strcmp(CombiSTEREO,'ON')
        % load STEREO
        savedir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/'];
        savefile = sprintf('%sksn2_InterpStereoMini_Max%.0feV2_Min%.2gfeV2.mat',...
            savedir,Maxm4Sq,0.1);
        if exist(savefile,'file')
            fprintf('file exist %s \n',savefile);
            dS = importdata(savefile);
        else
            return
        end
        chi2min_Combi   = zeros(numel(FixmNuSq_all),1);
        sin2T4_Combi_bf =  zeros(numel(FixmNuSq_all),1);
        mNu4Sq_Combi_bf =  zeros(numel(FixmNuSq_all),1);
    end
    
    %%
    for i=1:numel(FixmNuSq_all)
        FixmNuSq = FixmNuSq_all(i);
        A.ModelObj.mnuSq_i = FixmNuSq;
        S.LoadGridFile(S.LoadGridArg{:},'FixmNuSq',FixmNuSq);
        S.Interp1Grid('MaxM4Sq',Maxm4Sq);
        if strcmp(CombiSTEREO,'OFF')
            S.GridPlot('BestFit','ON','Contour','ON','SavePlot','png',...
                'ExtraStr',sprintf('_mNuSq%.2geV2',FixmNuSq)); close;
        else
            S.FindBestFit;
        end
        
        if S.chi2_bf<S.chi2_Null
            chi2min(i) = S.chi2_bf;
            mNu4Sq_bf(i) = S.mNu4Sq_bf;
            sin2T4_bf(i) = S.sin2T4_bf;
        else
            chi2min(i) = S.chi2_Null;
            mNu4Sq_bf(i) = 0;
            sin2T4_bf(i) = 0;
        end

        if strcmp(CombiSTEREO,'ON')
            chi2_KATRIN = S.chi2;
%             Deltachi2_Katrin = S.chi2-S.chi2_ref;
%             Deltachi2_Stereo = d.chi2Stereo-min(min(d.chi2Stereo_cut));
%             Deltachi2_Stereo(Deltachi2_Stereo<0) = 0;
%             Deltachi2Sum = Deltachi2_Katrin+Deltachi2_Stereo-min(min(Deltachi2Sum));
%             
            %% build shared chi^2 map
            % 2 regions:
            %1) KATRIN + STEREO
            Chi2Stereo = dS.chi2Stereo_cut;%-min(min(dStereo.chi2Stereo_cut));
            Chi2KATRIN  =  chi2_KATRIN(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol);
            Chi2Sum = Chi2Stereo + Chi2KATRIN;
            DeltaChi2Combi_Shared = Chi2Sum - min(min(Chi2Sum));
            
            % 2) KATRIN only
            chi2min_Combi(i) = min(min(chi2_KATRIN(~dS.InterIdx)));
            DeltaChi2Combi_KATRINonly = chi2_KATRIN-min(min(chi2_KATRIN(~dS.InterIdx)));
            
            % Combine 1)+2)
            DeltaChi2Combi = zeros(1e3);
            DeltaChi2Combi(dS.InterIdx) = DeltaChi2Combi_Shared+1e-05;
            DeltaChi2Combi(~dS.InterIdx) = DeltaChi2Combi_KATRINonly(~dS.InterIdx);
            
            %%
            % combine with STEREO
%             S.chi2_ref = chi2min(i);
%             Deltachi2Sum = d.chi2Stereo+chi2_Katrin;%-S.chi2_ref; % combined Delta-Chi^2 Maps
            S.chi2_ref = min(min(DeltaChi2Combi)); % 0 by construction
            S.chi2 = DeltaChi2Combi;%Deltachi2Sum;
            S.GridPlot('BestFit','ON','Contour','ON','CL',99.99,...
                'SavePlot','OFF',...
                'ExtraStr',sprintf('STEREOCombi_mNuSq%.2geV2',FixmNuSq));
            close;
       %     chi2min_Combi(i) = S.chi2_bf;
            mNu4Sq_Combi_bf(i) = S.mNu4Sq_bf;
            sin2T4_Combi_bf(i) = S.sin2T4_bf;
            %%
            
        end
    end
    save(savefileCombi,'chi2min','mNu4Sq_bf','sin2T4_bf','FixmNuSq_all')
    if strcmp(CombiSTEREO,'ON')
        save(savefileCombi,'chi2min_Combi','mNu4Sq_Combi_bf','sin2T4_Combi_bf','-append')
    end
end
%% some simple calculations
% interpolate
mNuSq_inter = linspace(min(FixmNuSq_all),max(FixmNuSq_all),1e3);
chi2_inter  = interp1(FixmNuSq_all,chi2min,mNuSq_inter,'spline');
% best fit
chi2bf      = min(chi2_inter);
Idx_bf      = find(chi2_inter==chi2bf);
mNuSq_bf    = mNuSq_inter(Idx_bf);
% uncertainties
mNuSqDown = interp1(chi2_inter(mNuSq_inter<mNuSq_bf),mNuSq_inter(mNuSq_inter<mNuSq_bf),chi2bf+1,'spline');
mNuSqUp   = interp1(chi2_inter(mNuSq_inter>mNuSq_bf),mNuSq_inter(mNuSq_inter>mNuSq_bf),chi2bf+1,'spline');
mNuSqErrDown = mNuSq_bf-mNuSqDown;
mNuSqErrUp = mNuSqUp-mNuSq_bf;
mNuSqErr   = 0.5.*(mNuSqErrDown+mNuSqErrUp);

if strcmp(CombiSTEREO,'ON')
    chi2_inter_Combi  = interp1(FixmNuSq_all,chi2min_Combi,mNuSq_inter,'spline');
    % best fit
    chi2bf_Combi      = min(chi2_inter_Combi);
    Idx_bf_Combi      = find(chi2_inter_Combi==chi2bf_Combi);
    mNuSq_bf_Combi    = mNuSq_inter(Idx_bf);
    % uncertainties
    mNuSqDown_Combi = interp1(chi2_inter_Combi(mNuSq_inter<mNuSq_bf_Combi),mNuSq_inter(mNuSq_inter<mNuSq_bf_Combi),chi2bf_Combi+1,'spline');
    mNuSqUp_Combi   = interp1(chi2_inter_Combi(mNuSq_inter>mNuSq_bf_Combi),mNuSq_inter(mNuSq_inter>mNuSq_bf_Combi),chi2bf_Combi+1,'spline');
    mNuSqErrDown_Combi = mNuSq_bf_Combi-mNuSqDown_Combi;
    mNuSqErrUp_Combi   = mNuSqUp_Combi-mNuSq_bf_Combi;
    mNuSqErr_Combi     = 0.5.*(mNuSqErrDown_Combi+mNuSqErrUp_Combi);
end
%% plot
GetFigure;
plot(mNuSq_inter,chi2bf.*ones(numel(mNuSq_inter),1),'k--','LineWidth',1);
hold on;
plot(mNuSq_inter,(chi2bf+1).*ones(numel(mNuSq_inter),1),'k--','LineWidth',1);
p1inter = plot(mNuSq_inter,smooth(chi2_inter,100),'-','LineWidth',3,'Color',rgb('SkyBlue'));
%p1 = plot(FixmNuSq_all,chi2min,'.','MarkerSize',15,'Color',rgb('DodgerBlue'));

if strcmp(CombiSTEREO,'ON')  
    p2Stereo = plot(mNuSq_inter,smooth(chi2_inter_Combi,100),':','LineWidth',3,'Color',rgb('Orange'));
   % p2 = plot(FixmNuSq_all,chi2min_Combi,'x','MarkerSize',10,'LineWidth',2,'Color',rgb('DarkOrange'));
end

dof = 28-4-2;% nqU-freeFitPar
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\chi^2 (%.0f dof)',dof));
PrettyFigureFormat('FontSize',22);

if strcmp(CombiSTEREO,'OFF') 
leg = legend(p1inter,sprintf('{\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2',mNuSq_bf,mNuSqErrUp,mNuSqErrDown),...
    'Location','northwest');
else
    leg = legend([p1inter,p2Stereo],...
        sprintf('{\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2 (KATRIN)',mNuSq_bf,mNuSqErrUp,mNuSqErrDown),...
        sprintf('{\\itm}_\\nu^2 = %.1f^{+%.2f}_{-%.2g} eV^2 (KATRIN + STEREO)',mNuSq_bf_Combi,mNuSqErrUp_Combi,mNuSqErrDown_Combi),...
     'Location','north');
end

PrettyLegendFormat(leg);

if strcmp(DataType,'Twin')
    t =  title('MC Twin');
    ylim([chi2bf-0.5 chi2bf+2]);

else
    t =  title('Data');
    ylim([chi2bf-0.3 chi2bf+4]);
end
t.FontWeight = 'normal'; t.FontSize = get(gca,'FontSize');

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep(strrep(savefileCombi,'results','plots'),'.mat','.png');
print(gcf,pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);
