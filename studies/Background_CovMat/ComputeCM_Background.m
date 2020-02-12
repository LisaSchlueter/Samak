% Compute Uncertainty on slope
function [CovMat, CovMatFrac, CovMatShape, CovMatFracShape] = ComputeCM_Background(varargin)
p=inputParser;
p.addParameter('StudyObject','',@(x)isa(x,'TBD'));
p.addParameter('plotFit','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('savePlot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nTrials',1000,@(x)isfloat(x));
p.addParameter('RecomputeFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
StudyObject = p.Results.StudyObject;
plotFit     = p.Results.plotFit;
savePlot    = p.Results.savePlot;
nTrials     = p.Results.nTrials;
RecomputeFlag = p.Results.RecomputeFlag;

TimeSec = StudyObject.TimeSec;
qU = StudyObject.qU(StudyObject.qU>StudyObject.Q);         % Background points
nqU = numel(qU);                                           % Number of Background points
qUfrac =  StudyObject.qUfrac(StudyObject.qU>StudyObject.Q);% Time fraction spent in background region
BkgRate_asimov = repmat(StudyObject.BKG_RateSec,nqU,1);     % Asimov Background Rate Model (flat)
qULong = StudyObject.qU;      % whole spectrum

cm_path = sprintf('../../../Samak2.0/inputs/Background/CM/');
cm_file = sprintf('BKG_%s_BKGTimeFrac%.3g_%.0fTrials.mat',StudyObject.TD,sum(qUfrac),nTrials);
if exist([cm_path,cm_file],'file') && strcmp(RecomputeFlag,'OFF')
    cm = importdata([cm_path,cm_file]);
    CovMat          = cm.CovMat;
    CovMatFrac      = cm.CovMatFrac;
    CovMatShape     = cm.CovMatShape;
    CovMatFracShape = cm.CovMatFracShape;
else % Compute CovMat
    
    if ~exist(cm_path,'dir')
        mkdir ../../../Samak2.0/inputs/Background/CM/
    end
    
    par = zeros(2,nTrials); % 2 fitted parameter: slope + offset
    err = zeros(2,nTrials);
    chi2min = zeros(nTrials,1);
    BkgRate_Fit = zeros(numel(qULong),nTrials); % Save fitted Background
    
    % Fit to Background
    BkgRateData = BkgRate_asimov +  sqrt(BkgRate_asimov./(qUfrac.*TimeSec)).*randn(nqU,nTrials); % Rate
    for i=1:nTrials
        Data = [qU,BkgRateData(:,i), sqrt(BkgRateData(:,i)./(qUfrac.*TimeSec))]; %Rate
        % init fit
        BKG_i = wmean(BkgRateData(:,i),1./(sqrt(BkgRateData(:,i)./(qUfrac.*TimeSec))).^2);
        Slope_i = 0;
        parInit = [BKG_i, Slope_i];
        tmparg = sprintf(['set pri -10;'...% Minuit Arguments
            'migrad minos'],'');
        Args = {parInit, Data, '-c', tmparg};
        [par(:,i), err(:,i),chi2min(i), ~] = fminuit('chi2BKG',Args{:});
        BkgRate_Fit(:,i) =  ComputeBkgSlope(par(:,i),qULong); % Rate
    end
    
    %Compute Covariance Matrix
    BkgRate_Fit = BkgRate_Fit(:,chi2min<6);
   % BkgRate_Fit = BkgRate_Fit(:,all(BkgRate_Fit>0,1)); %Background has to be positive
    Bkg_Fit = BkgRate_Fit.*(TimeSec.*StudyObject.qUfrac); % Covariance Matrix with counts (not rate)
    CovMat = cov(Bkg_Fit');
    BKG_asimov_long = StudyObject.BKG_RateSec.*TimeSec.*StudyObject.qUfrac;% flat background whole range
    %CovMatFrac = CovMat./BKG_asimov_long./BKG_asimov_long'
    %TEST
    CovMatFrac = bsxfun(@rdivide,CovMat,BKG_asimov_long);      %divides 1st row of CovMat by TBDIS(1), etc...
    CovMatFrac = bsxfun(@rdivide,CovMatFrac,BKG_asimov_long'); %divides columnwise
    %END TEST
    
    %Compute Shape Only
    BkgIS_Sample = mvnrnd(BKG_asimov_long,CovMat,nTrials)';
    BkgIS_Sample = BkgIS_Sample(:,all(BkgIS_Sample>0,1));
    BkgIS_SumExpected = sum(BKG_asimov_long);
    BkgIS_SumSample = sum(BkgIS_Sample,2);
    BkgIS_SampleNorm = BkgIS_Sample.*(BkgIS_SumExpected./BkgIS_SumSample);
    CovMatShape = cov(BkgIS_SampleNorm');
    %CovMatFracShape = CovMatShape./BKG_asimov_long./BKG_asimov_long';
    %TEST
    CovMatFracShape = bsxfun(@rdivide,CovMatShape,BKG_asimov_long);      %divides 1st row of CovMat by TBDIS(1), etc...
    CovMatFracShape = bsxfun(@rdivide,CovMatFracShape,BKG_asimov_long'); %divides columnwise
    %END TEST
    
    save([cm_path,cm_file],'CovMat','CovMatFrac','CovMatShape','CovMatFracShape',...
        'qUfrac','qU','TimeSec','BkgRate_asimov','Bkg_Fit','BkgRate_Fit','BkgRateData','par','err','chi2min');
    %% plot
    if strcmp(plotFit,'ON')
        % Background (Rate) Spectra
        f21 =  figure('Renderer','opengl');
        set(f21, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
        pdata = plot(qU-18575,1e3.*BkgRateData(:,1:100),'x');
        hold on
        c = {'GoldenRod','IndianRed','CadetBlue'};
        for i=2:1:100
%             pfit = plot(qULong-18575,1e3*BkgRate_Fit(:,i),'Color',rgb(c{i}),'LineWidth',5); %hold on;
%             pdata = errorbar(qU-18575,1e3.*BkgRateData(:,i),1e3*sqrt(BkgRateData(:,i)./(TimeSec.*qUfrac(i))),...
%                 'o','MarkerSize',12,'Color',rgb(c{i}),'MarkerFaceColor',rgb(c{i}),'MarkerEdgeColor',rgb(c{i}),'LineWidth',2);
            pfit = plot(qULong-18575,1e3*BkgRate_Fit(:,i),'LineWidth',5); %hold on;
            pdata = errorbar(qU-18575,1e3.*BkgRateData(:,i),1e3*sqrt(BkgRateData(:,i)./(TimeSec.*qUfrac(1))),...
                'o','MarkerSize',12,'LineWidth',2);
        end
        %xlim([-2 5])
        PrettyFigureFormat;
        set(gca,'FontSize',22);
        xlabel('retarding potential qU - 18575 (eV)');
        ylabel('background (mcps)');
        hold off;
        legend([pdata, pfit],'MC data','linear fit','Location','northwest');
        legend boxoff
        fig21_str = sprintf('CM_BKGFits_%s',StudyObject.TD);
        publish_figurePDF(f21,['./plots/CovMatInfo/pdf/',fig21_str,'.pdf']);
        print(f21,['./plots/CovMatInfo/png/',fig21_str,'.png'],'-dpng');
        savefig(f21,['./plots/CovMatInfo/fig/',fig21_str,'.fig'],'compact');
        %% CovMat
        %%
        f22 = figure('Renderer','opengl');
        set(f22, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.75]);
        subplot(2,2,[1,3]);
        %imagesc(CovMatFracShape);
        imagesc(CovMat);
        PrettyFigureFormat
        c = colorbar('northoutside');
        c.Label.String = 'fractional covariance (bkg)';
        c.FontSize = 18;
        colormap(parula);
        pbaspect([1 1 1])
        set(gca,'xtick',[1 StudyObject.nqU]),set(gca,'ytick',[])
        qUmin = sprintf('qU_{min} = E_0- %.0fV',abs(StudyObject.qU(1)-18575));
        qUmax = sprintf('qU_{max} = E_0+ %.0fV',StudyObject.qU(end)-18575);
        set(gca,'xticklabel',{qUmin,qUmax}); set(gca,'yticklabel',[]);
        set(gca,'FontSize',16)      
        
        subplot(2,2,2);
       [SysLine, SysArea] =  boundedline(qULong-18575,1e3.*mean(BkgRate_Fit,2),1e3.*std(permute(BkgRate_Fit,[2,1])));
       hold on;
       pStat = errorbar(qULong-18575,1e3.*mean(BkgRate_Fit,2),1e3.*sqrt(mean(BkgRate_Fit,2)./(TimeSec.*StudyObject.qUfrac)),...
           'o','Color',rgb('GoldenRod'),'MarkerFaceColor',rgb('RoyalBlue'),'MarkerEdgeColor',rgb('RoyalBlue'),'MarkerSize',2);
       SysArea.FaceColor = rgb('CadetBlue');
        SysLine.LineWidth = 5; SysLine.Color = rgb('RoyalBlue');
        SysArea.FaceAlpha = 0.5;
        xlabel('retarding potential qU - E_0 (V)');
        ylabel('background (mcps)');
        leg =legend([SysLine, SysArea,pStat],'mean','1 \sigma_{sys}','1 \sigma_{stat}','Location','northeast');
        leg.FontSize = 16;
        leg.NumColumns = 3;
        legend boxoff
        xlim([min(qULong) max(qULong)]-18575);     
        ylimdown = (1e3*mean(mean(BkgRate_Fit,2)))-1e3*min(sqrt(mean(BkgRate_Fit,2)./(TimeSec.*StudyObject.qUfrac)));%1e3*(mean(mean(BkgRate_Fit,2))-max(std(permute(BkgRate_Fit,[2,1]))));
        ylimup = 1e3*mean(mean(BkgRate_Fit,2))+(1e3.*max(sqrt(mean(BkgRate_Fit,2)./(TimeSec.*StudyObject.qUfrac))));%1e3*(mean(mean(BkgRate_Fit,2))+max(std(permute(BkgRate_Fit,[2,1]))));
        if ylimup<=1e3*max(max(mean(BkgRate_Fit,2)+std(permute(BkgRate_Fit,[2,1]))))
            ylimup =1e3*max(max(mean(BkgRate_Fit,2)+std(permute(BkgRate_Fit,[2,1]))));
            ylimdown =1e3*min(min(mean(BkgRate_Fit,2)-std(permute(BkgRate_Fit,[2,1]))));
        end
        ylim([0.99*ylimdown, 1.01*ylimup]);
        PrettyFigureFormat;
        set(gca,'FontSize',16);
        
        subplot(2,2,4);
        Nmax = max(size(BkgIS_SampleNorm));
        StepSize = 10;
        CovMatnTrials = cell(round(Nmax/StepSize),1);     % CovMat for different nTrials
        CovMatTrace = zeros(round(Nmax/StepSize),1);      % Trace of CovMats
        for x = StepSize:StepSize:Nmax %Compute CovMats für different nTrials
            xn = x/StepSize;
            CovMatnTrials{xn,1} = cov(BkgIS_SampleNorm(:,1:x)'); %CovMats
            CovMatTrace(xn) = norm(CovMatnTrials{xn,1},'fro');
        end
        
        while max(size([StepSize:StepSize:Nmax]))~=max(size(CovMatTrace))          
        if max(size([StepSize:StepSize:Nmax]))>=max(size(CovMatTrace))
            Nmax = Nmax-StepSize;  
        elseif max(size([StepSize:StepSize:Nmax]))<=max(size(CovMatTrace))
            Nmax = Nmax+StepSize;    
        end
        end
        
        plot([StepSize:StepSize:Nmax]',CovMatTrace,'Color',rgb('GoldenRod'),'LineWidth',4);
        xlabel('samples');
        ylabel('|| M ||');
        PrettyFigureFormat;
        set(gca,'FontSize',16);
        xlim([0 Nmax-2*StepSize]);
       
        
        f23 = figure('Renderer','opengl');
        set(f23, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,0.5]);
        plot([StepSize:StepSize:Nmax]',CovMatTrace,'Color',rgb('GoldenRod'),'LineWidth',4);
        xlabel('samples');
        ylabel('|| M ||');
        PrettyFigureFormat;
        set(gca,'FontSize',22);
        xlim([0 Nmax-2*StepSize]);
        
        f23_str = [sprintf('Convergence_%s_Bkg',StudyObject.TD)];
        publish_figurePDF(f23,['./plots/CovMatInfo/pdf/',f23_str,'.pdf']);
        print(f23,['./plots/CovMatInfo/png/',f23_str,'.png'],'-dpng');
        savefig(f23,['./plots/CovMatInfo/fig/',f23_str,'.fig'],'compact');
    end
    if strcmp(savePlot,'ON')
        fig100_str = sprintf('CM_Background_%s',StudyObject.TD);
        publish_figurePDF(f22,['./plots/CovMatInfo/pdf/',fig100_str,'.pdf']);
        print(f22,['./plots/CovMatInfo/png/',fig100_str,'.png'],'-dpng');
        savefig(f22,['./plots/CovMatInfo/fig/',fig100_str,'.fig'],'compact');
    end
end
end


