function [par err chi2min ndof] = KrLineFitallPixels(varargin)
%
%            FIT INTEGRAL SPECTRUM
%              Fit a 83mKr Lines
%          All Pixels Simultaneously
%
%  For each pixels
%  - Background fixed with pixel-wise fit : 148 par. fixed
%  - Amplitude fixed with pixel-wise fit : 148 par. fixed
%  A common fit for: E0 and W
%
%          Th. Lasserre - CEA Saclay
%                October 2017
%

    % Initialization
    clear par ; clear parnomix;
    addpath('fminuit'); addpath('krypton-data'); addpath('tools'); %digits(6);

    % Parser
    p = inputParser;
    p.addParameter('fign',1,@(x)isfloat(x) && x>0);
    p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('Mode','Data',@(x)ismember(x,{'Data','Sim'}));
    p.addParameter('CPS','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'PIXEL'}));
    p.addParameter('Pixel',0,@(x)isfloat(x) && x==0);
    p.addParameter('nPixels',111,@(x)isfloat(x) && x>0);
    p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
    p.addParameter('TD','KrL3_32',@(x)ismember(x,{'KrK32','KrL3_32','KrL3_32_HS'}));
    p.addParameter('DopplerEffectFlag','Conv',@(x)ismember(x,{'OFF','Conv','Voigt','realConv'}));
    p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
    p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV

p.parse(varargin{:});
    fign       =    p.Results.fign;
    pub        =    p.Results.pub;
    display    =    p.Results.display;
    Mode       =    p.Results.Mode;
    CPS        =    p.Results.CPS;
    FPD_Segmentation  = p.Results.FPD_Segmentation;
    Pixel      =    p.Results.Pixel;
    nPixels    =    p.Results.nPixels;
    nfit       =    p.Results.nfit;
    TD         =    p.Results.TD;
    DopplerEffectFlag = p.Results.DopplerEffectFlag;
    HVRipples  = p.Results.HVRipples;
    HVRipplesP2PV = p.Results.HVRipplesP2PV;

    % Parametrization: True Value
    global A ; 
    switch TD
        case {'KrK32','KrL3_32','KrL3_32_HS'}
            A=InitKrKATRIN_AllPixels(...
                'TD',TD,'nPixels',nPixels,'FPD_Pixel',Pixel,...
                'DopplerEffectFlag',DopplerEffectFlag,...
                'CPS',CPS,'FPD_Segmentation',FPD_Segmentation,...
                'HVRipples',HVRipples,'HVRipplesP2PV',HVRipplesP2PV);
    end
    
% Initialization with Pixel-by-Pixel Fit
    switch Mode
        case 'Data'
            switch TD
                case 'KrK32'
                    qumin=1;
                    krinit = importdata('Data_KrK32AllPixelsDopplerONRipplesON.mat');
                    A.SetLineInitParallPixels();
                    fprintf(2,'Initialize Line Papameters - KrK32 - %.0f Pixels\n',nPixels);
                    %disp(A.K32_Phi0allPixels);
                case 'KrL3_32'
                    qumin=1;
                    krinit = importdata('Data_KrL332AllPixelsDopplerONRipplesON.mat');
                    A.SetLineInitParallPixels();
                    fprintf(2,'Initialize Line Parameters - KrL3_32 - %.0f Pixels\n',nPixels);
                    %disp(A.L3_32_Phi0allPixels);
                case 'KrL3_32_HS'
                    qumin=100;
                    krinit = importdata('Data_KrL332HSAllPixelsDopplerONRipplesON_v3.mat');
                    A.SetLineInitParallPixels();
                    fprintf(2,'Initialize Line Papameters - KrL3_32_HS - %.0f Pixels\n',nPixels);
                    %disp(A.L3_32_Phi0allPixels);
            end
    end

    switch TD
       case 'KrK32'
            LineE_i    = A.K32_E_i;
            LineW_i    = A.K32_W_i;
        case {'KrL3_32','KrL3_32_HS'}
            LineE_i    = A.L3_32_E_i;
            LineW_i    = A.L3_32_W_i;
    end
    
    % Loop on fits
    npar       = 3;
    fit_p      = ones(npar,nfit);
    fit_xi2min = ones(1,nfit);
    % Init
    i_E        = 0; i_W        = 0; i_N = 0;
    
    % Data / Sim
    switch Mode
        case 'Data'
            [qU krdatamap] = readKrData('TD',TD);
    end
    
    % Number of Fits
    tic
    progressbar('Krypton Line Fits');
    for f=1:1:nfit
        progressbar(f/nfit);
        
        switch Mode
            case 'Data'
                % reshape(repmat(A.qU,1,2),[],1)
                % repmat(A.qU,A.nPixels,1)
                Data = [reshape(repmat(A.qU,1,A.nPixels),[],1),...
                    reshape(permute(krdatamap(1:A.nPixels,1,qumin:end),[3 2 1]),[],1),...
                    1.*reshape(permute(krdatamap(1:A.nPixels,2,qumin:end),[3 2 1]),[],1)];
            case 'Sim'
                A.ComputeKrDSallPixels(); A.ComputeKrISallPixels(); %A.AddStatFluctKrISallPixels();
                Data = [repmat(A.qU,A.nPixels,1),...
                    reshape(cell2mat(A.KrISallPixels),[],1),...
                    reshape(cell2mat(A.KrISEallPixels),[],1)];
        end
                
        % Initializing Fit
        ParIni = [i_E i_W i_N];
        parnames = ['E W N'];
        tmparg = sprintf(['set pri -10 ; fix 3'...
            'set now; min ; imp' ]);
        Args = {ParIni, Data, '-c',tmparg};
        [par, err, chi2min, errmat] = fminuit('KrLineChi2allPixels',Args{:}); 
        fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
        
        switch display
            case 'ON'
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf('  Processing \t= %g %% \n',f/nfit*100);
                fprintf(2,'--------------------------------------------------------------\n');
                fprintf(2,'  E \t= %.3f \t± \t %g \t eV \n',  (LineE_i+par(1)),err(1));
                fprintf(2,'  W \t= %.3f \t± \t %g \t eV \n', (LineW_i+par(2)),err(2));
                fprintf(2,'  N \t= %.3f \t± \t %g \t  \n', par(3),err(2));
                ndof = A.nqU*A.nPixels-3;
                fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
                fprintf(2,'--------------------------------------------------------------\n');
        end
    end
    toc
        
    %% Plot Results
    figure(fign+999)
    subplot(2,1,1)
    hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2),Data(1:A.nqU*A.nPixels,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    hold on
    %     hfit1 = plot(A.Te,KrLineModelDiff4par(par)./trapz(A.Te,KrLineModelDiff4par(par)).*trapz(A.qU,KrLineModel4par(par))/4,...
    %         'LineWidth',1,'LineStyle','--','Color','Black');
    model = KrLineModelallPixels(par);
    
    for (i=0:1:A.nPixels-1)
        hfit1 = plot(Data(i*A.nqU+1:i*A.nqU+A.nqU,1),model(i*A.nqU+1:i*A.nqU+A.nqU),...
            'Color','Red','LineWidth',1,'LineStyle','-');
    end
    
    hold off
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Counts','FontSize',10);
    switch Mode
        case {'DATA','Sim'}
            title(sprintf('KATRIN Krypton gas - %s - %s - Pixel 1 to %.0f - Doppler %s',A.TD,Mode,Pixel,A.DopplerEffectFlag));
    end
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',LineE_i,LineW_i);
    myfit = sprintf('Fit: E=%.3f \\pm %.3f eV - W=%.3f \\pm %.3f eV',LineE_i+par(1),err(1),LineW_i+par(2),err(2));
    mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU*A.nPixels-3);
    legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
    axis([min(A.qU) max(A.qU)+1 0.7*min(Data(:,2)) max(Data(:,2))*1.2])
    
    subplot(2,1,2)
    hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels),Data(1:A.nqU*A.nPixels,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    hold on
    line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Red');
    hold off
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Residuals','FontSize',10);
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    axis([min(A.qU) max(A.qU)+1 min(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2 max(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2])

    return;

    pub = 'YES';
    switch pub
        case 'YES'
            switch Mode
                case {'DATA','Sim'}
                    myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
                        TD,Pixel);
                    fprintf(2,'publish: %s\n',myname);publish_figure(999+fign,myname);
            end
            
            figure(fign+1000)
            [h cor] = corplot(errmat(1:4,1:4));
            xticklabels({'Position','Width','Amplitude','Background'})
            yticklabels({'Position','Width','Amplitude','Background'})
            title(sprintf('KATRIN Krypton - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
            fprintf(2,'Error Matrix \n');
            disp(cor);
            
            switch Mode
                case {'DATA','Sim'}
                    myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
                        TD,Pixel);
                    fprintf(2,'publish: %s\n',myname);publish_figure(999+fign,myname);
            end
                        
    if nfit<5
        % For Saving:
        par(1) = LineE_i         + par(1);
        par(2) = LineW_i         + par(2);
    end
            
return;

    figure(fign)
    subplot(2,1,1)
    hfit1 = plot(Data(:,1),KrLineModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hold on;
    hfit2 = plot(Data(:,1),KrLineModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    %errorbar_tick(hdata,200);
    hold off
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Counts','FontSize',10);
    title(sprintf('KATRIN Integral - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i));
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    mydata = sprintf('Data: E=%.2f eV \n',LineE_i);
    myfit = sprintf('Fit: E=%.2f \\pm %.2f eV',LineE_i+par(1),err(1));
    mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
    legend([hdata hfit1 hfit2],mydata,myfit,mychi2,'Location','NorthEast') ;% legend(a,'boxoff');
    axis([min(A.qU) max(A.qU)+1 min(Data(:,2)) max(Data(:,2))*1.1])
    
    subplot(2,1,2)
    parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = 0; parnomix(2) = 0;
    hfit = plot(Data(:,1),...
        KrLineModel4par(par)./KrLineModel4par(parnomix)-1,...
        'Color','Black','LineWidth',1,'LineStyle','-');
    hold on;
    hdata = errorbar(Data(:,1),...
        Data(:,2)./KrLineModel4par(parnomix)-1,...
        Data(:,3)./KrLineModel4par(parnomix),...
        'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
    %errorbar_tick(hdata,200);
    hold off;
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Spectral distorsion','FontSize',10);
    set(gca,'FontSize',12);
    a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','NorthWest');
    %legend(a,'boxoff');
    %axis([min(A.qU) max(A.qU)+1 min((Data(:,2)./KrLineModel4par(parnomix)-1)) max(Data(:,3)./KrLineModel4par(parnomix))*3])

    switch pub
        case 'YES'
           myname = sprintf('./figures/f1-krkatrin_%s_pixel%.0f.eps',...
                TD,Pixel); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
    end
    
    figure(fign+1)
    subplot(2,2,1);
    title('K32_E','FontSize',10)
    nhist(fit_p(1,:)*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
    subplot(2,2,2);
    title('K32_W','FontSize',12)
    nhist((fit_p(2,:))*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
    subplot(2,2,3);
    title('Normalzation','FontSize',12)
    nhist(fit_p(3,:),'text','pdf','color','sequential');xlabel('no unit','FontSize',8); grid on
    subplot(2,2,4);
    title('Background','FontSize',12)
    nhist(((fit_p(4,:))),'text','pdf','color','sequential');xlabel('mcps','FontSize',8); grid on
    switch pub
        case 'ON'
           myname = sprintf('./figures/f2-krkatrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,0); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+1,myname);
    end
    
    figure(fign+2)
    ndhist(fit_p(1,:)*1e3,((fit_p(2,:)))*1e3);
    colorbar
    ylabel('E (meV)','FontSize',10);
    xlabel('W (meV)','FontSize',10);
    R1 = corrcoef(fit_p(1,:),((fit_p(2,:))));
    switch pub
        case 'ON'
           myname = sprintf('./figures/f3-krkatrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,0); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+2,myname);
    end
    
    % Display
    % A.DisplayKrInfo();
                
end
