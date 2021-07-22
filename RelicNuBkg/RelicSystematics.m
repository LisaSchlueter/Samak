function RelicSystematics(saveplot,DataType) 
    startup;
    % Recompute = 'OFF';
    % B=RelicNuAnalysis('Params','KNM1');
     matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
    % B.SystBreakdown('TwinBias_mnuSq',0,'Plot','OFF','Recompute',Recompute);
    % B.SystBreakdown('TwinBias_mnuSq',1,'Plot','OFF','Recompute',Recompute);
    % B.SystBreakdown('TwinBias_mnuSq',5,'Plot','OFF','Recompute',Recompute);
    % % B.SystBias('TwinBias_mnuSq',0,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
    % % B.SystBias('TwinBias_mnuSq',1,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
    % % B.SystBias('TwinBias_mnuSq',5,'Recompute',Recompute,'Plot','ON','DeltaChi2',1);
    load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq0_range40_RealData.mat']);
    YD=Y;
    %load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq0_range40.mat']);
    %Y0=Y;
    load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq1_range40.mat']);
    Y1=Y;
    %load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq5_range40.mat']);
    %Y5=Y;
    load([matFilePath,'FakeSensitivityBreakdown_TDR_mnuSq0_range40.mat']);
    YF=Y;
    if strcmp(DataType,'Twin2')
        load([matFilePath,'SensitivityBreakdown_KNM2_Prompt_mnuSq1_range40.mat']);
    elseif strcmp(DataType,'real2')
        load([matFilePath,'SensitivityBreakdown_KNM2_Prompt_mnuSq0_range40_RealData.mat']);
    end
    Y2=Y;
    %load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq0.mat']);
    %Bias0=Y;Bias0(3)=Bias(4);Bias0(4)=Bias(1);Bias0(5)=Bias(2);Bias0(6)=Bias(3);Bias0(7)=Bias(6);Bias0(8)=Bias(5);Bias0(9)=Bias(7);
    %load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq1.mat']);
    %Bias1=Y;Bias1(3)=Bias(4);Bias1(4)=Bias(1);Bias1(5)=Bias(2);Bias1(6)=Bias(3);Bias1(7)=Bias(6);Bias1(8)=Bias(5);Bias1(9)=Bias(7);
    %load([matFilePath,'SystEffectOnSensitivity_KNM1_mnuSq5.mat']);
    %Bias5=Y;Bias5(3)=Bias(4);Bias5(4)=Bias(1);Bias5(5)=Bias(2);Bias5(6)=Bias(3);Bias5(7)=Bias(6);Bias5(8)=Bias(5);Bias5(9)=Bias(7);

    % fig1=figure(1);
    % bar(X,[Y0;Y1;Y5]);
    % ax=gca;
    % ax.YScale='log';
    % ylabel('\eta');
    % legend('m_{\nu}^{2}=0 eV^{2}','m_{\nu}^{2}=1 eV^{2}','m_{\nu}^{2}=5 eV^{2}');
    % legend boxoff;
    % PrettyFigureFormat;
    % 
    % fig10=figure(10);
    % Y0rel = Y0./Y0(1);
    % Y1rel = Y1./Y1(1);
    % Y5rel = Y5./Y5(1);
    % bar(X,[Y0rel;Y1rel;Y5rel]);
    % ax=gca;
    % ax.YScale='log';
    % legend('m_{\nu}^{2}=0 eV^{2}','m_{\nu}^{2}=1 eV^{2}','m_{\nu}^{2}=5 eV^{2}');
    % legend boxoff;
    % ylabel('Fraction of total');
    % PrettyFigureFormat;

    fig100=figure('Renderer','painters');
    set(fig100, 'Units', 'normalized', 'Position', [0.001, 0.001,0.5, 0.5]);
    X = categorical({'Total','Statistical','Final-state distribution','Response function','Scan fluctuations','Stacking','Detector efficiency','Theoretical corrections','Bkg slope','Bkg rate'});
    switch DataType
        case 'real'
            X = reordercats(X,{'Response function','Theoretical corrections','Scan fluctuations','Detector efficiency','Stacking','Final-state distribution','Bkg slope','Bkg rate','Statistical','Total'});
        case 'Twin'
            X = reordercats(X,{'Bkg slope','Final-state distribution','Theoretical corrections','Scan fluctuations','Response function','Stacking','Detector efficiency','Bkg rate','Statistical','Total'});
        case 'Twin2'
            X= categorical({'Total','Statistical','Final-state distribution','Response function','Scan fluctuations','Stacking','Detector efficiency','Theoretical corrections','Bkg slope','Bkg rate','Plasma','Penning trap'});
            X= reordercats(X,{'Detector efficiency','Penning trap','Theoretical corrections','Response function','Final-state distribution','Scan fluctuations','Stacking','Plasma','Bkg slope','Bkg rate','Statistical','Total'});
        case 'real2'
            X= categorical({'Total','Statistical','Final-state distribution','Response function','Scan fluctuations','Stacking','Detector efficiency','Theoretical corrections','Bkg slope','Bkg rate','Plasma','Penning trap'});
            X= reordercats(X,{'Bkg slope','Scan fluctuations','Theoretical corrections','Stacking','Detector efficiency','Final-state distribution','Response function','Plasma','Penning trap','Bkg rate','Statistical','Total'});
        case 'Fake'
            X = reordercats(X,{'Bkg rate','Detector efficiency','Scan fluctuations','Stacking','Bkg slope','Theoretical corrections','Response function','Final-state distribution','Statistical','Total'});
    end
    bsingle = cell(numel(X),1);

    PlotColor = {rgb('White'),rgb('Navy'),rgb('GoldenRod'),rgb('PowderBlue'),...
                    rgb('CadetBlue'),rgb('DarkOrange'),rgb('FireBrick'),rgb('DarkSlateGray'),...
                    rgb('YellowGreen'),rgb('Magenta'),...
                    rgb('SeaGreen'),rgb('DodgerBlue'),rgb('DarkGreen'),rgb('Red')};

    hold on;
    for i=1:numel(X)
        switch DataType
            case 'real'
                bsingle{i}  = barh(X(i),YD(i));
            case 'Twin'
                bsingle{i}  = barh(X(i),Y1(i));
            case 'Twin2'
                bsingle{i}  = barh(X(i),Y2(i));
            case 'real2'
                bsingle{i}  = barh(X(i),Y2(i));
            case 'Fake'
                bsingle{i}  = barh(X(i),YF(i));
        end
        btmp= bsingle{i};
        btmp.LineStyle = 'none';
        btmp.FaceColor = PlotColor{i+1}; btmp.LineStyle ='none';
    end
    %B=barh(X,Y1,'green');
    ax=gca;
    ax.XScale='log';
    xlabel('\eta');

    % SysUncertainties = bsingle{1}.YData;
    % a=annotation('textbox', [0.14 0.1 1 0.12], ...
    %     'String', SysUncertainties, ...
    %     'EdgeColor', 'none', ...
    %     'HorizontalAlignment', 'left');
    % a.FontSize=16;a.FontWeight='normal';

    effs = ones(1,numel(bsingle));
    vals = ones(1,numel(bsingle));
    labels = string(zeros(10,1));
    for i=1:numel(bsingle)
        effs(i) = bsingle{i}.YEndPoints;
        vals(i) = bsingle{i}.XEndPoints;
        labels(i)=sprintf('%.2g',bsingle{i}.YData);
        if effs(i)<1e9 && ~strcmp(DataType,'Fake') && ~strcmp(DataType,'Twin2') && ~strcmp(DataType,'real2')
            effs(i)=1e9;
            labels(i)='<1e9';
        end
        if effs(i)<1e6 && strcmp(DataType,'Fake') && i==5
            effs(i)=1e6;
            labels(i)='<1e6';
        end
        if effs(i)<1e8 && (strcmp(DataType,'Twin2') || strcmp(DataType,'real2'))
            effs(i)=1e8;
            labels(i)='<1e8';
        end
    end
    text(effs,vals,labels,'HorizontalAlignment','left',...
       'VerticalAlignment','middle','FontSize',16);
    switch DataType
       case 'real'
            xlim([1e9,1e12]);
       case 'Twin'
           xlim([1e9,1e12]);
       case 'Twin2'
           xlim([1e8,5e11]);
       case 'real2'
            xlim([1e8,5e11]);
       case 'Fake'
           ylim(categorical({'Scan fluctuations','Total'}));
           xlim([1e6 2e10]);
    end
    PrettyFigureFormat;

    % fig2=figure(2);
    % bar(X,[Y0;Bias0]);
    % ylabel('\eta');
    % legend('m_{\nu}^{2}=0 ev^{2}, CM','m_{\nu}^{2}=0 ev^{2}, envelope');
    % legend boxoff;
    % PrettyFigureFormat;
    % 
    % fig3=figure(3);
    % bar(X,[Y1;Bias1]);
    % ylabel('\eta');
    % legend('m_{\nu}^{2}=1 ev^{2}, CM','m_{\nu}^{2}=1 ev^{2}, envelope');
    % legend boxoff;
    % PrettyFigureFormat;
    % 
    % fig4=figure(4);
    % bar(X,[Y5;Bias5]);
    % ylabel('\eta');
    % legend('m_{\nu}^{2}=5 ev^{2}, CM','m_{\nu}^{2}=5 ev^{2}, envelope');
    % legend boxoff;
    % PrettyFigureFormat;

    if strcmp(saveplot,'ON')
        SaveDir = [getenv('SamakPath'),sprintf('RelicNuBkg/Plots/FinalPlots/')];
        MakeDir(SaveDir);
        SaveName = sprintf('SysBD_%s.pdf',DataType);
        export_fig(fig100,[SaveDir,SaveName]);
    end