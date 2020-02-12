    p.addParameter('displayHist','ON',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('displayFit','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('displayAll','ON',@(x)ismember(x,{'ON','OFF'}));

index = obj.ModelObj.Q_i-obj.ModelObj.qU;
            
            % Time
            T     = obj.SingleRunData.TimeSec;
            TC    = cumsum(T); % comulated time
            switch obj.chi2
                case 'chi2Stat'
                    % Endpoint
                    E0    = obj.SingleRun_FitResults.chi2Stat.E0;
                    E0Err = obj.SingleRun_FitResults.chi2Stat.E0Err;
                    % Background
                    B    = obj.SingleRun_FitResults.chi2Stat.B;
                    BErr = obj.SingleRun_FitResults.chi2Stat.BErr;
                    % Normalizationn
                    N      = obj.SingleRun_FitResults.chi2Stat.N;%.*cell2mat(cellfun(@(x) x.NormFactorTBDDS,obj.SingleRunObj,'UniformOutput',0));
                    NErr   = obj.SingleRun_FitResults.chi2Stat.NErr;%.*cell2mat(cellfun(@(x) x.NormFactorTBDDS,obj.SingleRunObj,'UniformOutput',0));
                    % p-value
                    pValue = obj.SingleRun_FitResults.chi2Stat.pValue;
                    Chi2   = obj.SingleRun_FitResults.chi2Stat.chi2min;
                    dof    = obj.SingleRun_FitResults.chi2Stat.dof;
                case 'chi2CM'
                    % Endpoint
                    E0    = obj.SingleRun_FitResults.chi2CMall.E0;
                    E0Err = obj.SingleRun_FitResults.chi2CMall.E0Err;
                    % Background
                    B    = obj.SingleRun_FitResults.chi2CMall.B;
                    BErr = obj.SingleRun_FitResults.chi2CMall.BErr;
                    % Normalizationn

                    N      = obj.SingleRun_FitResults.chi2CMall.N;%.*cell2mat(cellfun(@(x) x.NormFactorTBDDS,obj.SingleRunObj,'UniformOutput',0));
                    NErr   = obj.SingleRun_FitResults.chi2CMall.NErr;%.*cell2mat(cellfun(@(x) x.NormFactorTBDDS,obj.SingleRunObj,'UniformOutput',0));

                    % p-value
                    pValue = obj.SingleRun_FitResults.chi2CMall.pValue;
                    Chi2   = obj.SingleRun_FitResults.chi2CMall.chi2min;
                    dof    = obj.SingleRun_FitResults.chi2CMall.dof;
            end
            % RhoD
            rhoD    = obj.SingleRunData.WGTS_CD_MolPerCm2;
            rhoDErr = obj.SingleRunData.WGTS_CD_MolPerCm2*0;
            % TT
            tt      = obj.SingleRunData.WGTS_MolFrac_TT;
            ttErr   = obj.SingleRunData.WGTS_MolFrac_TT*0;
            % DT
            dt      = obj.SingleRunData.WGTS_MolFrac_DT;
            dtErr   = obj.SingleRunData.WGTS_MolFrac_DT*0;
            % HT
            ht      = obj.SingleRunData.WGTS_MolFrac_HT;
            htErr   = obj.SingleRunData.WGTS_MolFrac_HT*0;
            % Rate
            r1qU    = obj.SingleRunData.TBDIS(1,:)./sum(obj.SingleRunData.qUfrac(1,:,:),3)/numel(obj.PixList);
            r1qUErr = sqrt(r1qU)./sum(obj.SingleRunData.qUfrac(1,:,:),3)/numel(obj.PixList);
                        
            % Labels
            r='';
            for i=1:numel(obj.RunList)
                r    = [r num2str(obj.RunList(i))];
                r    = string(r);
            end
            
            FitResults = struct('E0',E0,'E0Err',E0Err,'N',N,'NErr',NErr,'B',B,'BErr',BErr,...
                'pValue',pValue,'chi2min',Chi2,'dof',dof);
            x   = linspace(1,numel(obj.RunList),numel(obj.RunList));
           
            

       
             
            %% Plot
            switch displayHist
                case 'ON'
                    switch obj.DataType
                        case 'Real'
                            myMainTitle=[sprintf('KATRIN Run-wise Fit - %d Runs (%s) - %.1f',...
                                numel(obj.RunList),obj.chi2,round(index(obj.exclDataStart),1)),'eV below E_0 - ','Run-wise Data'];
                            DataTypeLabel = '';
                        case  'Twin'
                            myMainTitle=[sprintf('KATRIN Run-wise Fit - %d Runs (%s) - %.1f',...
                                numel(obj.RunList),obj.chi2,round(index(obj.exclDataStart),1)),'eV below E_0',' MC Sim'];
                            DataTypeLabel = 'Twin_';
                        case 'FitriumTwin'
                            myMainTitle=[sprintf('KATRIN Run-wise Fit - %d Runs (%s) - %.1f',...
                                numel(obj.RunList),obj.chi2,round(index(obj.exclDataStart),1)),'eV below E_0',' MC Sim'];
                            DataTypeLabel = 'FitriumTwin_';
                        case 'KafitTwin'
                            myMainTitle=[sprintf('KATRIN Run-wise Fit - %d Runs (%s) - %.1f',...
                                numel(obj.RunList),obj.chi2,round(index(obj.exclDataStart),1)),'eV below E_0',' MC Sim'];
                            DataTypeLabel = 'KafitTwin_';
                    end
                    %% Fit Parameters
                    maintitle=myMainTitle;
                    savedir = [getenv('SamakPath'), obj.DataSet,'_FitAllRuns/','plots/'];
                    if ~exist(savedir,'dir')
                        system(['mkdir ',savedir]);
                        system(['mkdir ',[savedir, '/plots/']]);
                    end
                    savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-1.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));

                    fig1 = figure('Name','MRA Run-wise Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
                    a=annotation('textbox', [0 0.91 1 0.1], ...
                        'String', maintitle, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center');
                    a.FontSize=20;a.FontWeight='bold';
                    subplot(2,1,1)
                    bar(x,E0-mean(E0),'facecolor',rgb('IndianRed'));
                    hold on
%                    errorb(x,E0-mean(E0),E0Err,'s','MarkerSize',8,'FaceColor',rgb('DarkGray'),'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',0.5)
                    errorb(x,E0-mean(E0),E0Err,'LineWidth',1,'Color',rgb('DarkSlateGray'));
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel(sprintf('E_0-%.1f (eV)',-mean(E0)));
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp1 = sprintf('<E0>=%.2f eV \\pm %.2f eV (std)',mean(E0),std(E0));title(sp1)
                    grid on
                    %PrettyFigureFormat
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        MyMarkerSize=24;
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        MyMarkerSize=18;
                        set(gca,'FontSize',10);
                        xtickangle(90);
                        MyMarkerSize=12;
                    else
                        xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    subplot(2,1,2)
                    bar(x,B-mean(B),'facecolor',rgb('IndianRed'));
                    hold on
                    errorb(x,B-mean(B),BErr,'LineWidth',1,'Color',rgb('DarkSlateGray'));
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('Background (cps)')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp2 = sprintf('<B>=%.2f cps \\pm %.2f cps (std)',mean(B),std(B));title(sp2)
                    grid on
                    %PrettyFigureFormat
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    %publish_figurePDF(gcf,savefile);
                     export_fig(gcf,[savedir,savefile],'-q101','-m3');
                     
                    maintitle=myMainTitle;

                    savedir = [getenv('SamakPath'),obj.DataSet,'_FitAllRuns/','plots/'];
                  %  savedir = './plots/';

                    if ~exist(savedir,'dir')
                        system(['mkdir ',savedir]);
                    end
                    

                    savefile= sprintf('MRA_RunWiseFitResults%d_%s_%.0feVbelowE0-2.png',...
                                     numel(obj.RunList),obj.chi2,index(obj.exclDataStart));
                                 

                    fig1 = figure('Name','MRA Run-wise Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                    a=annotation('textbox', [0 0.91 1 0.1], ...
                        'String', maintitle, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center');
                    a.FontSize=24;a.FontWeight='bold';
                    subplot(2,1,1)
                    bar(x,N-mean(N),'facecolor',rgb('IndianRed'));
                    hold on
                    %errorbar(x,N,NErr,'s','MarkerSize',MyMarkerSize,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',0.5)
                    errorb(x,N-mean(N),NErr,'LineWidth',1,'Color',rgb('DarkSlateGray'));

                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('N-<N>')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp3 = sprintf('<N-1>=%.2f \\pm %.2f  (std)',mean(N),std(N));title(sp3)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                        xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    %PrettyFigureFormat
                    subplot(2,1,2)
                    bar(x,pValue,'facecolor',rgb('IndianRed'))
                    hold on
                    line([min(x)-0.5,max(x)+0.5],[0.05 0.05],'Color','Red','LineWidth',1,'LineStyle','--')
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('p-value')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp3 = sprintf('<p-value>=%.2f \\pm %.2f (std)',mean(pValue),std(pValue));title(sp3)
                    set(gca,'yscale','log');
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end

                    %publish_figurePDF(gcf,[savedir,savefile]);
                    export_fig(gcf,[savedir,savefile],'-q101','-m3');
                    
                    %% Slow Control Parameters - Colmun Density
                   maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters - Column Density - %s',obj.DataType);
                      
                    savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-3.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));

                    fig1 = figure('Name','MRA Run-wise Column Density','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                    a=annotation('textbox', [0 0.91 1 0.1], ...
                        'String', maintitle, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center');
                    a.FontSize=24;a.FontWeight='bold';
                    subplot(2,1,1)
                    errorbar(x,rhoD./mean(rhoD),rhoDErr./mean(rhoD),'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                    xticks(x)
                    xticklabels(r)
                    ylabel('\rho d / <\rho d>')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp1 = ['<\rho d>',sprintf('=%g mol/cm^2 \\pm %g mol/cm^2 (std)',mean(rhoD),std(rhoD))];title(sp1)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    subplot(2,1,2)
                    errorbar(x,r1qU,r1qUErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                    xticks(x)
                    xticklabels(r)
                    ylabel(sprintf('qU_{%0.1f} Rate / <qU_{%0.1f} Rate>',obj.ModelObj.qU(1)-obj.ModelObj.Q_i,obj.ModelObj.qU(1)-obj.ModelObj.Q_i));
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp1 = ['<qU Rate> ',sprintf('=%g cps \\pm %g  (std) at qU=%0.1f eV below E_0',...
                        mean(r1qU),std(r1qU),obj.ModelObj.qU(1)-obj.ModelObj.Q_i)];title(sp1)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end

                    %publish_figurePDF(gcf,savefile);
export_fig(gcf,[savedir,savefile],'-q101','-m3');

                    %% Slow Control Parameters - Isotopologue
                    maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters - Isotopologues - %s',obj.DataType); 
                    savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-4.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));

                    fig1 = figure('Name','MRA Run-wise Isotopologue Fractions','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                    a=annotation('textbox', [0 0.91 1 0.1], ...
                        'String', maintitle, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center');
                    a.FontSize=24;a.FontWeight='bold';
                    subplot(3,1,1)
                    errorbar(x,tt./mean(tt),dtErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('[TT]/<[TT]>')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp2 = sprintf('<TT>=%g  \\pm %.g (std)',mean(tt),std(tt));title(sp2)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                        xticks(x(1:10:end)); xticklabels(r(1:10:end));
                    end
                    %publish_figurePDF(gcf,savefile);
export_fig(gcf,[savedir,savefile],'-q101','-m3');
subplot(3,1,2)
                    errorbar(x,ht./mean(ht),htErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('[HT]/<[HT]>')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp2 = sprintf('<DT>=%g  \\pm %.g (std)',mean(ht),std(ht));title(sp2)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                        xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    subplot(3,1,3)
                    errorbar(x,dt./mean(dt),dtErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('[DT]/<[DT]>')
                    xlabel('run');
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp2 = sprintf('<DT>=%g  \\pm %.g (std)',mean(dt),std(dt));title(sp2)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    %publish_figurePDF(gcf,savefile);
                    export_fig(gcf,[savedir,savefile],'-q101','-m3');
 
                    %% Data Taking Time Parameters
                    maintitle=sprintf('KATRIN  First Tritium - Slow Control Parameters - Data TakingTime Summary - %s',obj.DataType);
                    savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-5.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));
                    %x   = linspace(1,numel(obj.RunList),numel(obj.RunList));
                    fig1 = figure('Name','MRA Data Taking Time (Runs)','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
                    a=annotation('textbox', [0 0.91 1 0.1], ...
                        'String', maintitle, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center');
                    a.FontSize=24;a.FontWeight='bold';
                    subplot(2,1,1)
                    bar(x,T,'facecolor',rgb('IndianRed'));
                    xticks(x)
                    xticklabels(r)
                    ylabel('Time (sec)')
                    xlabel('run'); 
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp1 = sprintf('<Run Time>=%g seconds \\pm %g eV (std)',mean(T),std(T));title(sp1)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                        xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    subplot(2,1,2)
                    bar(x,TC./86400,'facecolor',rgb('IndianRed'));
                    hold off
                    xticks(x)
                    xticklabels(r)
                    ylabel('Cumulative Time (day)')
                    xlabel('run'); 
                    xlim([0.5 numel(obj.RunList)+0.5])
                    sp2 = sprintf('Total Time %g days',TC(end)./86400);title(sp2)
                    grid on
                    if numel(obj.RunList)<60
                        set(gca,'FontSize',12);
                        xtickangle(45);
                    elseif numel(obj.RunList)<180
                        set(gca,'FontSize',10);
                        xtickangle(90);
                    else
                       xticks(x(1:5:end)); xticklabels(r(1:5:end));
                        xtickangle(45);
                        MyMarkerSize=8;
                    end
                    %publish_figurePDF(gcf,savefile);
export_fig(gcf,[savedir,savefile],'-q101','-m3');                    
                    switch displayAll
                        case 'ON'
                            %% Fit Parameters Verus rhoD
                            maintitle=myMainTitle;
                            savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-6.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));
                            fig1 = figure('Name','MRA Run-wise Fits Versus Column Density','NumberTitle','off','rend','painters','pos',[10 10 1400 900]);
                            a=annotation('textbox', [0 0.91 1 0.1], ...
                                'String', maintitle, ...
                                'EdgeColor', 'none', ...
                                'HorizontalAlignment', 'center');
                            a.FontSize=24;a.FontWeight='bold';
                            subplot(2,2,1)
                            errorbar(rhoD,E0,E0Err,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel(sprintf('E_0-%.1f (eV)',obj.ModelObj.Q_i));
                            xlabel('\rho d (mol/cm^2)');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,2)
                            errorbar(rhoD,B,BErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('B (cps)')
                            xlabel('\rho d (mol/cm^2)');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,3)
                            errorbar(rhoD,N,NErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('N')
                            xlabel('\rho d (mol/cm^2)');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,4)
                            errorbar(rhoD,pValue,pValue.*0,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('p-value')
                            xlabel('\rho d (mol/cm^2)');xtickangle(45);
                            set(gca,'yscale','log');
                            PrettyFigureFormat
                            %publish_figurePDF(gcf,savefile);
                            export_fig(gcf,[savedir,savefile],'-q101','-m3');
                            
                            %% Fit Parameters Verus DT
                            maintitle=myMainTitle;
                            savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-7.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));
                            fig1 = figure('Name','MRA Run-wise Fits Versus [DT]','NumberTitle','off','rend','painters','pos',[10 10 1400 900]);
                            a=annotation('textbox', [0 0.91 1 0.1], ...
                                'String', maintitle, ...
                                'EdgeColor', 'none', ...
                                'HorizontalAlignment', 'center');
                            a.FontSize=24;a.FontWeight='bold';
                            subplot(2,2,1)
                            errorbar(dt,E0,E0Err,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel(sprintf('E_0-%.1f (eV)',obj.ModelObj.Q_i));
                            xlabel('[DT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,2)
                            errorbar(dt,B,BErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('B (cps)')
                            xlabel('[DT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,3)
                            errorbar(dt,N,NErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('N')
                            xlabel('[DT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,4)
                            errorbar(dt,pValue,pValue.*0,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('p-value')
                            xlabel('[DT]');xtickangle(45);
                            set(gca,'yscale','log');
                            PrettyFigureFormat
                            %publish_figurePDF(gcf,savefile);
                            export_fig(gcf,[savedir,savefile],'-q101','-m3');
                            
                            %% Fit Parameters Verus TT
                            maintitle=myMainTitle;
                            savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-8.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));
                            fig1 = figure('Name','MRA Run-wise Fits Versus [TT]','NumberTitle','off','rend','painters','pos',[10 10 1400 900]);
                            a=annotation('textbox', [0 0.91 1 0.1], ...
                                'String', maintitle, ...
                                'EdgeColor', 'none', ...
                                'HorizontalAlignment', 'center');
                            a.FontSize=24;a.FontWeight='bold';
                            subplot(2,2,1)
                            errorbar(tt,E0,E0Err,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel(sprintf('E_0-%.1f (eV)',obj.ModelObj.Q_i));
                            xlabel('[TT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,2)
                            errorbar(tt,B,BErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('B (cps)')
                            xlabel('[TT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,3)
                            errorbar(tt,N,NErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('N')
                            xlabel('[TT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,4)
                            errorbar(tt,pValue,pValue.*0,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('p-value')
                            xlabel('[TT]');xtickangle(45);
                            set(gca,'yscale','log');
                            PrettyFigureFormat
                            %publish_figurePDF(gcf,savefile);
                            export_fig(gcf,[savedir,savefile],'-q101','-m3');
                            
                            %% Fit Parameters Verus HT
                            maintitle=myMainTitle;
                            savefile=sprintf('MRA_RunWiseFitResults%d_%s%s_%.0feVbelowE0-9.png',numel(obj.RunList),DataTypeLabel,obj.chi2,index(obj.exclDataStart));
                            fig1 = figure('Name','MRA Run-wise Fits Versus [HT]','NumberTitle','off','rend','painters','pos',[10 10 1400 900]);
                            a=annotation('textbox', [0 0.91 1 0.1], ...
                                'String', maintitle, ...
                                'EdgeColor', 'none', ...
                                'HorizontalAlignment', 'center');
                            a.FontSize=24;a.FontWeight='bold';
                            subplot(2,2,1)
                            errorbar(ht,E0,E0Err,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel(sprintf('E_0-%.1f (eV)',obj.ModelObj.Q_i));
                            xlabel('[HT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,2)
                            errorbar(ht,B,BErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('B (cps)')
                            xlabel('[HT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,3)
                            errorbar(ht,N,NErr,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('N')
                            xlabel('[HT]');xtickangle(45);
                            PrettyFigureFormat
                            subplot(2,2,4)
                            errorbar(ht,pValue,pValue.*0,'ks','MarkerSize',MyMarkerSize,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1);
                            ylabel('p-value')
                            xlabel('[HT]');xtickangle(45);
                            set(gca,'yscale','log');
                            PrettyFigureFormat
                            %publish_figurePDF(gcf,savefile);
                            export_fig(gcf,[savedir,savefile],'-q101','-m3');
                            
                            %% Build pdf file with all fits
                            %                     if ~isunix %does not work for Lisa :(
                            %                         PATH = getenv('PATH');
                            %                         setenv('PATH', [PATH ':/usr/local/bin/']);
                            %                         cd plots/runwise
                            %                         command = 'gs -sDEVICE=pdfwrite -sOutputFile="run.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f KATRIN_FT_AllPixels_Samak_*.pdf -c quit';
                            %                         unix(command);
                            %                         unix('rm KATRIN_FT_AllPixels_Samak_*.pdf');
                            %                         mycommand1 = sprintf('mv run.pdf ../MRA-RunWiseFits-%s-exclDataStart%g-%s-samak.pdf',chi2,exclDataStart,datestr(now,'dd-mmm-yyyy'));
                            %                         unix(mycommand1);
                            %                         mycommand2 = sprintf('open ../MRA-RunWiseFits-%s-exclDataStart%g-%s-samak.pdf',chi2,exclDataStart,datestr(now,'dd-mmm-yyyy'));
                            %                         unix(mycommand2);
                            %                         cd ../..
                            %                         unix('rm -rf plots/tmp');
                            %                     end
                            
                            %save fit results in txt file
                            cd plots
                            filename = sprintf('MRA_FitResults_%s_%.0fbelowE0_%uRuns_%u_%u.txt',obj.chi2,index(obj.exclDataStart),numel(obj.RunList),obj.RunList(1),obj.RunList(end));
                            unix(['touch ',filename]);
                            Endpoint = [FitResults.E0+obj.ModelObj.Q_i; FitResults.E0Err];
                            Background = [FitResults.B; FitResults.BErr];
                            Normalization = [FitResults.N+1; FitResults.NErr];
                            
                            fileID = fopen(filename,'w');
                            mytitle = sprintf('Fit to %u Runs: %u-%u \n',numel(obj.RunList),obj.RunList(1),obj.RunList(end));
                            fprintf(fileID,mytitle);
                            fprintf(fileID,'mean E0 =%.3f \n',mean(FitResults.E0+obj.ModelObj.Q_i));
                            fprintf(fileID,'mean E0err = %.3f \n', mean(FitResults.E0Err));
                            fprintf(fileID,'err of mean (E0) = %.3f \n', sqrt(var(FitResults.E0)));
                            fprintf(fileID,'chi2/dof \n');
                            fprintf(fileID,'%.2f / %.0f \n',[FitResults.chi2min; FitResults.dof]);
                            fprintf(fileID,'p-value \n');
                            fprintf(fileID,'%.3f \n',FitResults.pValue);
                            fprintf(fileID,'Endpoint (eV) \n');
                            fprintf(fileID,'%.3f +/- %.3f \n',Endpoint);
                            fprintf(fileID,'Backgroud (mps) \n');
                            fprintf(fileID,'%.3f +/- %.3f \n',Background);
                            fprintf(fileID,'Normalization \n');
                            fprintf(fileID,'%.4f +/- %.4f \n',Normalization);
                            fclose(fileID);
                            cd ..
                    end
            end
            