clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));


RingList = 1:11;


MC = 10000;

% mnu, e0, bck, norm

titlelist = {'neutrino mass squared','endpoint','background','normalization'};
xnames = {'neutrino mass squared (eV)', 'endpoint bias E_0 - 18573.7 (eV)',...
    'background (cps)','normalization'};
legendnames = {'','eV','cps',''};
savelist = {'mnu','e0','b','n'};

for ring = 1:length(RingList)
    R = load(['data/FT_MC_RINGSTAT',num2str(ring),'MC10000','.mat']);
    Result = squeeze(R.ResultsTableRing(:,ring,:));
    Result(Result(:,1) < -300,:) = [];
    Result(:,2) = Result(:,2) - 18573.7;
    Result(:,4) = Result(:,4);
    Result(:,3) = Result(:,3);
    close all;
    for dis = 1:4
        figdis = figure(dis);
        pos = get(figdis,'position');
        set(figdis,'position',[pos(1:2)/3.5 pos(3:4)*1.5])
        %set(figdis,'Visible','off');
        set(figdis,'color','white');
        %screensize = (get(0, 'Screensize'));
        %set(figbar, 'Position', [screensize(1:2)*0.95,screensize(4),screensize(4)*0.95]);
        if dis < 5
            dis_temp = dis;
        else
            dis_temp = 9;
        end
        
        fittedDist = fitdist(Result(:,dis_temp),'Normal');
        mu_hist = fittedDist.mu;
        
        sigma_hist = fittedDist.sigma;
        hold on
        [~,Nbins,Wbins] = nhist(Result(:,dis_temp));
        binWidth = Wbins(2) - Wbins(1);
        x_N = linspace(min(Result(:,dis_temp)),max(Result(:,dis_temp)),...
            length(Nbins)*50);
        plot(x_N,length(Result(:,dis_temp))*binWidth*normpdf(x_N,mu_hist,sigma_hist),'LineWidth',4);
        hold off
        
        
        titlehandle = title([titlelist{dis},' (MC trials = ',num2str(MC),')']);
        xlabelhandle = xlabel(xnames{dis});
        %titlehandle.FontSize = 7;
        PrettyFigureFormat;
        %sysunc = 0.017;
        sysunc = 0.0;
        if dis == 1
            totalunc = sqrt(sysunc.^2 + sigma_hist.^2);
            legend({sprintf('mean = %0.5g eV^2 \n \\sigma_{stat} = %0.3g eV^2',...
                mu_hist,sigma_hist)},'interpreter','tex');
            mu_mc(ring) = mu_hist; sigma_mc(ring) = sigma_hist;
        else
            legend(sprintf(['mean = %0.5g ',legendnames{dis},'\n \\sigma_{stat} = %0.3g ',legendnames{dis}],...
                mu_hist,sigma_hist));
        end
        savename{dis} = ['plots/FT','RINGSTAT',savelist{dis},num2str(MC),'FTdist.pdf'];
        if ispc; export_fig(savename{dis},'-pdf');end
    end
    %append_pdfs(['plots/All.pdf'],savename{:});
    %delete(savename{:});
end
save('data/SensitivityRINGSTAT');
% figure(5)
% hold on
% errorbar(MC,mu_mc,sigma_mc,'ko',...
%     'MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% %scatter(MC,mu_mc,'ko','MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% plot(MC,zeros(1,length(MC)),'Color','red','LineStyle','--','LineWidth',3);
% hold off

% 'facecolor',rgb('CadetBlue')
%
%     annotation('textbox',[0.6 0.5 0.2 0.2],'String',...
%         {['\sigma = ',num2str(sigma_N),' eV'],['mean = ',num2str(mu_N),' eV']},'FitBoxToText','on');

