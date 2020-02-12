RunList = 'KNM1_m149mvRW';
save_file = [getenv('SamakPath'),'knm1ana/RateAt200eV/results/',sprintf('RateAt200eV_%s.mat',RunList)];
if exist(save_file,'file')
    load(save_file)
else
%% Get the real data points
D =  MultiRunAnalysis('RunList',RunList,'DataType','Real');
%% Compute Rate by changing column density, fix tritium puritiy to mean value
WGTS_CD = linspace(min(D.SingleRunData.WGTS_CD_MolPerCm2),max(D.SingleRunData.WGTS_CD_MolPerCm2),20);
TBDIS   = zeros(D.ModelObj.nqU,numel(WGTS_CD));
RF      = zeros(D.ModelObj.nTe, D.ModelObj.nqU,numel(WGTS_CD));

D.ModelObj.TimeSec = mean(D.SingleRunData.TimeSec); % take average time
for i=1:numel(WGTS_CD)
    progressbar(i/numel(WGTS_CD));
    
    % loop over column densities
    D.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD(i);
    
    % Compute new RF and adjust normalization
    D.ModelObj.AdjustRF; 
    D.ModelObj.ComputeTBDDS; D.ModelObj.ComputeTBDIS;
    
    %save
    TBDIS(:,i) = D.ModelObj.TBDIS;
    RF(:,:,i) = D.ModelObj.RF; 
end

SimActivity = WGTS_CD.*(D.RunData.WGTS_MolFrac_TT + ...
    0.5*D.RunData.WGTS_MolFrac_DT + 0.5*D.RunData.WGTS_MolFrac_HT);
SimRate  = TBDIS(1,:)./(D.RunData.qUfrac(1).*D.RunData.TimeSec);
SimErr   = sqrt(TBDIS(1,:))./(D.RunData.qUfrac(1).*D.RunData.TimeSec);
SimMean  = mean(SimRate);
%% Compute rate with twins: The same tritium purity and column density values as in real data
M =  MultiRunAnalysis('RunList',RunList,'DataType','Twin');
%%
DataActivity = D.SingleRunData.WGTS_CD_MolPerCm2.*...
     (D.SingleRunData.WGTS_MolFrac_TT + 0.5*D.SingleRunData.WGTS_MolFrac_DT + 0.5*D.SingleRunData.WGTS_MolFrac_HT);
DataRate = D.SingleRunData.TBDIS(1,:)./(D.SingleRunData.qUfrac(1,:).*D.SingleRunData.TimeSec);
DataErr  = sqrt(D.SingleRunData.TBDIS(1,:))./(D.SingleRunData.qUfrac(1,:).*D.SingleRunData.TimeSec);
DataMean =  mean(DataRate);

SimActivity_Twin = M.SingleRunData.WGTS_CD_MolPerCm2.*...
     (M.SingleRunData.WGTS_MolFrac_TT + 0.5*M.SingleRunData.WGTS_MolFrac_DT + 0.5*M.SingleRunData.WGTS_MolFrac_HT);
SimRate_Twin  = M.SingleRunData.TBDIS(1,:)./(M.SingleRunData.qUfrac(1,:).*M.SingleRunData.TimeSec);
SimErr_Twin   = sqrt(M.SingleRunData.TBDIS(1,:))./(M.SingleRunData.qUfrac(1,:).*M.SingleRunData.TimeSec);
SimMean_Twin  = mean(SimRate_Twin);
Runs = D.RunList;
if ~exist([getenv('SamakPath'),'knm1ana/RateAt200eV/results/'],'dir')
    system(['mkdir ',getenv('SamakPath'),'knm1ana/RateAt200eV/results/'])
end
save(save_file,'DataActivity','DataRate','DataMean','DataErr','SimActivity','SimRate','SimMean','SimErr',...
    'SimActivity_Twin','SimRate_Twin','SimMean_Twin','SimErr_Twin','Runs');
end
%% plot
close
fig5 = figure('Renderer','openGL');
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);

pdata = errorbar(DataActivity,DataRate./DataMean,DataErr./DataMean,...
    'o','MarkerSize',8,'Color',rgb('RoyalBlue'),'MarkerFaceColor',rgb('CadetBlue'));
hold on;
psim_twin  = errorbar(SimActivity_Twin,SimRate_Twin./SimMean_Twin, SimErr_Twin./SimMean_Twin,'s',...
       'MarkerSize',8,'Color',rgb('DarkOrange'),'MarkerFaceColor',rgb('Orange'));
psim  = errorbar(SimActivity,SimRate./SimMean, SimErr./SimMean,'d',...
       'MarkerSize',8,'Color',rgb('DarkSlateGray'),'MarkerFaceColor',rgb('SlateGray'));   
hold off;
leg = legend([pdata, psim_twin,psim],'Data','MC Twin','MC const. tritium purity','Location','northwest'); legend boxoff
xlabel('~ tritium activity (mol/cm^2)');
ylabel('rel. rate @ -200V');
PrettyFigureFormat;
set(gca,'FontSize',18);
grid on;
print(fig5,'./plots/m200eVRate','-dpng','-r450')