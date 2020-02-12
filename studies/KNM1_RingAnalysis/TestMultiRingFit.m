mypath = [getenv('SamakPath'),'/studies/KNM1_RingAnalysis'];
range = 30;
if range==30
    excludedataStart = 17;
elseif range==90
    excludedataStart = 1;
end
RingList = 2:8;
qUOffsetMode = 'interpRF';%'reinitRF';%'interpolTBDIS';%
savename = sprintf('TestMultiRingFit_range%.0f_ringlist%s_%s.mat',range,strrep(join(string(RingList)),' ',''),qUOffsetMode);

if exist([mypath,'/results/',savename],'file')
    load([mypath,'/results/',savename]);
else
    %%
    M =  MultiRunAnalysis('RunList','KNM1_m149mvRW','AnaFlag','Ring','RingList',RingList,'PullFlag',5);
    FSDStart = 2.*M.ModelObj.nPixels+3;
    FSDStop  = 2.*M.ModelObj.nPixels+3+5;
    M.fixPar = ['1 ',char(strjoin(string(FSDStart:FSDStop)))];
   tic; M.Fit; toc;
    M.PlotFit;
    M.exclDataStart = excludedataStart;
    %%
    A =  MultiRunAnalysis('RunList','KNM1_m149mvRW','AnaFlag','StackPixel','RingList',RingList,'PullFlag',5);
    A.fixPar = '1 5 6 7 8 9 10 11';
    A.exclDataStart = excludedataStart;
    R = RingAnalysis('RunAnaObj',A,'RingList',RingList);
    R.FitRings;
    save([mypath,'/results/',savename],'M','A','R');
end
 %% compare
  FSDStart = 2.*M.ModelObj.nPixels+3;
    FSDStop  = 2.*M.ModelObj.nPixels+3+5;
 fig1 = figure('Renderer','openGL');
 set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
 
 eMultiRing = errorbar(RingList,M.FitResult.par(FSDStop+1:end),M.FitResult.err(FSDStop+1:end),'o-','LineWidth',3,'MarkerSize',8,'Color',rgb('SeaGreen'));
 hold on;
 SingleRingmeanE0 = mean(R.FitResult.par(:,2));
 eSingleRing = errorbar(RingList,R.FitResult.par(:,2)-SingleRingmeanE0,R.FitResult.err(:,2),'s-','LineWidth',3,'MarkerSize',8,'Color',rgb('DarkOrange'));
 
 PrettyFigureFormat;
 xlabel('ring');
 ylabel('qU offset (multiring) OR relative endpoint (single ring) (eV) ');
 set(gca,'FontSize',18);
 
leg = legend([eMultiRing,eSingleRing],'multi-ring fit','single ring fits'); legend boxoff
title(sprintf('%s: stacked runs %.0f - %.0f',strrep(M.ModelObj.TD,'_m',' -'),M.RunList(1),M.RunList(end)));
print([mypath,'/plots/',strrep(savename,'.mat','.png')],'-dpng','-r450')
%%
fprintf(' multi ring fit: endpoint %.2f +/- %.2f eV \n',M.FitResult.par(2)+M.ModelObj.Q_i,M.FitResult.err(2));
fprintf(' ringwise fit: endpoints  %.2f +/- %.2f eV \n',[R.FitResult.par(:,2)+M.ModelObj.Q_i,R.FitResult.err(:,2)]');


