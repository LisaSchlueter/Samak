M =  MultiRunAnalysis('RunList','KNM1','AnaFlag','StackPixel','fixPar','5 6 7 8 9 10');
TimeSec = 60*60*2.*(250:10:360);

mNuSq90_90eV = zeros(numel(TimeSec),1);
mNuSq90_30eV = zeros(numel(TimeSec),1);

progressbar('calculating sensitivity ....');
for i=1:numel(TimeSec)
    progressbar(i/numel(TimeSec));
    M.ModelObj.TimeSec = TimeSec(i);
    M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
    
    M.RunData.TBDIS  = M.ModelObj.TBDIS;
    M.RunData.TBDISE = M.ModelObj.TBDIS.^0.5;
    
    M.exclDataStart = 1;
    M.Fit;
    mNuSq90_90eV(i) = M.FitResult.err(1)*1.64;
    
    M.exclDataStart = 17;
    M.Fit;
    mNuSq90_30eV(i) = M.FitResult.err(1)*1.64;
    
end

mNu90_90eV = sqrt(mNuSq90_90eV);
mNu90_30eV = sqrt(mNuSq90_30eV);

%save([genpath('SamakPath'),'/knm1ana/knm1sensitivity/results/StatSensitivityKnm1.mat]'],'mNu90_90eV','mNu90_30eV');
%%
range = 30;

fig1 = figure('Renderer','opengl');
set(fig1,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);  

switch range
    case 90
 plot(TimeSec./(2*60*60),mNu90_90eV,...
    '-','LineWidth',4,'Color',rgb('SkyBlue'));
    case 30
 plot(TimeSec./(2*60*60),mNu90_30eV,...
    '-','LineWidth',4,'Color',rgb('FireBrick'));
end

leg = legend([num2str(range),' eV range']); legend boxoff

xlabel('nRuns (a 2 hours)')
ylabel(sprintf('m_\\nu sensitivity 90%% C.L. (eV)'));
grid on;

xlim([250,360]);
PrettyFigureFormat;
set(gca,'FontSize',20);
title('KNM1 settings')

print(fig1,[getenv('SamakPath'),'/knm1ana/knm1sensitivity/plots/StatSensitivity',num2str(range),'eV'],'-dpng','-r450');
