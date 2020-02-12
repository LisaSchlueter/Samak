%clear; close all;

% Add path to samak folder
addpath(genpath('../../../Samak2.0'));

RunList = load('RunListFTFullCD.mat');
RunList = RunList.RunListFTFullCD;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];

mpix = 'mpix'; % mpix
ringCutFlag = ''; %ex2
N = 10000;

% Choose type of fit (chi2)
chi2name = 'chi2CM';

% Choose data range to analyze
dataStart = 1; % 9 = -200 eV; 1 = -1600 eV

% Choose fitter
fitter = 'minuit';

DataEffCorr = 'OFF';

load('FAKEFTALLex2');

PixelList = 1;

RunList = 'StackCD100all';

% Choose fixed parameters (optional)

fixpar = num2str([2+2*length(PixelList)+1 2+2*length(PixelList)+2]);
%     


for r = 1:1
%     try
    TBDIS = cell(N,1);
    for ring = 1:11
    A_MC = ref_RunAnalysis('FTALL',mpix,ringCutFlag,'FPD_Segmentation','RING',...
        'nTeBinningFactor',20,'FPD_Pixel',1:148,'FPD_Ring',ring);
    A_MC.ComputeTBDDS(); A_MC.ComputeTBDIS();
    
    MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
    'chi2','chi2CM','RunList',RunList,...
    'fitter',fitter,'exclDataStart',dataStart,...
    'fixPar',fixpar,...
    'DataEffcorr','OFF','ring',ring); %'ROI+PileUp
    
    saveTBDIS = A_MC.TBDIS;
    %qU = A_MC.qU;
    %WGTS_CD_MolPerCm2 = A_MC.WGTS_CD_MolPerCm2;
    
    
%     TBDIS_NO = TBDIS;
    for mc = 1:N
     A_MC.AddStatSystFluctTBDIS('CM',MRA.FitCM);
       % A_MC.AddStatFluctTBDIS();
        TBDIS{mc, ring} = A_MC.TBDIS;
        % Save stacked pixels excluding two outer rings
        %TBDIS_NO{mc} = sum(TBDIS{mc}(:,1:124),2);
        
        A_MC.TBDIS = saveTBDIS;
        disp(mc);
    end
    end
    
%     save(['fake-tritium-data/FAKE','FTALL',mpix,ringCutFlag,'.mat'],...
%         'TBDIS','qU','TimeSec','qUfrac',...
%         'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
%         'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
%         'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun')
    

%     qU = mean(qU(:,1:124),2);
%     TBDISmpix = TBDIS;
%     TBDIS = TBDIS_NO;

    A_MC = ref_RunAnalysis('FTALL',mpix,ringCutFlag,'FPD_Segmentation','MULTIPIXEL',...
        'nTeBinningFactor',5,'FPD_Pixel',1:148,'FPD_Ring',ring);
    qU = A_MC.qU;
    
    save(['fake-tritium-data/FAKE','FTALL','ringStat.mat'],... %num2str(RunList(r))
        'TBDIS','qU','TimeSec','qUfrac',...
        'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
        'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
        'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun')
%     catch
%         didntwork(a) = RunList(r);
%         a = a + 1;
%     end
    
end