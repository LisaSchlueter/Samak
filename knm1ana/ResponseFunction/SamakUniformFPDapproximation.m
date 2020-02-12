% Main: Compute RF run-wise and difference to stacked/mean for several qU
SamakCheckUniformFPDapprox('MyqUIndex',2);
SamakCheckUniformFPDapprox('MyqUIndex',10);
SamakCheckUniformFPDapprox('MyqUIndex',15);
SamakCheckUniformFPDapprox('MyqUIndex',20);
SamakCheckUniformFPDapprox('MyqUIndex',25);
SamakCheckUniformFPDapprox('MyqUIndex',30);
SamakCheckUniformFPDapprox('MyqUIndex',35);

function SamakCheckUniformFPDapprox(varargin)
%
% Load All Samak RF for Uniform FPD
% Superimpose RFs
% Compare to Stacked RF
%
% Th. Lasserre
% Last Modified, June 5 2019

p = inputParser;
p.addParameter('MyqUIndex',2);
p.parse(varargin{:});
MyqUIndex    = p.Results.MyqUIndex;

% Load MultiRunAnalysis Object for KNM1
if ~exist('Real','var')
options = {...
    'RunList','KNM1',...
    'exclDataStart',2,...
    'fixPar','1 5 6 7 8 9 10 11',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',1,...
    'NonPoissonScaleFactor',1,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]};
Real=MultiRunAnalysis('DataType','Real',options{:});
end

% Working Retarding Potential
MyqU      = Real.RunData.qU(MyqUIndex);
fprintf('Build plots/KNM1_RunWiseRF_Mean_Stacked_qU%.0f',MyqU);

% Load Run-Wise Response Functions
filelist = [];
if isempty(filelist)
tmpdir = dir([getenv('SamakPath'),sprintf('/inputs/ResponseFunction/')]);
filelist = arrayfun(@(x) x.name,tmpdir,'UniformOutput',0);
filelist(~contains(filelist,'.mat'))=[];
filelist(~contains(filelist,'samakRF_Uniform'))=[];
filelist(~contains(filelist,'Bm4.23T_Bs2.52T_Ba6.3'))=[];
filelist(~contains(filelist,'_Bin100_KatrinD2'))=[];
filelist(~contains(filelist,'Temin18371'))=[];
filelist(~contains(filelist,'Temax1862'))=[];
%filelist =  extractBefore(filelist,'.mat');
end

%% Recover Data & Compute Mean Te & RF
meanRF = zeros(2500,1);
meanTe = zeros(2500,1);
meanqU = 0;
for counter=1:1:numel(filelist)
    rf(counter) = importdata( [getenv('SamakPath') '/inputs/ResponseFunction/' char(filelist(counter))]);
    meanRF = meanRF + rf(counter).RF(:,MyqUIndex);
    meanTe = meanTe + rf(counter).Te;
    meanqU = meanqU + rf(counter).qU(MyqUIndex);
    %fprintf('counter = %d\n',counter);
end
meanRF  = meanRF ./numel(filelist);
meanTe  = meanTe ./numel(filelist);
meanqU  = meanqU ./numel(filelist);

%% Plot All RF + Mean + Stacked
figure('Units', 'pixels', ...
    'Position', [0 0 1200 800]);

s1 = subplot(3,1,[1 2]);
for counter=1:1:numel(filelist)
    r=plot(rf(counter).Te-rf(counter).qU(MyqUIndex),rf(counter).RF(:,MyqUIndex),'LineWidth',5,'Color',rgb('DarkGray'));
    hold on
end
m=plot(meanTe-meanqU,meanRF,'LineWidth',1,'Color',rgb('Red'),'LineStyle','--');
s=plot(Real.ModelObj.Te-Real.ModelObj.qU(MyqUIndex),Real.ModelObj.RF(:,MyqUIndex),'LineWidth',1,'Color',rgb('SteelBlue'),'LineStyle','-.');
hold off
ylabel('transmission','FontSize',18);
grid on            
leg=legend([r m s],'274 KNM1 runs','mean response function','stacked response function','Location','southeast');
leg.Color = 'none'; legend boxoff; leg.FontSize = 20;
title(sprintf('KNM1 Run-Wise Response Functions - qU=%.1f eV',MyqU));
PrettyFigureFormat
set(gca,'FontSize',20);
set(gca,'xscale','log');

% RF - mean(RF)
s2 = subplot(3,1,3);
for counter=1:1:numel(filelist)
    r=plot(rf(counter).Te-rf(counter).qU(MyqUIndex),rf(counter).RF(:,MyqUIndex)-Real.ModelObj.RF(:,MyqUIndex),'LineWidth',15,'Color',rgb('DarkGray'));
    hold on
end
m=plot(meanTe-meanqU,meanRF-Real.ModelObj.RF(:,MyqUIndex),'LineWidth',1,'Color',rgb('Red'),'LineStyle','--');
s=plot(Real.ModelObj.Te-Real.ModelObj.qU(MyqUIndex),Real.ModelObj.RF(:,MyqUIndex)-Real.ModelObj.RF(:,MyqUIndex),'LineWidth',1,'Color',rgb('SteelBlue'),'LineStyle','-.');
hold off
xlabel('surplus energy (eV)','FontSize',18);
ylabel('difference to stacked','FontSize',18);
PrettyFigureFormat
set(gca,'FontSize',20);
linkaxes([s1,s2],'x');
xlim([1e-1 max(max(Real.ModelObj.Te-Real.ModelObj.qU(MyqUIndex)))])
set(gca,'xscale','log');
export_fig(gcf,sprintf('plots/KNM1_RunWiseRF_Mean_Stacked_qU%.0f',MyqU),'-m3');

end
