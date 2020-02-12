%% load benchmark
Fm1eV = importdata('differential_spectrum_nuMassSqNeg.txt'); % Fitrium with mNuSq  -1eV
F1eV =  importdata('differential_spectrum_nuMassSqPos.txt');  % Fitrium with mNuSq  +1eV
Te = Fm1eV(:,1);

F1eVN  = F1eV(:,2)./sum(F1eV(:,2));    % Normalized
Fm1eVN = Fm1eV(:,2)./sum(Fm1eV(:,2));  % Normalized
%%
Q_i = 18575; % endpoint

% compute model for mNu -1eV
Tm1eV = ref_TBD_DScomparison('Q_i',Q_i,'mNuSq_i',-1);  
Tm1eV.ComputeTBDDS;
Sm1eV  = Tm1eV.TBDDS;                 % Samak Differential Spectrum
Sm1eVN = Sm1eV./sum(Sm1eV);                             % Normalized

% compute model for mNu +1eV
T1eV = ref_TBD_DScomparison('Q_i',Q_i,'mNuSq_i',1);    
T1eV.ComputeTBDDS; 
S1eV = T1eV.TBDDS;%  interp1(T1eV.Te,,Te,'lin');         % Samak Differential Spectrum
S1eVN = S1eV./sum(S1eV);                                % Normalized

if ~all(T1eV.Te==Te)
    fprintf(2,'Error: Not the same Te! \n')
end
%% compute integral: just integral whole energy range and add background
if ~exist('M','var')
    M = MultiRunAnalysis('RunList','KNM1');
end
Tactivity = 2.0167*1e2;% Tritum Activity 
Time = M.RunData.qUfrac(14).*M.RunData.TimeSec;                 %time spent approx 40eV below endpoint

qu_i =[1:10:numel(Te)-1]; 
qu = Te(qu_i)-Q_i;
nqu = numel(qu_i);

SISm1eV = zeros(nqu,1); SIS1eV = zeros(nqu,1); FISm1eV = zeros(nqu,1); FIS1eV = zeros(nqu,1);
savedir = [getenv('SamakPath'),'knm1ana/knm1_DScomparison/plots/'];

%integrate and add tritium activity
for i=1:nqu
SISm1eV(i) = simpsons(Tm1eV.Te(qu_i(i):end),Tactivity.*Sm1eVN(qu_i(i):end));
SIS1eV(i) =  simpsons(T1eV.Te(qu_i(i):end),Tactivity.*S1eVN(qu_i(i):end));
FISm1eV(i) = simpsons(Te(qu_i(i):end),Tactivity.*Fm1eVN(qu_i(i):end));
FIS1eV(i) =  simpsons(Te(qu_i(i):end),Tactivity.*F1eVN(qu_i(i):end));
end

% add background and convert into counts
Bkg = 0.290;
SISm1eV = (SISm1eV+Bkg).*Time; SIS1eV = (SIS1eV+Bkg).*Time;
FISm1eV = (FISm1eV+Bkg).*Time; FIS1eV = (FIS1eV+Bkg).*Time;

fig1 = figure('Renderer','openGL');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
stairs(qu,(SISm1eV-FISm1eV)./sqrt(SISm1eV),'LineWidth',3);
hold on;
stairs(qu,(SIS1eV-FIS1eV)./sqrt(SIS1eV),'LineWidth',3);
xlabel('retarding energy -18575 (eV)');
ylabel('norm. residuals');
leg = legend(sprintf('m_\\nu^2= -1eV^2'),sprintf('m_\\nu^2= 1eV^2')); legend boxoff;
leg.Location='southwest';
PrettyFigureFormat;
%title('integrated differential spectrum: Samak vs Benchmark')
fprintf('1 - ratio integral Samak/Benchmark (1eV) = %.14g \n',1-SIS1eV/FIS1eV);
fprintf('1 - ratio integral Samak/Benchmark (-1eV) = %.14g \n',1-SISm1eV/FISm1eV);
print(gcf,[savedir,'DScomparison_IntSpec.png'],'-dpng','-r450');
%% plot overlay

close all
fig1 = figure('Renderer','openGL');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
stairs(Te-18575,S1eVN,'LineWidth',3);
hold on;
stairs(Te-18575,F1eVN,'LineWidth',3);
hold off;
PrettyFigureFormat;
set(gca,'YScale','log')
xlabel('kin. energy - 18575 (eV)')
ylabel(sprintf('d\\Gamma/dE (cps/eV)'))
leg = legend('Samak','Benchmark'); legend boxoff
print(fig1,[savedir,'DScomparison_Overlay_1eV.png'],'-dpng','-r450');
%% plot residuals
fig2 = figure('Renderer','openGL');
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
plot(Te-18575,Sm1eVN-Fm1eVN ,...
    'LineWidth',3);
PrettyFigureFormat;
xlabel('retarding energy - 18575 (eV)')
ylabel(sprintf('\\Delta differential spectra'))
leg = legend('Samak - Benchmark'); legend boxoff
print(fig2,[savedir,'DScomparison_Residuals.png'],'-dpng','-r450');
%% plot ratio:
fig3 = figure('Renderer','openGL');
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
Ratio = (S1eVN)./(F1eVN);%
Ratio(isinf(Ratio))=1;
plot(Te-18575,Ratio,...
    'LineWidth',3);
PrettyFigureFormat;
leg = legend('Samak/Benchmark'); legend boxoff
xlabel('retarding energy - 18575 (eV)')
ylabel(sprintf('ratio differential spectra'))
%ylim([1-5e-06,1+2e-06])
xlim([-40,-4])
grid on
print(fig3,[savedir,'DScomparison_1eVRatio_Zoom.png'],'-dpng','-r450');

