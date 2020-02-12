%
% Get Feldman Cousins Confidence Interval
% 
% T. Lasserre
% Last Updated: 18/07/2019
%

CentralValue  = -0.98;
OneSigmaError = 0.975;
%CentralValue  = -0.257;
%OneSigmaError = 0.94;

% FC Data https://arxiv.org/pdf/physics/9711021.pdf - TABLE XII.
FCsigma     = -3.0:0.1:3.0;
FCsigma_fine= -3.0:0.04:3.0;
FC90low    = [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.02 0.12 0.22 0.31 0.38 0.45 0.51 0.58 0.65 0.72 0.79 0.87 0.95 1.02 1.11 1.19 1.28 1.37];
FC90up     = [0.26 0.27 0.28 0.29 0.30 0.32 0.33 0.34 0.36 0.38 0.40 0.43 0.45 0.48 0.52 0.56 0.60 0.64 0.70 0.75 0.81 0.88 0.95 1.02 1.10 1.18 1.27 1.36 1.45 1.55 1.64 1.74 1.84 1.94 2.04 2.14 2.24 2.34 2.44 2.54 2.64 2.74 2.84 2.94 3.04 3.14 3.24 3.34 3.44 3.54 3.64 3.74 3.84 3.94 4.04 4.14 4.24 4.34 4.44 4.54 4.64];
FC90LowI   = @(x) interp1(FCsigma,FC90low,x);
FC90UpI    = @(x) interp1(FCsigma,FC90up,x);

% Conversion to FC Intervals on m 90%CL
R      = CentralValue./OneSigmaError;
Low90  = FC90LowI(R);
High90 = FC90UpI(R);

% Neutrino Mass Feldman-Cousins Limit 
LowerLimitL90m2 =  @(r,sigma) (FC90LowI(r).*sigma);
UpperLimitL90m2 =  @(r,sigma) (FC90UpI(r).*sigma);

% Neutrino Mass Feldman-Cousins Limit 
LowerLimitL90m =  @(r,sigma) sqrt(FC90LowI(r).*sigma);
UpperLimitL90m =  @(r,sigma) sqrt(FC90UpI(r).*sigma);

%% Belt m2
myMainTitle=[sprintf('KATRIN KNM1 Feldman Cousins Confidence Interval 90%%')];
maintitle=myMainTitle;
savefile=sprintf('plots/FC_KNM1.png');
fig1 = figure('Name','KATRIN KNM1 Feldman Cousins Confidence Interval 90%%','NumberTitle','off','rend','painters','pos',[10 10 1200 800]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
plot(FCsigma_fine,LowerLimitL90m2(FCsigma_fine,OneSigmaError),'color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',7);
hold on
plot(FCsigma_fine,UpperLimitL90m2(FCsigma_fine,OneSigmaError),'color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',7);
lm=line([CentralValue CentralValue],[0 5],'Color',rgb('SteelBlue'),'LineWidth',2,'LineStyle','--');
lu=line([-3 3],[UpperLimitL90m2(R,OneSigmaError) UpperLimitL90m2(R,OneSigmaError)],'Color',rgb('SteelBlue'),'LineWidth',2,'LineStyle','--');
hold off
xlim([-2 3]);
xlabel('Measured m^2 (eV^2)');
ylabel('True m^2 (eV^2)');
grid on
PrettyFigureFormat
set(gca,'FontSize',28);

%% Belt m
myMainTitle=[sprintf('KATRIN KNM1: Feldman Cousins Confidence Interval')];
maintitle=myMainTitle;
savefile=sprintf('plots/FC_KNM1.png');
fig1 = figure('Name','KATRIN KNM1 Feldman Cousins Confidence Interval 90%%','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=30;a.FontWeight='bold';
h=plot(FCsigma_fine,LowerLimitL90m(FCsigma_fine,OneSigmaError),'color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',7);
hold on
h=plot(FCsigma_fine,UpperLimitL90m(FCsigma_fine,OneSigmaError),'color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',7);
lm=line([CentralValue CentralValue],[0 UpperLimitL90m(R,OneSigmaError) ],'Color',rgb('SteelBlue'),'LineWidth',4,'LineStyle','--');
lu=line([-3 CentralValue],[UpperLimitL90m(R,OneSigmaError) UpperLimitL90m(R,OneSigmaError)],'Color',rgb('SteelBlue'),'LineWidth',4,'LineStyle','--');
hold off
xlim([-2 3]);
ylim([0 2.5]);
xlabel('Measured m^2 (eV^2)');
ylabel('True m (eV)');
leg = legend([h lu],...
    '90% Confidence Interval',...
    sprintf('m_\\beta < %.2f eV',UpperLimitL90m(R,OneSigmaError)));
leg.Color = 'none'; legend boxoff; leg.FontSize = 30;
grid on
PrettyFigureFormat
set(gca,'FontSize',28);
export_fig(gcf,'plots/KNM1_BlindedFSD_FCBelt90CL','-q101','-m3');