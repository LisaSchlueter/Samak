
function is_Pv = ISProbInterp(varargin)

p = inputParser;
p.addParameter('WGTS_CD_MolPerCm2',4.2*1e17,@(x)isfloat(x));
p.addParameter('MACE_Bmax_T',4.23,@(x)isfloat(x));
p.addParameter('WGTS_B_T',2.52,@(x)isfloat(x));
p.addParameter('CDsigma',5*1e17,@(x)isfloat(x));
p.addParameter('NIS',7,@(x)isfloat(x));
p.addParameter('ISXsection',3*64e-22,@(x)isfloat(x));
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});
WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
ISXsection        = p.Results.ISXsection;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
WGTS_B_T          = p.Results.WGTS_B_T;
NIS               = p.Results.NIS;
SanityPlot        = p.Results.SanityPlot;

if WGTS_CD_MolPerCm2<(2*1e17)
    DataSet = 'Knm1';
else 
    DataSet = 'Knm2';
end
    
[RhoDSigma,Theta,ISProb] = Compute_InitISProbMeshGrid('DataSet',DataSet);


is_Pv = zeros(NIS+1,numel(squeeze(ISXsection)));
[X,Y] = meshgrid(RhoDSigma,Theta);
ThetaFun = @(bmax,bs)  asin(sqrt(bs./bmax));

for i=1:NIS+1
    is_Pv(i,:) = interp2(X,Y,squeeze(ISProb(i,:,:))',...
        ISXsection.*WGTS_CD_MolPerCm2,...
        ThetaFun(MACE_Bmax_T,WGTS_B_T),...
        'spline');
end

if strcmp(SanityPlot,'ON')
    
    Scattering = 0;
    
    figure('Units','normalized','Position',[0.1,0.1,0.45,0.63]);
    s = surf(X.*1e4,Y,squeeze(ISProb(Scattering+1,:,:))');
    s.LineStyle = 'none';
    xlabel(sprintf('\\rhod\\sigma'));
    ylabel(sprintf('\\theta_{max}'));
    zlabel(sprintf('%.0f Inel. scattering prob. (%%)',Scattering));
    PrettyFigureFormat('FontSize',24);
    c = colorbar;
    c.Label.String = sprintf('%.0f inel. scattering prob. (%%)',Scattering);
    c.Label.FontSize = get(gca,'FontSize');
    view(90,90);
    xlim([min(min(X)),max(max(X))].*1e4);
    ylim([min(min(Y)),max(max(Y))]);
    grid off
    set(gca,'XMinorTick','off');
    set(gca,'YMinorTick','off');
end
end