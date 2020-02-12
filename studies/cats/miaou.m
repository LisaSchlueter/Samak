function [stats] = miaou(varargin)

close all

% Parser
p = inputParser;
p.addParameter('RunList',0);
p.addParameter('chi2','chi2CM',@(x)ismember(x,{'chi2Stat', 'chi2CM', 'chi2CMFrac','chi2CMShape', 'chi2P','chi2Nfix'}));
p.addParameter('fitter','minuit',@(x)ismember(x,{'minuit','matlab'}));
p.addParameter('fixPar','');
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('displayFit','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

RunList       = p.Results.RunList;
chi2          = p.Results.chi2;
fitter        = p.Results.fitter;
fixPar        = p.Results.fixPar;
exclDataStart = p.Results.exclDataStart;
displayFit    = p.Results.displayFit;

%% 
if RunList == 0
[ n, datafiles , RunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,'string');
RunList=RunList(RunList>40531);
RunList=RunList(RunList<40693);
end
 
%% Fit All RunList
fprintf(2,'processing run...\n')
close all
%First Fit : Fixing Neutrino Mass
FT = MultiRunAnalysis('RunList',RunList,'chi2',chi2,'exclDataStart',exclDataStart,'fixPar','1 5 6');
FT.ComputeCM('Mode','Read');
FT.Fit('CATS','ON');
%Second Fit : Getting input from first fit and fixing All but endpoint
%FT = MultiRunAnalysis('RunList',RunList,'chi2',chi2,'exclDataStart',exclDataStart,'fixPar','1 3 4',...
%     'i_mnu',FT.FitResult.par(1),'i_Q',0,'i_B',FT.FitResult.par(3),'i_N',FT.FitResult.par(4));
%FT.Fit();

switch displayFit
    case 'ON'
        %P = PLOTC('Xdata',FT.RunData.qU,'Ydata',FT.RunData.TBDIS,'ModelObj',FT.ModelObj,'CovMat',FT.FitCM,'saveplot','ON');
        %P.PlotSpectrumAndResiduals();
        FT.PlotFit('saveplot','ON','Mode','Count','ResidualsFlag','Norm');
end












% % Indexing
% index =  -(FT.ModelObj.qU-FT.ModelObj.Q)';
% 
% % Define Number of Free Parameters
% switch  fixPar
%     case ''
%         nFreeFitPar        = 4;
%         MaskFreeFitPar     = [1:4];
%     case '1'
%         nFreeFitPar        = 3;
%         MaskFreeFitPar     = [2:4];
%     case '1 3 4'
%         nFreeFitPar        = 1;
%         MaskFreeFitPar     = [2];
% end


% 
% % Load Data
% D=FT.RunData.TBDIS(exclDataStart:end);
% 
% % Load Model
% ModBF=TestFuncJacob(FT.FitResult.par,FT.ModelObj);
% ModBF=ModBF(exclDataStart:end);
% 
% % Fit Parameters at Best Fit
% BestPar = FT.FitResult.par;
% 
% % Design Matrix
% myf = @(p) TestFuncJacob(p,FT.ModelObj);
% Jaco = ComputeDesignMatrix(myf,BestPar','epsilon',1e-3);
% Jaco = Jaco((exclDataStart:end),(MaskFreeFitPar));
% 
% % Plot Design Matrix
% figure(1); ribbon(Jaco);colorbar
% 
% % Plot Design Matrix Columns
% figure(11)
% fig11 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1000 800]);
% a=annotation('textbox', [0 0.91 1 0.1], ...
%     'String', '', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center');
% a.FontSize=24;a.FontWeight='bold';
% countsplot=0;
% for i=1:nFreeFitPar
% countsplot=countsplot+1;
% subplot(numel(BestPar(MaskFreeFitPar)),1,countsplot)
%     stairs(Jaco((exclDataStart:end),i)','LineWidth',2,'Color',rgb('CadetBlue'));
%     yleg=sprintf('par %g',MaskFreeFitPar(i));
%     xleg='qU bin';
%     xlabel(xleg);
%     ylabel(yleg);
%     set(gca,'yscale','log');
% PrettyFigureFormat;
% end
% 
% % Fit Covariance Matrix
% % CovMat    = diag(D);
% CovMat      = diag(D) + FT.FitCM_Obj.CovMat(exclDataStart:end,exclDataStart:end);
% 
% % Call Cats
% stats = cats( D - ModBF + Jaco*BestPar(MaskFreeFitPar)' , Jaco , 'Vinv' , inv(CovMat));
% 
% % Residual Plot
% figure(2)
% resplot(stats);
% PrettyFigureFormat;
% 
% % Leverage Plot
% figure(3)
% levplot(stats);
% PrettyFigureFormat;
% 
% % Cook Plot
% figure(4)
% cookdplot(stats);
% PrettyFigureFormat;
% 
% % summary plot
% figure(6)
% stdrlevplot(stats)
% PrettyFigureFormat;
