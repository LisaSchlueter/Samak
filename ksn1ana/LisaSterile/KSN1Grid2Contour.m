function [mnu4sq_contour, sin2t4_contour] = KSN1Grid2Contour(mnu4sq,sin2t4,chi2,chi2_ref,CL,varargin)
%% calculate contour from msq4 - sin2t4 grid at confidence level cl
% define delta chi2 for given confidence level
p = inputParser;
p.addParameter('Mode','New',@(x)ismember(x,{'New','Old'}));
p.addParameter('Npar',2,@(x)isfloat(x));
p.parse(varargin{:});
Mode = p.Results.Mode;
Npar = p.Results.Npar;
DeltaChi2 = GetDeltaChi2(CL,Npar);
if chi2_ref>min(min(chi2))
    chi2_ref = min(min(chi2));
end
%% new
switch Mode
    case 'New'
        nInter = 1e3;
        [X,Y] = meshgrid(mnu4sq(:,1),sin2t4(1,:));
        mNugrid    = repmat(logspace(log10(min(min(mnu4sq))),log10(max(max(mnu4sq))),nInter),nInter,1);
        sin2t4grid = repmat(logspace(log10(min(min(sin2t4))),log10(max(max(sin2t4))),nInter),nInter,1)';
        chi2grid = reshape(interp2(X,Y,chi2,mNugrid,sin2t4grid,'spline'),nInter,nInter );

        hold on;
        [M,c]=contour(sin2t4grid,mNugrid,chi2grid-chi2_ref,[DeltaChi2 DeltaChi2],...
            'LineWidth',3,'LineColor','k');
        c.delete;
        
        % remove information about contour
        ExclIndex = find(M(1,:)==DeltaChi2);
        nContour = numel( ExclIndex);
        
        sin2t4_contour = cell(nContour,1);
        mnu4sq_contour = cell(nContour,1);
        for i=1:nContour
            startI = ExclIndex(i)+1;
            if i==nContour % last contour
                stopI = size(M,2);
            else
                stopI = ExclIndex(i+1)-1;
            end
            sin2t4_contour{i} = M(1,startI:stopI);
            mnu4sq_contour{i} = M(2,startI:stopI);
        end
        
    case 'Old'
        %% old
        nGridStep = size(mnu4sq,1);
        sin2t4_contour = zeros(nGridStep,1);
        mnu4sq_contour = mnu4sq(:,1);
        
        
        for i=1:nGridStep
            if max(chi2(:,i)-chi2_ref)<DeltaChi2 || min(chi2(:,i)-chi2_ref)>DeltaChi2
                sin2t4_contour(i) = NaN; % if contour point is not in grid of interest
            else
                InterpIndex = find(chi2(:,i)-chi2_ref<DeltaChi2,1,'last'); % last crossing of delta chi^2
                if InterpIndex>=3
                    InterpIndex = InterpIndex-2; % add to more points
                end
                sin2t4_contour(i) = interp1(chi2(InterpIndex:end,i)-chi2_ref,sin2t4(i,InterpIndex:end),DeltaChi2,'spline');
            end
        end
        
        
        ExclLogic = ( abs(sin2t4_contour)>=1) | isnan(sin2t4_contour);
        sin2t4_contour(ExclLogic) = [];
        mnu4sq_contour(ExclLogic) = [];
        
        nGridStep = numel(mnu4sq_contour);
        %% sanity check
        SanityPlot = 'OFF';%OFF';
        if strcmp(SanityPlot,'ON')
            GetFigure;
            mnuIndex = 15;
            tStr1 = 'problem m4 =';
            tStr2 = 'good m4 =';
            for i=1:nGridStep
                if any((chi2(:,i)-chi2_ref)<-0.1)
                    PlotArg = {'-.','Color',rgb('FireBrick'),'LineWidth',2};
                    tStr1 = [tStr1,sprintf(' %.0f ,',sqrt(mnu4sq_contour(i)))];
                else
                    PlotArg = {'-','Color',rgb('DimGray'),'LineWidth',2};
                    tStr2 = [tStr2,sprintf(' %.0f ,',sqrt(mnu4sq_contour(i)))];
                end
                p = plot(sin2t4(i,:),chi2(:,i)-chi2_ref, PlotArg{:});
                xlabel(sprintf('|U_{e4}|^2'));
                ylabel(sprintf('\\Delta \\chi^2'));
                PrettyFigureFormat;
                hold on;
            end
            tStr = sprintf('%s eV \n%s eV',tStr1,tStr2);
            title(tStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
            ylim([-20 10])
            xlim([0 0.2])
        end
end
end