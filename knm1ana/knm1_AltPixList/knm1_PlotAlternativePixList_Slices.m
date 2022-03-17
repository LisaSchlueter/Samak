% plot

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
chi2 = 'chi2Stat';
NP = 1.064;
RecomputeFlag = 'OFF';
AltPixList ='Slice3_1';  % defines alternative pixel list
PlotFPD = 'OFF';
% label
savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
savename = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,AltPixList,DataType,strrep(freePar,' ',''),range,chi2,NP);


if exist(savename,'file')
    d = importdata(savename);
    fprintf('load %s \n',savename)
else
    fprintf('file not found %s \n',savename)
    return
end

AngleDeg_m = mean(d.SliceAngPos,2);
if strcmp(AltPixList,'Slice3_1')
    AngleDeg_m(1) = 0;
end

if any(d.mNuSqErr==0)
    %sometimes asymmetric uncertainty may not work properly (convergence problems...), use some average uncertainty....
    fprintf(2,'WARNING %.0f slices have convergence problems (->tiny uncertainty) \n',sum(d.mNuSqErr==0));
    d.mNuSqErr(d.mNuSqErr==0)= median(d.FitResult.err(:,1));
end
%% exlucde slices that didn't converge
InclIdx = d.mNuSqErr<2.5*median(d.mNuSqErr); %& d.mNuSqErr>0;

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);

%% FPD viewer: slices (sanity plot)
if strcmp(PlotFPD ,'ON')
    close all
    Mode = 'Num';
    
    Slices_plt = NaN.*ones(148,1);
    for i=1:size(d.SliceAngPos,1)
        if strcmp(Mode,'Deg')
            Slices_plt(d.PixList{i}) = AngleDeg_m(i);
        elseif strcmp(Mode,'Num')
            Slices_plt(d.PixList{i}) = i;
        end
    end
    [plotHandle, cbHandle] = FPDViewer(Slices_plt,'Label','ON','ReDrawSkeleton','ON');
    
    
    if strcmp(Mode,'Deg')
        colormap(hot)
        cbHandle.Label.String = sprintf('Mean azimuth angle (degree)');
        cbHandle.Label.FontSize = get(gca,'FontSize')+3;
    elseif strcmp(Mode,'Num')
        colormap(jet);
        cbHandle.delete;
        Theta_polar = -AngleDeg_m+90;
        Theta_polar(Theta_polar<0) = Theta_polar(Theta_polar<0)+360;
        for i=1:numel(Theta_polar)
            t = text(4.95*cos(deg2rad(Theta_polar(i))),4.95*sin(deg2rad(Theta_polar(i))),sprintf('%.0f^\\circ',AngleDeg_m(i)),'FontSize',16,'Color',rgb('Silver'),'FontWeight','bold','HorizontalAlignment','center');
        end
    end
    
    pname1 = sprintf('%sknm1_FPD_%s_%s.pdf',pltdir,AltPixList,Mode);
    export_fig(pname1);
    fprintf('save plot to %s \n',pname1);
end
%% FPD viewer: mNuSq
if strcmp(PlotFPD ,'ON')
    close all
    mNuSq_Pix = NaN.*ones(148,1);
    for i=1:size(d.SliceAngPos,1)
        % if InclIdx(i)
        mNuSq_Pix(d.PixList{i}) = d.mNuSq(i);
        % else
        %      mNuSq_Pix(d.PixList{i}) = -inf;
        % end
    end
    [plotHandle, cbHandle,AngleDeg] = FPDViewer(mNuSq_Pix,'Label','OFF','ReDrawSkeleton','ON');
    colormap(hot)
    
    cbHandle.Label.String = sprintf('{\\itm}_\\nu^2 (eV^{ 2})');
    cbHandle.Label.FontSize = get(gca,'FontSize')+3;
    cbHandle.Position(1) = 0.86;
    
    
    Theta_polar = -AngleDeg_m+90;
    Theta_polar(Theta_polar<0) = Theta_polar(Theta_polar<0)+360;
    for i=1:numel(Theta_polar)
        t = text(4.95*cos(deg2rad(Theta_polar(i))),4.95*sin(deg2rad(Theta_polar(i))),sprintf('%.0f^\\circ',AngleDeg_m(i)),'FontSize',16,'Color',rgb('Silver'),'FontWeight','bold','HorizontalAlignment','center');
    end
    
    pname2 = sprintf('%sknm1_FPD_%s_mNuSq.pdf',pltdir,AltPixList);
    export_fig(pname2);
    fprintf('save plot to %s \n',pname2);
end
%% FPD viewer: E0
% close all
% E0_Pix = NaN.*ones(148,1);
% for i=1:size(d.SliceAngPos,1)
%     if InclIdx(i)
%        E0_Pix(d.PixList{i}) = d.E0(i)-18574;
%     end
% end
% FPDViewer(E0_Pix,'Label','ON','ReDrawSkeleton','ON');
% colormap(hot)
%% m^2 as a function of pixel position
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);

%Diff = [diff(AngleDeg_m)',360-AngleDeg_m(end)]';

wmeanAll = wmean(d.mNuSq(InclIdx),1./d.mNuSqErr(InclIdx).^2);

%if strcmp(AltPixList,'Slice3')
Threshold = 128;
wmean1 = wmean(d.mNuSq(AngleDeg_m<=Threshold & InclIdx),1./d.mNuSqErr(AngleDeg_m<=Threshold & InclIdx).^2);
wmean2 = wmean(d.mNuSq(AngleDeg_m>Threshold & InclIdx),1./d.mNuSqErr(AngleDeg_m>Threshold & InclIdx).^2);
Errwmean1 = mean(d.mNuSqErr(AngleDeg_m<=Threshold & InclIdx))./sqrt(sum(AngleDeg_m<=Threshold & InclIdx));
Errwmean2 = mean(d.mNuSqErr(AngleDeg_m>Threshold & InclIdx)./sqrt(sum(AngleDeg_m>Threshold & InclIdx)));
if ~strcmp(AltPixList,{'Slice','Slice3_1'})
    pm1 = plot([min(min(d.SliceAngPos(find(AngleDeg_m<=Threshold),:))) max(max(d.SliceAngPos(find(AngleDeg_m<=Threshold),:)))],...
        wmean1.*ones(2,1),'LineWidth',3,'Color',rgb('Gold'));
    hold on;
    pm2 = plot([min(min(d.SliceAngPos(find(AngleDeg_m>Threshold),:))) max(max(d.SliceAngPos(find(AngleDeg_m>Threshold),:)))],wmean2.*ones(2,1),'LineWidth',3,'Color',rgb('red'));
elseif strcmp(AltPixList,'Slice3_1')
    %     pm1 = plot([0 105+7.5], wmean1.*ones(2,1),'LineWidth',3,'Color',rgb('Gold'));
    %     hold on;
    %     plot([345-7.5 360], wmean1.*ones(2,1),'LineWidth',3,'Color',rgb('Gold'));
    %     hold on;
    %     pm2 = plot([min(min(d.SliceAngPos(find(AngleDeg_m>Threshold),:)))-7.5 7.5+max(max(d.SliceAngPos(find(AngleDeg_m>Threshold),:)))],wmean2.*ones(2,1),'LineWidth',3,'Color',rgb('red'));
    
    fUniform = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/',...
        sprintf('knm1_AltRunList_%s_NP%2g_KNM1.mat',chi2,NP)];
    du = importdata(fUniform);
    
    wMeanAll =  wmean(d.mNuSq,1./d.mNuSqErr.^2);
    ErrMeanAll =  mean(d.mNuSqErr)./sqrt(numel(d.mNuSqErr));
    pu =   plot([-10 370], du.FitResult.par(1).*ones(2,1),':','LineWidth',3,'Color',rgb('Silver'));
    hold on;
end
%end

e1 = errorbar(AngleDeg_m(InclIdx),d.mNuSq(InclIdx),d.mNuSqErr(InclIdx),'k.','MarkerSize',20,'CapSize',0,'LineWidth',2);%'Color',rgb('DodgerBlue'));
hold on
if any(InclIdx==0)
    angle = AngleDeg_m(~InclIdx);
    mnu = d.mNuSq(~InclIdx);
    mnuerr = d.mNuSqErr(~InclIdx);
    ym = ylim;
    for i=1:numel(mnu)
        %  e1 = errorbar(AngleDeg_m(~InclIdx),d.mNuSq(~InclIdx),d.mNuSqErr(~InclIdx),'.','Color',rgb('Silver'),'MarkerSize',20,'CapSize',0,'LineWidth',1);
        t = text(5,ym(1)+5*i,sprintf('Slice at \\langle%.0f^\\circ\\rangle: {\\itm}^2_\\nu = (%.0f\\pm%.0f) eV^2',angle(i),mnu(i),mnuerr(i)),...
            'FontSize',16,'Color',rgb('Silver'));
    end
end
xlabel('Mean azimuth angle');
ylabel(sprintf('{\\itm}_\\nu^2 (eV^{ 2})'));
PrettyFigureFormat('FontSize',18);

xlim([-9 360]);
xtickformat('%.0f°')
xticks([0:45:360]);

if ~strcmp(AltPixList,{'Slice','Slice3_1'})
    leg = legend([pm1,pm2],sprintf('\\langle{\\itm}_\\nu^2\\rangle = %.1f \\pm %.1f eV^2',wmean1,Errwmean1),...
        sprintf('\\langle{\\itm}_\\nu^2\\rangle = %.1f \\pm %.1f eV^2',wmean2,Errwmean2));
    PrettyLegendFormat(leg);
    leg.NumColumns =2;
    leg.FontSize = get(gca,'FontSize')+2;
    leg.ItemTokenSize = [40,18];
    leg.Location = 'north';
elseif strcmp(AltPixList,{'Slice3_1'})
    leg = legend([e1,pu],sprintf('Pseudo-slice-wise fits'),...
        sprintf('Uniform fit'));
    PrettyLegendFormat(leg);
    leg.NumColumns =2;
    leg.FontSize = get(gca,'FontSize')+4;
    leg.ItemTokenSize = [25,18];
    leg.Location = 'north';
end

if strcmp(AltPixList,'Slice3')
    ylim([-17,10]);
elseif strcmp(AltPixList,'Slice3_1')
    e1.Color = rgb('DodgerBlue');
    ylim([-16,8]);
    leg.Location = 'southwest';
elseif strcmp(AltPixList,'Slice4')
    ylim([-10,10]);
elseif strcmp(AltPixList,'Slice2')
    ylim([-15,13]);
end



pname3 = sprintf('%sknm1_AltPixList_%s_mNuSq.pdf',pltdir,AltPixList);
export_fig(pname3);
fprintf('save plot to %s \n',pname3);

%% significance:
Diff = (wmean2-wmean1);
ErrDiff = sqrt(Errwmean1^2+Errwmean2^2);
%Sigma1 = abs(Diff)./mean([Errwmean1,Errwmean2]); % wrong
Sigma2 = abs(Diff)./ErrDiff;


%% Significance among slices
if strcmp(AltPixList,'Slice3_1')
    wMean123   =  wmean(d.mNuSq(1:3),1./d.mNuSqErr(1:3).^2);
    ErrMean123 =  mean(d.mNuSqErr(1:3))./sqrt(numel(d.mNuSqErr(1:3)));
    
    wMean23    =  wmean(d.mNuSq(2:3),1./d.mNuSqErr(2:3).^2);
    ErrMean23  =  mean(d.mNuSqErr(2:3))./sqrt(numel(d.mNuSqErr(2:3)));
    
    wMean4End    =  wmean(d.mNuSq(4:end),1./d.mNuSqErr(4:end).^2);
    ErrMean4End  =  mean(d.mNuSqErr(4:end))./sqrt(numel(d.mNuSqErr(4:end)));
   
    % AngleDeg_m;
    
    % difference and significance of difference
    Diff1    = abs(wMean123-wMean4End);
    Diff1Err = sqrt(ErrMean123.^2+ErrMean4End^2);
    Sig1 = Diff1./Diff1Err;
    
    fprintf(' --------------------------------------------------\n');
    fprintf(' --------------------------------------------------\n');
    fprintf('difference between pseudo-slices \n');
    fprintf('(1:3) vs. (4:8) = (%.1f +-%.1f) eV --> %.1fsigma \n',Diff1,Diff1Err,Sig1);
    fprintf(' --------------------------------------------------\n');
    Diff2    = abs(wMean23-wMean4End);
    Diff2Err = sqrt(ErrMean23.^2+ErrMean4End^2);
    Sig2 = Diff2./Diff2Err;
    fprintf('difference between pseudo-slices \n');
    fprintf('(2:3) vs. (4:8) = (%.1f +-%.1f) eV --> %.1fsigma \n',Diff2,Diff2Err,Sig2);
    fprintf(' --------------------------------------------------\n');
    fprintf(' --------------------------------------------------\n');
    
    %% chi2 for constant
     chi2const = sum((wmeanAll-d.mNuSq).^2./d.mNuSqErr.^2);
     pconst = 1-chi2cdf(chi2const,numel(d.mNuSq)-1);
end

%),d.mNuSqErr(InclIdx)
