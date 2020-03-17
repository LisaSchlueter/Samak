% plot alternative pixel list results
nFits = 500;
freePar = 'E0 Bkg Norm';
DataType = 'Real';
Parameter = 'E0';
RunList = 'KNM2_Prompt';
range = 40;                % fit range in eV below endpoint
AltPixLists = {'Half','Azi'};
savedir = [getenv('SamakPath'),'knm2ana/knm2_AlternativeRunLists/results/'];

%% load random half
savenameRand = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nFits);
if exist(savenameRand,'file')
    load(savenameRand);
else
    fprintf('random half pix list not computed yet \n')
end
%% load Alt run lists
ParAlt = cell(numel(AltPixLists),1);
nSeg   = zeros(numel(AltPixLists),1);

for i=1:numel(AltPixLists)
    savename = sprintf('%sknm2_PixListAlt_%s_%s_%s_%.0feV.mat',...
        savedir,AltPixLists{i},DataType,strrep(freePar,' ',''),range);
    
    if exist(savename,'file')
        d = importdata(savename);
        switch Parameter
            case 'mNuSq'
                ParAlt{i} = d.mNuSq;
            case 'E0'
                ParAlt{i} = d.E0;
        end
        nSeg(i) = numel(d.E0);
    else
        fprintf('alternative pix list %s not computed yet \n',AltPixLists{i})
    end
end
%% get parameters
switch Parameter
    case 'mNuSq'
        y    = mNuSq;
        yErr = mNuSqErr;
        xStr = sprintf('{\\itm}_\\nu^2');
        Unit = sprintf('eV^2');
    case 'E0'
        y    = E0+Q_i;
        yErr = E0Err;
        xStr = sprintf('{\\itE}_0^{fit}');
        Unit = sprintf('eV');
end
%% plot random half
f1 = figure('Units','normalized','Position',[0.1,0.1,0.7,0.5]);
meanY =wmean(y,1./yErr.^2);
h1 =histogram(y-meanY,'FaceColor',rgb('SkyBlue'),'FaceAlpha',1);
xlabel(sprintf('%s - \\langle%s\\rangle (%s)',xStr,xStr,Unit));
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)
t = title(sprintf('%.0f stacked runs - uniform fit to %.0f random pixels',numel(RunList),size(PixList,2)),...
    'FontWeight','normal');
legStrRand = sprintf('\\langle%s\\rangle = %.2f %s , \\sigma =  %.1e %s',xStr,meanY,Unit,std(y),Unit);
hold on;
yLimits = ylim;
%% plot alternative 
y = linspace(0,1000,20);
legStr = {legStrRand};
for i=1:numel(AltPixLists)
    AltY = ParAlt{i};
    for n=1:nSeg(i)
        plot(AltY(n).*ones(20,1)-meanY,y,'-','LineWidth',2) 
    end
    if strcmp(AltPixLists{i},'Half')
        legStr = [legStr,'Inner half','Outer half'];
    elseif strcmp(AltPixLists{i},'Azi')
        legStr = [legStr,'Azi 1','Azi 2','Azi 3','Azi 4'];
    end
end

leg = legend(legStr{:});
leg.EdgeColor = rgb('Silver');
ylim(yLimits);
