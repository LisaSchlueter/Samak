% Script to Test Shape Only Analysis with Simulated Data
RunList = 'StackCD100all';
exclDataStart = 7;
nSamples = 1000;

R = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','exclDataStart',exclDataStart,'chi2','chi2CM');
RShape = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','exclDataStart',exclDataStart,'chi2','chi2CMShape');


save_name = sprintf('./results/Results_TestShapeOnly_%s_%.0feVrange_%.0fSamples.mat',...
    RunList,18575-R.ModelObj.qU(exclDataStart),nSamples);
if exist(save_name,'file')
    %Check if already computed
    load(save_name);
else
    %Do Computation
    TBDIS_Sim = mvnrnd(R.ModelObj.TBDIS,R.ModelObj.TBDIS',nSamples); %data with stat fuct
    parCM = zeros(6,nSamples);
    errCM = zeros(6,nSamples);
    chi2minCM = zeros(nSamples,1);
    parCMShape = zeros(6,nSamples);
    errCMShape = zeros(6,nSamples);
    chi2minCMShape = zeros(nSamples,1);
    
    progressbar('Computing...')
    for i=1:nSamples
        progressbar(i/nSamples);
        R.RunData.TBDIS = TBDIS_Sim(i,:)';
        R.Fit;
        parCM(:,i) = R.FitResult.par;
        errCM(:,i) = R.FitResult.err;
        chi2minCM(i) = R.FitResult.chi2min;
        
        RShape.RunData.TBDIS = TBDIS_Sim(i,:)';
        RShape.Fit;
        parCMShape(:,i) = RShape.FitResult.par;
        errCMShape(:,i) = RShape.FitResult.err;
        chi2minCMShape(i) = RShape.FitResult.chi2min;
    end
    dof = R.FitResult.dof;
    Q_i =RShape.ModelObj.Q_i;
    BKG_RateSec_i = R.ModelObj.BKG_RateSec_i;
    
    save(save_name,'chi2minCM', 'chi2minCMShape','parCM','errCM',...
        'parCMShape','parCM','nSamples','R','RShape','dof','BKG_RateSec_i','Q_i');
end
%% plot
f18= figure(18);
set(f18, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
subplot(2,4,1);
h1 = histfit(parCM(2,:)+Q_i,10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',wmean(parCM(2,:)+Q_i,errCM(2,:))));
xlabel('Endpoint CM (eV)');
PrettyFigureFormat;

subplot(2,4,5);
h1 = histfit(parCMShape(2,:)+Q_i,10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',wmean(parCMShape(2,:)+Q_i,errCMShape(2,:))))
xlabel('Endpoint Shape (eV)');
PrettyFigureFormat;

subplot(2,4,2);
h1 = histfit(1e3*(parCM(3,:)+BKG_RateSec_i),10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',1e3*(wmean(parCM(3,:),errCM(3,:))+BKG_RateSec_i)));
xlabel('Background CM (mcps)')
PrettyFigureFormat;

subplot(2,4,6);
h1 = histfit(1e3*(parCMShape(3,:)+BKG_RateSec_i),10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',1e3*(wmean(parCMShape(3,:),errCMShape(3,:))+BKG_RateSec_i)))
xlabel('Background Shape (mcps)')
PrettyFigureFormat;

subplot(2,4,3);
h1 = histfit(parCM(4,:)+1,10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',wmean(parCM(4,:),errCM(3,:))+1));
xlabel('Normalization CM (mcps)')
PrettyFigureFormat;

subplot(2,4,7);
h1 = histfit(parCMShape(3,:)+1,10,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('mean = %.3f',wmean(parCMShape(4,:),errCMShape(3,:))+1));
xlabel('Normalization Shape (mcps)')
PrettyFigureFormat;

subplot(2,4,4);
h1 = histogram(chi2minCM,'FaceColor',rgb('CadetBlue'));
hold on;
%x =(min(chi2minCM):max(chi2minCM))';
%y = repmat(dof,numel(x),1);
%pchi2 = plot(x,chi2pdf(x,y)*h1.BinWidth*nSamples,'LineWidth',3,'Color',rgb('Goldenrod'));
hold off;
xlabel(sprintf('\\chi^2 CM (%.0f dof)',dof));
legend(sprintf('mean = %.3f',mean(chi2minCM)));
PrettyFigureFormat;

subplot(2,4,8);
h1 = histogram(chi2minCMShape,'FaceColor',rgb('CadetBlue'));
%hold on;
%x =(min(chi2minCMShape):max(chi2minCMShape))';
%y = repmat(dof,numel(x),1);
%pchi2 = plot(x,chi2pdf(x,y)*h1.BinWidth*nSamples,'LineWidth',3,'Color',rgb('Goldenrod'));
%hold off;
xlabel(sprintf('\\chi^2 Shape (%.0f dof)',dof));
legend(sprintf('mean = %.3f',mean(chi2minCMShape)));
PrettyFigureFormat;

maintitle = (sprintf('Test of Shape Only Analysis with Simulation (stat Fluct.) \n MTD: %s',RunList));
 a=annotation('textbox', [0 0.91 1 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=18;a.FontWeight='bold';
    

if ~exist('./plots/pdf','dir')
    mkdir ./plots/pdf/
    mkdir ./plots/png/
end
print(f18,strrep(save_name,'results/','plots/png/'),'-dpng');
publish_figurePDF(f18,strrep(save_name,'results/','plots/pdf/'));