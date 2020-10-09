WGTS_CD_MolPerCm2 = 5e17.*(0:0.1:2);


savename = [getenv('SamakPath'),'RelicNuBkg/EffTritiumMass/results/ISCProb_RhoD.mat'];
A = ref_RelicNuBkg_TDR;
if exist(savename,'file')
   load(savename);
else
    P0 = zeros(numel(WGTS_CD_MolPerCm2),1);
    for i=2:numel(WGTS_CD_MolPerCm2)
        Prob = A.ComputeISProb('WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2(i));
        P0(i) = Prob(1)/1e2;
    end
    MakeDir([getenv('SamakPath'),'RelicNuBkg/EffTritiumMass/results/']);
    save(savename,'P0','WGTS_CD_MolPerCm2');
end

if ~exist('mt_total','var')
    eps_Acceptance = (1-cos(asin(sqrt(A.WGTS_B_T/A.MACE_Bmax_T))))./2;
    Nt_total =2.*WGTS_CD_MolPerCm2.*pi*A.WGTS_FTR_cm^2*A.WGTS_epsT; % number of tritium atoms in source
    % N_t = Nt_total.*P0'.*eps_Acceptance;
    
    mt_total      = Nt_total.*A.M/2*1e3;
    mt_acceptance = mt_total.*eps_Acceptance;
    mt_scat       = mt_total.*eps_Acceptance.*P0';
    m_t_eff       = mt_total.*eps_Acceptance.*P0'.*A.TTNormGS;
    save(savename,'m_t_eff','mt_scat','mt_acceptance','mt_total','Nt_total','eps_Acceptance','-append');
end
%% plot
GetFigure;
colormap('winter')
ptot  = plot(WGTS_CD_MolPerCm2./5e17,mt_total*1e6,':','LineWidth',3,'Color',rgb('Red'));
hold on;
pa    = plot(WGTS_CD_MolPerCm2./5e17,mt_acceptance*1e6,'--','LineWidth',3,'Color',rgb('ForestGreen'));
pscat = plot(WGTS_CD_MolPerCm2./5e17,mt_scat*1e6,'-.','LineWidth',3,'Color',rgb('Orange'));
peff  = plot(WGTS_CD_MolPerCm2./5e17,m_t_eff*1e6,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
hold on;
plot(ones(100,1),linspace(1e-6,1e3,100),'k:','LineWidth',2);
PrettyFigureFormat('FontSize',20);
xlabel(sprintf('Rel. column density (5\\cdot10^{17} molecules\\cdotcm^{-2})'))
ylabel(sprintf('Effective Tritium mass (\\mug)'))
leg = legend([ptot,pa,pscat,peff],...
    sprintf(' {\\itm}_{tot}'),...
    sprintf(' {\\itm}_{tot}\\cdot \\epsilon_{acpt}'),...
    sprintf(' {\\itm}_{tot}\\cdot \\epsilon_{acpt}\\cdot\\epsilon_{scat}'),...
    sprintf(' {\\itm}_{tot}\\cdot \\epsilon_{acpt}\\cdot\\epsilon_{scat}\\cdot\\epsilon_{fsd}'),...
    'Location','southeast','EdgeColor',rgb('Silver'));
ylim([1 1e3])
xlim([0.1, max(WGTS_CD_MolPerCm2./5e17)])
set(gca,'YScale','log')
savepathplot = [getenv('SamakPath'),'RelicNuBkg/EffTritiumMass/plots/'];
MakeDir(savepathplot);
plotname = sprintf('%sISCProb_RhoD_log.png',savepathplot);
print(gcf,plotname,'-dpng','-r400');

