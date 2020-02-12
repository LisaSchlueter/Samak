addpath(genpath('../../../Samak2.0'));

stepsize = [1,0.5,0.1,0.05];
%stepsize = 0.1;
inicio = [(1:18)*1e3,18500,18540,18570];
inicio = 18545;
CF = InitKATRINE0_1pixeq();





for kk = 1:length(inicio)


for ii = 1:length(stepsize)
    CF.Te = (0.01:stepsize(ii):18580)';
    CF.nTe = length(CF.Te);
    CF.SetKinVariables();
    CF.ComputeTBDDS();
    FullDS = CF.TBDDS;
    FullDSInt(ii,kk) = simpsons(CF.Te,FullDS);
    FullDSIntT(ii,kk) = trapz(CF.Te,FullDS);
    TBDDSf = @(e) interp1(CF.Te,CF.TBDDS,e,'spline');
    [FullDSIntG(ii,kk),err(ii,kk)] = quadgk(TBDDSf,min(CF.Te),max(CF.Te),'RelTol',1.e-10);
    
    CF.Te = (inicio(kk):stepsize(ii):18580)';
    CF.nTe = length(CF.Te);
    CF.SetKinVariables();
    CF.ComputeTBDDS();
    
    E0DS = CF.TBDDS;
    E0DSInt(ii,kk) = simpsons(CF.Te,E0DS);
    E0DSIntT(ii,kk) = trapz(CF.Te,E0DS);
    [E0DSIntG(ii,kk),errrr(ii,kk)] = quadgk(TBDDSf,min(CF.Te),max(CF.Te),'RelTol',1.e-10);
    
    CumFrac(ii,kk) = E0DSInt(ii,kk)/FullDSInt(ii,kk);
    CumFracT(ii,kk) = E0DSIntT(ii,kk)/FullDSIntT(ii,kk);
    CumFracG(ii,kk) = E0DSIntG(ii,kk)/FullDSIntG(ii,kk);


end
end
