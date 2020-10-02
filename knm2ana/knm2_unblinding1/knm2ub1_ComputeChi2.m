function chi2 = knm2ub1_ComputeChi2(data,model,covmat,exclDataStart)

chi2 = ((data(exclDataStart:end) - model(exclDataStart:end))')* (covmat(exclDataStart:end,exclDataStart:end)  \ (data(exclDataStart:end) - model(exclDataStart:end)));
end