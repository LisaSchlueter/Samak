function zwi = p_par_square(z,E_source,BSOURCE,theta_source)

ME      = 511000.0;           % eV

zwi = (E_source*E_source + 2.0*E_source*ME)*(1-sin(theta_source).*sin(theta_source).*bfield(z)./BSOURCE);
% zwi += deltaU(z)*deltaU(z) - 2*deltaU(z)*(E_source+ME);

end
