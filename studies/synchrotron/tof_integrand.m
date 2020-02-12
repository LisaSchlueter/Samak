function r=tof_integrand(z,E_source,BSOURCE,theta_source)

ME      = 510998.0;           % eV
CLIGHT  = 299792458;                % m per second

p_par2 = p_par_square(z,E_source,BSOURCE,theta_source);
b      = bfield(z);

if (p_par2>0)
    one_over_v = (E_source + ME)./sqrt(p_par2)/CLIGHT;
    
    deltaE_synchrotron = 0.39.*b.*b.*b./BSOURCE .* sin(theta_source).*sin(theta_source).*E_source;
    
    %fprintf('%f %f %f %f %f\n',z,p_par2, E_source, deltaE_synchrotron,one_over_v*deltaE_synchrotron);
    
    r = one_over_v.*deltaE_synchrotron;
    
    return
    
else
    %fprintf('%f %f %f\n',z,p_par2, E_source);
    return
end

end
