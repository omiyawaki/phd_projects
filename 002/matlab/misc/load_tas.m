function tas = load_tas(srfc, type, par) 
 
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tas = srfc.t2m;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        tas = srfc.tas;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
    elseif contains(type, 'hahn')
        tas = srfc.TREFHT;
    end
    
end
