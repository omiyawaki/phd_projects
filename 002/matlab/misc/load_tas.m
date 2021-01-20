function tas = load_tas(srfc, type, par) 
 
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        tas = srfc.t2m;
    elseif strcmp(type, 'merra2')
        tas = srfc.T2M;
    elseif strcmp(type, 'jra55')
        tas = srfc.tas; 
    elseif strcmp(type, 'gcm')
        tas = srfc.tas;
    elseif contains(type, 'echam')
        tas = srfc.temp2;
    end
    
end
