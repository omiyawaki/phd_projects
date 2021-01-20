function olrcs = load_olrcs(radcs, type, par) 
 
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        olrcs = radcs.ttrc;
    elseif strcmp(type, 'merra2')
        olrcs = srfc.T2M;
    elseif strcmp(type, 'jra55')
        olrcs = srfc.olrcs; 
    elseif strcmp(type, 'gcm')
        olrcs = srfc.olrcs;
    elseif contains(type, 'echam')
        olrcs = radcs.traf0;
    end
    
end
