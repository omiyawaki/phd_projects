function ts = load_ts(srfc, type, par) 
 
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ts = srfc.t2m;
    elseif strcmp(type, 'merra2')
        ts = srfc.T2M;
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        ts = srfc.ts;
    elseif contains(type, 'echam')
        ts = srfc.tsurf;
    elseif contains(type, 'hahn')
        ts = srfc.TS;
    end
    
end
