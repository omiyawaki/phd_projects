function ps = rename_ps(type, srfc) 

    if any(strcmp(type,{'era5', 'era5c', 'erai'}))
        ps = srfc.sp;
    elseif any(strcmp(type,{'hahn'}))
        ps = srfc.PS;
    elseif any(strcmp(type,{'merra2', 'merra2c'}))
        ps = srfc.PS;
    elseif any(strcmp(type, {'gcm', 'jra55'}))
        ps = srfc.ps;
    elseif strcmp(type, 'echam')
        ps = srfc.aps;
    end

end
