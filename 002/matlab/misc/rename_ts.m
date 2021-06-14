function tas = rename_tas(type, srfc) 

    if any(strcmp(type,{'era5', 'era5c', 'erai'}))
        tas = srfc.skt;
    elseif any(strcmp(type,{'hahn'}))
        tas = srfc.TS;
    elseif any(strcmp(type,{'merra2', 'merra2c'}))
        tas = srfc.TS;
    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
        tas = srfc.ts;
    elseif strcmp(type, 'echam')
        tas = srfc.tsurf;
    end

end
