function tas = rename_tas(type, srfc) 

    if any(strcmp(type,{'era5', 'era5c', 'erai'}))
        tas = srfc.t2m;
    elseif any(strcmp(type,{'hahn'}))
        tas = srfc.TREFHT;
    elseif any(strcmp(type,{'merra2', 'merra2c'}))
        tas = srfc.T2M;
    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
        tas = srfc.tas;
    elseif strcmp(type, 'echam')
        tas = srfc.temp2;
    end

end
