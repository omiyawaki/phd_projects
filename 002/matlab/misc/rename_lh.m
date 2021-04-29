function lh = rename_lh(type, stf)

    if any(strcmp(type,{'era5', 'era5c', 'erai'}))
        lh = -stf.slhf;
    elseif any(strcmp(type,{'hahn'}))
        lh = stf.LHFLX;
    elseif any(strcmp(type,{'merra2', 'merra2c'}))
        lh = stf.EFLUX;
    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
        lh = stf.hfls;
    elseif strcmp(type, 'echam')
        lh = -stf.ahfl;
    end

end
