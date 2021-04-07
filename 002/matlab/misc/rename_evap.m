function evap = rename_evap(type, hydro)

    if any(strcmp(type,{'era5', 'era5c', 'erai'}))
        evap = -hydro.e;
    elseif any(strcmp(type,{'merra2'}))
        evap = hydro.EVAP;
    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
        evap = hydro.evspsbl;
    elseif strcmp(type, 'echam')
        evap = -hydro.evap;
    end

end
