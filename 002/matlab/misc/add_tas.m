function temp = add_tas(temp, lo, la, mo, type, srfc) 

    % adds surface data to vertical temperature profile

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        temp(end) = srfc.t2m(lo,la,mo); % surface is the last element in era plev
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        temp(1) = srfc.T2M(lo,la,mo);
    elseif strcmp(type, 'jra55')
        temp(end) = srfc.tas(lo,la,mo); % surface is last element in jra
    elseif strcmp(type, 'gcm')
        temp(1) = srfc.tas(lo,la,mo);
    elseif contains(type, 'hahn') % surface is last element in hahn data
        temp(end) = srfc.TREFHT(lo,la,mo);
    elseif contains(type, 'echam')
        temp(1) = srfc.temp2(lo,la,mo);
    end

end
