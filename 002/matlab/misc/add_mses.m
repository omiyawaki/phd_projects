function temp = add_mses(temp, lo, la, mo, type, mses) 

    % adds surface data to vertical temperature profile

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        temp(end) = mses(lo,la,mo); % surface is the last element in era plev
    elseif strcmp(type, 'merra2')
        temp(1) = mses(lo,la,mo);
    elseif strcmp(type, 'jra55')
        temp(end) = mses(lo,la,mo); % surface is last element in jra
    elseif strcmp(type, 'gcm')
        temp(1) = mses(lo,la,mo);
    elseif contains(type, 'hahn') % surface is last element in hahn data
        temp(end) = mses(lo,la,mo);
    elseif contains(type, 'echam')
        temp(1) = mses(lo,la,mo);
    end

end
