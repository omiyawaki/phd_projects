function var_vec = make_varvec(type, fw)

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c');
        % var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        var_vec = {'slhf', 'sshf', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
    elseif strcmp(type, 'hahn')
        var_vec = {'LHFLX', 'SHFLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        % var_vec = {'EFLUX', 'HFLUX', 'PRECCON', 'PRECTOT', 'EVAP', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        var_vec = {'EFLUX', 'HFLUX', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
    elseif any(strcmp(type, {'jra55'}))
        var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
    elseif any(strcmp(type, {'gcm'}))
        var_vec = {'hfls', 'hfss', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
    elseif strcmp(type, 'echam')
        var_vec = {'ahfl', 'ahfs', 'aprc', 'aprl', 'evap', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
    end

end