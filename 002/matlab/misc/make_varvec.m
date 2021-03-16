function var_vec = make_varvec(type, fw)

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c');
        if any(strcmp(fw, {'mse_old', 'dse_old'}))
            var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        else
            var_vec = {'sshf', 'slhf', 'cp', 'lsp', 'e', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        end
    elseif strcmp(type, 'hahn')
        var_vec = {'LHFLX', 'SHFLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    elseif strcmp(type, 'merra2')
        var_vec = {'EFLUX', 'HFLUX', 'PRECCON', 'PRECTOT', 'EVAP', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    elseif any(strcmp(type, {'jra55'}))
        var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    elseif any(strcmp(type, {'gcm'}))
        if any(strcmp(fw, {'mse_old', 'dse_old'}))
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
        else
            var_vec = {'hfls', 'hfss', 'pr', 'evspsbl', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc', 'tend'};
        end
    elseif strcmp(type, 'echam')
        var_vec = {'ahfl', 'ahfs', 'aprc', 'aprl', 'evap', 'lw', 'sw', 'rtoa', 'olr', 'lwsfc', 'swsfc'};
    end

end