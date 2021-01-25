function fw_vec = assign_fw(type, par)

    if any(strcmp(type, {'era5', 'era5c', 'erai'}));
        fw_vec = par.era.fw;
    else
        fw_vec = par.(type).fw;
    end

end
