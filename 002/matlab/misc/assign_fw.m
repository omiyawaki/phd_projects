function fw_vec = assign_fw(type, par)

    if any(strcmp(type, {'era5', 'era5c', 'erai'})); fw_vec = par.era.fw;
    elseif strcmp(type, 'jra55'); fw_vec = par.jra55.fw;
    elseif strcmp(type, 'merra2'); fw_vec = par.merra2.fw;
    elseif strcmp(type, 'gcm'); fw_vec = par.gcm.fw;
    elseif contains(type, 'echam'); fw_vec = par.echam.fw; end;

end
