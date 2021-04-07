function read_snevap(type, par)
% read sea ice depth
    if strcmp(type, 'echam')

        file=dir(sprintf('/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/%s/BOT_%s_0020_39.nc', par.echam.clim, par.echam.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        alb=double(ncread(fullpath, 'snevap'));
        newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam/%s', par.echam.clim);
        if ~exist(newdir, 'dir'); mkdir(newdir); end
        filename='alb.mat';
        save(sprintf('%s/%s', newdir, filename), 'alb');

    elseif any(strcmp(type, {'era5', 'era5c', 'erai'}))

        % the data is originally reported as accumulated m (depth) over a day, so
        % divide by 86400 s and multiply by 1000 kg/m^3 to get the
        % conventional kg/m^2/s mass flux over the full day
        snevap = 1e3/86400*double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/snevap/%s_snevap_%s.ymonmean.nc', type, type, par.(type).yr_span), 'es'));
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/snevap.mat', type, par.(type).yr_span), 'snevap');

    end
end
