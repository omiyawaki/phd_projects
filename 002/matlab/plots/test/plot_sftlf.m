function plot_sftlf(type, par)
% plot sicedow depth
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read clear sky albedo data

    sftlf = squeeze(nanmean(sftlf,1)); % zonal mean

    % mon x lat of sicedow depth
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, sftlf, '-k');
    xlabel('Latitude (deg)'); ylabel('Land fraction (\%)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'tickdir', 'in');
    print(sprintf('%s/sftlf/sftlf_lat', par.plotdir), '-dpng', '-r300');
    close;

end
