function plot_tas(type, par)
% plot tasow depth
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
    load(sprintf('%s/srfc.mat', prefix));

    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        tas = srfc.t2m; clear srfc; % read tas
    elseif strcmp(type, 'gcm')
        tas = srfc.tas; clear srfc; % read tas
    elseif strcmp(type, 'echam')
        tas = srfc.temp2; clear srfc; % read tas
        if strcmp(par.echam.clim, '20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
    end

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    tas = squeeze(nanmean(tas,1)); % zonal mean

    % mon x lat of tasow depth
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    % colormap(cmp);
    [C,h] = contour(mesh_lat, mesh_mon, tas, [100:10:400]);
    % [C,h] = contour(mesh_lat, mesh_mon, tas, [270 270], 'linecolor', 'w');
    clabel(C, h, [100 150 200 250 270 300], 'fontsize', 6, 'interpreter', 'latex');
    ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), echamtext));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    if strcmp(type, 'echam') & ~contains(par.echam.clim, 'echr')
        caxis([150 270]);
        cb = colorbar('limits', [150 270], 'ytick', [150:10:270], 'location', 'eastoutside');
    else
        cb = colorbar('location', 'eastoutside');
    end
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$T_{2\\,\\mathrm{m}}$ (K)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/tas/tas_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
