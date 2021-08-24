function plot_alb(type, par)
% plot albedoow depth
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, ~] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/albedo.mat', prefix)); % read clear sky albedoedo data

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);
    plotdir = make_plotdir(type, par);

    albedo = squeeze(nanmean(albedo,1)); % zonal mean

    % mon x lat of albedoow depth
    figure(); clf; hold all;
    % cmp = colCog(20);
    colormap('gray');
    [C,h] = contourf(mesh_lat, mesh_mon, albedo, [0:0.05:1], 'linecolor', 'none');
    [C,h] = contour(mesh_lat, mesh_mon, albedo, [0.5:0.1:0.9], 'linecolor', 'w');
    clabel(C, h, [0.5:0.1:0.9], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([0 1]);
    ylabel('Latitude (deg)');
    make_title_type(type, par)
    % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    %     title(sprintf('%s', upper('ERA5')));
    % elseif strcmp(type, 'erai')
    % elseif strcmp(type, 'echam')
    %     title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s', par.model));
    % end
    cb = colorbar('limits', [0 1], 'ytick', [0:0.1:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Surface albedo (unitless)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/alb/albedo_mon_lat', plotdir), '-dpng', '-r300');
    close;

end
