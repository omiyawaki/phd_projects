function plot_ga_malr_diff_mon_lat(type, par)
% plot inversion strength
make_dirs_si_bl(type, par)

% load data
% [~, ~, ~, lat, par] = load_flux(type, par);
if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
elseif strcmp(type, 'merra2')
    par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
elseif strcmp(type, 'gcm')
    par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
elseif strcmp(type, 'echam')
    par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
    par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
end
load(sprintf('%s/%s/si_bl_%g/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp, par.si_bl));
load(sprintf('%s/%s/si_bl_%g/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp, par.si_bl));

[mesh_lat, mesh_mon] = meshgrid([1:12], lat);

for l = {'lo', 'l', 'o'}; land = l{1};
    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
    elseif strcmp(land, 'l'); land_text = 'Land';
    elseif strcmp(land, 'o'); land_text = 'Ocean';
    end

    % mon x lat of diff
    figure(); clf; hold all;
    cmp = colCog(12);
    colormap(flipud(cmp));
    contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, 'linecolor', 'none');
    [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, '-w');
    contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
    clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-60 60]);
    xlabel('Month'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper(type), land_text));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, %s', par.model, land_text));
    end
    cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat', par.plotdir, par.si_bl, land), '-dpng', '-r300');
    close;

    % ANOT mon x lat of diff
    figure(); clf; hold all;
    cmp = colCog(12);
    colormap(flipud(cmp));
    contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, 'linecolor', 'none');
    contour(mesh_lat, mesh_mon, ga_malr_diff.(land), par.ga_thresh*[1 1], 'color', 'r', 'linewidth', 2);
    [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, '-w');
    contour(mesh_lat, mesh_mon, ga_malr_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
    clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-60 60]);
    xlabel('Month'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper(type), land_text));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, %s', par.model, land_text));
    end
    cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_malr_diff_mon_lat_anot', par.plotdir, par.si_bl, land), '-dpng', '-r300');
    close;

    % DALR BL mon x lat of diff
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(flipud(cmp));
    contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, 'linecolor', 'none');
    [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, '-w');
    contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
    clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-100 100]);
    xlabel('Month'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper(type), land_text));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, %s', par.model, land_text));
    end
    cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_dalr_bl_diff_mon_lat', par.plotdir, par.si_bl, land), '-dpng', '-r300');
    close;

    % ANOT DALR BL mon x lat of diff
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(flipud(cmp));
    contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, 'linecolor', 'none');
    contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), par.ga_bl_thresh*[1 1], 'color', 'r', 'linewidth', 2);
    [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, '-w');
    contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), [0 0], 'color', 0.75*[1 1 1]);
    clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-100 100]);
    xlabel('Month'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2') | strcmp(type, 'era5c')
        title(sprintf('%s, %s', upper(type), land_text));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, %s', par.model, land_text));
    end
    cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/si_bl_%g/%s/ga_dalr_bl_diff_mon_lat_anot', par.plotdir, par.si_bl, land), '-dpng', '-r300');
    close;

end
end
