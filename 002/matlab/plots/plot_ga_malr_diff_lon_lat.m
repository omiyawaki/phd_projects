function plot_ga_malr_diff_lon_lat(type, par)
% plot inversion strength
    make_dirs_si_bl(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/%s/ga_malr_diff_si_lon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_dalr_bl_diff_si_lon_lat.mat', prefix_proc, par.lat_interp));
    if strcmp(type, 'gcm')
        load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    else
        landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;
    end

    [mesh_lat, mesh_lon] = meshgrid(grid.dim3.lon, lat);
    [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim3.lon, lat);

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};

            % MALR FT lon x lat of diff
            figure(); clf; hold all;
            cmp = colCog(12);
            colormap(cmp);
            contourf(mesh_lat, mesh_lon, ga_malr_diff_t.(land).(time)', -100:10:100, 'linecolor', 'none');
            % [C,h]=contour(mesh_lat, mesh_lon, ga_malr_diff_t.(land).(time)', -100:10:100, '-k');
            % clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex');
            if strcmp(type, 'gcm')
                contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
            else
                contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
            end
            caxis([-60 60]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
            end
            cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
            set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ga_malr_diff/%s/%s/ga_malr_diff_lon_lat', par.plotdir, land, time), '-dpng', '-r300');
            close;

            % DALR BL lon x lat of diff
            figure(); clf; hold all;
            cmp = colCog(10);
            colormap(cmp);
            contourf(mesh_lat, mesh_lon, ga_dalr_bl_diff_t.(land).(time)', -100:20:100, 'linecolor', 'none');
            % [C,h]=contour(mesh_lat, mesh_lon, ga_dalr_bl_diff_t.(land).(time)', -100:20:100, '-w');
            % clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
            if strcmp(type, 'gcm')
                contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
            else
                contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
            end
            caxis([-100 100]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
            end
            cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
            set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ga_malr_diff/%s/%s/ga_dalr_bl_diff_lon_lat', par.plotdir, land, time), '-dpng', '-r300');
            close;

        end

    end
end
