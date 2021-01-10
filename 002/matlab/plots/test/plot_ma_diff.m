function plot_ma_diff(type, par)
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        load(sprintf('/project2/mas1/miyawaki/projects/002/dama/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for plev_eval = 500; % evaluation plev in hPa
            diff = permute(ta.(land) - ma.(land).ta, [3 1 2]); % bring plev to front
            diff = squeeze(interp1(plev, diff, plev_eval)); % evaluate difference at plev_eval
            % mon x lat of temperature difference between model and moist adiabat
            figure(); clf; hold all;
            cmp = colCog(40);
            colormap(cmp);
            imagesc([1:12], [lat(1) lat(end)], diff);
            caxis([-20 20]);
            xlabel('Month'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-20 20], 'ytick', [-20:5:20], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$(T - T_m)_{%g \\,\\mathrm{hPa}}$ (K)', plev_eval));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ma_diff/plev_%g/%s/ma_diff_mon_lat', par.plotdir, plev_eval, land), '-dpng', '-r300');
            close;
        end
    end
end
