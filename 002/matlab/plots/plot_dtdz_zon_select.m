function plot_dtdz_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/dtdz_mon_lat.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/malr_mon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/dtdz_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/malr_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/dtdz_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/malr_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/dtdz_mon_lat.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/malr_mon_lat.mat', type, par.lat_interp));
        plev = 1:47;
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 7]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==7; mon_str = 'July'; end;

            dtdz_mon.(land) = squeeze(dtdz.(land)(:,month,:));
            dtmdz_si_mon.(land) = squeeze(dtmdz_si.(land)(:,month,:));

            dtdz_sp = interp1(grid.dim3.lat, dtdz_mon.(land), -85); % sounding at -85 S
            dtdz_np = interp1(grid.dim3.lat, dtdz_mon.(land), 85); % sounding at 85 N
            dtdz_smid = interp1(grid.dim3.lat, dtdz_mon.(land), -45); % sounding at -45 S
            dtdz_nmid = interp1(grid.dim3.lat, dtdz_mon.(land), 45); % sounding at 45 N
            dtdz_eq = interp1(grid.dim3.lat, dtdz_mon.(land), 0); % sounding at equator

            dtmdz_si_sp = interp1(grid.dim3.lat, dtmdz_si_mon.(land), -85); % sounding at -85 S
            dtmdz_si_np = interp1(grid.dim3.lat, dtmdz_si_mon.(land), 85); % sounding at 85 N
            dtmdz_si_smid = interp1(grid.dim3.lat, dtmdz_si_mon.(land), -45); % sounding at -45 S
            dtmdz_si_nmid = interp1(grid.dim3.lat, dtmdz_si_mon.(land), 45); % sounding at 45 N
            dtmdz_si_eq = interp1(grid.dim3.lat, dtmdz_si_mon.(land), 0); % sounding at equator

            figure(); clf; hold all;
            h_np = plot(dtdz_np, grid.dim3.si, 'color', par.blue);
            h_nmid = plot(dtdz_nmid, grid.dim3.si, 'color', par.orange);
            h_nmid_ma = plot(dtmdz_si_nmid, grid.dim3.si, ':', 'color', par.orange);
            h_eq = plot(dtdz_eq, grid.dim3.si, 'color', par.maroon);
            h_eq_ma = plot(dtmdz_si_eq, grid.dim3.si, ':', 'color', par.maroon);
            h_smid = plot(dtdz_smid, grid.dim3.si, '--', 'color', par.orange);
            h_smid_ma = plot(dtmdz_si_smid, grid.dim3.si, ':', 'color', par.orange);
            h_sp = plot(dtdz_sp, grid.dim3.si, '--', 'color', par.blue);
            xlabel('$\Gamma$ (K/km)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s, %s', upper(type), land_text, mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s, %s', par.model, land_text, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s, %s', upper(type), land_text, mon_str));
            end
            legend([h_np h_nmid h_eq h_smid h_sp], '85 N', '45 N', 'Equator', '45 S', '85 S', 'location', 'northwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [-inf inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/dtdz_zon_sel/%s/%g/all', par.plotdir, land, month), '-dpng', '-r300');
            close;

        end
    end
end
