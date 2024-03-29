% old functions
function plot_temp_zon_select_old(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
        plev = 1:47;
    end

    lat_pole = 85;
    lat_mid = 50;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 7]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==7; mon_str = 'July'; end;

            tasi_mon.(land) = squeeze(tasi.(land)(:,month,:));
            masi_mon.(land) = squeeze(masi.(land).ta(:,month,:));

            tasi_sp = interp1(lat, tasi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tasi_np = interp1(lat, tasi_mon.(land), lat_pole); % sounding at lat_pole N
            tasi_smid = interp1(lat, tasi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tasi_nmid = interp1(lat, tasi_mon.(land), lat_mid); % sounding at lat_mid N
            tasi_eq = interp1(lat, tasi_mon.(land), 0); % sounding at equator

            masi_sp = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
            masi_np = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
            masi_smid = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
            masi_nmid = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
            masi_eq = interp1(lat, masi_mon.(land), 0); % sounding at equator

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np, grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m == 7
                h_np = plot(tasi_np, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid, grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid, grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_np h_nmid], sprintf('%g N', lat_pole), sprintf('%g N', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np, tasi_nmid]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp, grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_sp h_smid], sprintf('%g S', lat_pole), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp, tasi_smid]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;
        end
    end
end
function plot_va(type, par)
    make_dirs(type, par)

    % load data
    if strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        lat = grid.dim3.lat;
        plev = grid.dim3.plev/100;
        var = 'va';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_%s_r1i1p1_*.ymonmean.nc', par.model, var, par.model, par.gcm.clim));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        va = ncread(fullpath, var);
    end

    % zonal mean
    va_z = squeeze(nanmean(va,1));
    % annual mean
    va_zt = squeeze(nanmean(va_z,3));

    % meridional velocity, lat x height
    [mesh_p, mesh_lat] = meshgrid(lat, plev);
    figure(); clf; hold all;
    cmp = colCog(30);
    colormap(cmp);

    [C, h] = contour(mesh_p, mesh_lat, va_zt', -5:0.25:5);
    clabel(C, h, [-2:0.5:2], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-1 1]);
    xlabel('latitude (deg)'); ylabel('pressure (hPa)');
    title(sprintf('Northward Velocity (m/s)'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ydir', 'reverse', 'yscale', 'log', 'ylim', [100 1000], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/va/va_lat_p', par.plotdir), '-dpng', '-r300');
    close;
end
function plot_ga_diff(type, par)
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    plev = par.pa/100;
    % load(sprintf('%s/dtdz.mat', prefix)); % read model lapse rate
    % load(sprintf('%s/dtdzi.mat', prefix)); % read model lapse rate
    load(sprintf('%s/dtdzz.mat', prefix)); % read model lapse rate
    % load(sprintf('%s/dtmdz.mat', prefix)); % read moist adiabatic lapse rate
    load(sprintf('%s/dtmdzz.mat', prefix)); % read moist adiabatic lapse rate
    % load(sprintf('%s/dtmdz_zt.mat', prefix)); % read moist adiabatic lapse rate

    dtdz_zt = squeeze(nanmean(nanmean(dtdzz, 1), 4));
    dtmdz_zt = squeeze(nanmean(nanmean(dtmdzz, 1), 4));
    dtdz_jan = squeeze(nanmean(dtdzz(:,:,:,1), 1));
    dtmdz_jan = squeeze(nanmean(dtmdzz(:,:,:,1), 1));
    dtdz_jul = squeeze(nanmean(dtdzz(:,:,:,7), 1));
    dtmdz_jul = squeeze(nanmean(dtmdzz(:,:,:,7), 1));

    diff = (dtmdzz - dtdzz)./dtmdzz * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates

    % LAT x PLEV
    diff_zt = squeeze(nanmean(nanmean(diff, 1), 4)); % take zonal and annual mean
    diff_jan = squeeze(nanmean(diff(:,:,:,1), 1)); % take zonal and annual mean
    diff_jul = squeeze(nanmean(diff(:,:,:,7), 1)); % take zonal and annual mean

    % vertically average
    dtdz_ztv = permute(dtdz_zt, [2 1]); % bring height front
    dtdz_ztv = interp1(par.pa, dtdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtdz_ztv = squeeze(nanmean(dtdz_ztv,1)); % take vertical average
    dtmdz_ztv = permute(dtmdz_zt, [2 1]); % bring height front
    dtmdz_ztv = interp1(par.pa, dtmdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtmdz_ztv = squeeze(nanmean(dtmdz_ztv,1)); % take vertical average
    dtdz_janv = permute(dtdz_jan, [2 1]); % bring height front
    dtdz_janv = interp1(par.pa, dtdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtdz_janv = squeeze(nanmean(dtdz_janv,1)); % take vertical average
    dtmdz_janv = permute(dtmdz_jan, [2 1]); % bring height front
    dtmdz_janv = interp1(par.pa, dtmdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtmdz_janv = squeeze(nanmean(dtmdz_janv,1)); % take vertical average
    dtdz_julv = permute(dtdz_jul, [2 1]); % bring height front
    dtdz_julv = interp1(par.pa, dtdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtdz_julv = squeeze(nanmean(dtdz_julv,1)); % take vertical average
    dtmdz_julv = permute(dtmdz_jul, [2 1]); % bring height front
    dtmdz_julv = interp1(par.pa, dtmdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    dtmdz_julv = squeeze(nanmean(dtmdz_julv,1)); % take vertical average

    % lat of climatological and moist adiabatic lapse rate
    figure(); clf; hold all;
    plot(grid.dim3.lat(1:6:end), dtdz_ztv(1:6:end), 'xk', 'markersize', 4);
    plot(grid.dim3.lat(1:6:end), dtmdz_ztv(1:6:end), '^k', 'markersize', 4);
    xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    print(sprintf('%s/ga_diff/ga_comp_lat', par.plotdir), '-dpng', '-r300');
    close;

    figure(); clf; hold all;
    plot(grid.dim3.lat(1:6:end), dtdz_janv(1:6:end), 'xk', 'markersize', 4);
    plot(grid.dim3.lat(1:6:end), dtmdz_janv(1:6:end), '^k', 'markersize', 4);
    xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    print(sprintf('%s/ga_diff/ga_comp_lat_jan', par.plotdir), '-dpng', '-r300');
    close;

    figure(); clf; hold all;
    plot(grid.dim3.lat(1:6:end), dtdz_julv(1:6:end), 'xk', 'markersize', 4);
    plot(grid.dim3.lat(1:6:end), dtmdz_julv(1:6:end), '^k', 'markersize', 4);
    xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    print(sprintf('%s/ga_diff/ga_comp_lat_jul', par.plotdir), '-dpng', '-r300');
    close;

    [mesh_p, mesh_lat] = meshgrid(grid.dim3.lat, plev);

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_p, mesh_lat, diff_zt', [-50 -20 -10 0 10 20 50 100]);
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [200 1000], 'ytick', [200 300 400:200:1000], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_plev', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, diff_zt', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_plev_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, diff_jan', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_plev_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, diff_jul', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_plev_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of model lapse rate
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_p, mesh_lat, dtdz_zt', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [200 1000], 'ytick', [200 300 400:200:1000], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_plev', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of model lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, dtdz_zt', [0:9], 'color', 'k');
    clabel(C, h, [0 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_plev_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of model lapse rate 10 S to 75 N JANUARY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, dtdz_jan', [-5:9], 'color', 'k');
    clabel(C, h, [-3 -1 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_plev_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of model lapse rate 10 S to 75 N JULY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, dtdz_jul', [0:9], 'color', 'k');
    clabel(C, h, [0 1 2 3 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_plev_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of moist adiabatic lapse rate
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_p, mesh_lat, dtmdz_zt', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [200 1000], 'ytick', [200 300 400:200:1000], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_m_lat_plev', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of moist adiabatic lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_p, mesh_lat, dtmdz_zt', [5 6 7 8 9], 'color', 'k');
    clabel(C, h, [5 6 7 8 9], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [250 1000], 'ytick', [300:100:1000], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_m_lat_plev_sc', par.plotdir), '-dpng', '-r300');
    close;

    % MON X LAT
    diff = permute(diff, [3 1 2 4]); % bring height front
    diff_v = interp1(par.pa, diff, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    diff_vz = squeeze(nanmean(nanmean(diff_v,1),2)); % take vertical and zonal average

    % mon x lat of temperature difference between model and moist adiabat
    figure(); clf; hold all;
    cmp = colCog(30);
    colormap(cmp);
    imagesc([1:12], [grid.dim3.lat(1) grid.dim3.lat(end)], diff_vz);
    caxis([-60 60]);
    xlabel('Month'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-60 60], 'ytick', [-60:20:60], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_ga_si_diff(type, par) % sigma coordinates
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzzsi.mat', prefix)); % read model lapse rate
    load(sprintf('%s/dtmdzzsi.mat', prefix)); % read moist adiabatic lapse rate

    dtdz_zt = squeeze(nanmean(nanmean(dtdzzsi, 1), 4));
    dtmdz_zt = squeeze(nanmean(nanmean(dtmdzzsi, 1), 4));
    dtdz_jan = squeeze(nanmean(dtdzzsi(:,:,:,1), 1));
    dtmdz_jan = squeeze(nanmean(dtmdzzsi(:,:,:,1), 1));
    dtdz_jul = squeeze(nanmean(dtdzzsi(:,:,:,7), 1));
    dtmdz_jul = squeeze(nanmean(dtmdzzsi(:,:,:,7), 1));

    diff = (dtmdzzsi - dtdzzsi)./dtmdzzsi * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates

    % LAT x PLEV
    diff_zt = squeeze(nanmean(nanmean(diff, 1), 4)); % take zonal and annual mean
    diff_jan = squeeze(nanmean(diff(:,:,:,1), 1)); % take zonal and annual mean
    diff_jul = squeeze(nanmean(diff(:,:,:,7), 1)); % take zonal and annual mean

    % % vertically average
    % dtdz_ztv = permute(dtdz_zt, [2 1]); % bring height front
    % dtdz_ztv = interp1(par.pa, dtdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_ztv = squeeze(nanmean(dtdz_ztv,1)); % take vertical average
    % dtmdz_ztv = permute(dtmdz_zt, [2 1]); % bring height front
    % dtmdz_ztv = interp1(par.pa, dtmdz_ztv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_ztv = squeeze(nanmean(dtmdz_ztv,1)); % take vertical average
    % dtdz_janv = permute(dtdz_jan, [2 1]); % bring height front
    % dtdz_janv = interp1(par.pa, dtdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_janv = squeeze(nanmean(dtdz_janv,1)); % take vertical average
    % dtmdz_janv = permute(dtmdz_jan, [2 1]); % bring height front
    % dtmdz_janv = interp1(par.pa, dtmdz_janv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_janv = squeeze(nanmean(dtmdz_janv,1)); % take vertical average
    % dtdz_julv = permute(dtdz_jul, [2 1]); % bring height front
    % dtdz_julv = interp1(par.pa, dtdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtdz_julv = squeeze(nanmean(dtdz_julv,1)); % take vertical average
    % dtmdz_julv = permute(dtmdz_jul, [2 1]); % bring height front
    % dtmdz_julv = interp1(par.pa, dtmdz_julv, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % dtmdz_julv = squeeze(nanmean(dtmdz_julv,1)); % take vertical average

    % % lat of climatological and moist adiabatic lapse rate
    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_ztv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_ztv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat', par.plotdir), '-dpng', '-r300');
    % close;

    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_janv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_janv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s, January', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s, January', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat_jan', par.plotdir), '-dpng', '-r300');
    % close;

    % figure(); clf; hold all;
    % plot(grid.dim3.lat(1:6:end), dtdz_julv(1:6:end), 'xk', 'markersize', 4);
    % plot(grid.dim3.lat(1:6:end), dtmdz_julv(1:6:end), '^k', 'markersize', 4);
    % xlabel('Latitude (deg)'); ylabel('Lapse rate (K/km)')
    % legend('$\Gamma$', '$\Gamma_m$', 'location', 'southwest');
    % if strcmp(type, 'era5') | strcmp(type, 'erai')
    %     title(sprintf('%s, July', upper(type)));
    % elseif strcmp(type, 'gcm')
    %     title(sprintf('%s, July', par.model));
    % end
    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    % set(gca, 'xlim', [-10 80], 'xtick', [0:15:75], 'ylim', [0 9], 'ytick', [0:2:8]);
    % print(sprintf('%s/ga_diff/ga_comp_lat_jul', par.plotdir), '-dpng', '-r300');
    % close;

    [mesh_si, mesh_lat] = meshgrid(grid.dim3.lat, grid.dim3.si);

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_zt', [-50 -20 -10 0 10 20 50 100]);
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_zt', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_jan', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_jul', [-50 -20 -10 0 10 20 50 100], 'color', 'k');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_diff_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_zt', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_zt', [0:9], 'color', 'k');
    clabel(C, h, [0 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N JANUARY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_jan', [-5:9], 'color', 'k');
    clabel(C, h, [-3 -1 2 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate 10 S to 75 N JULY
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtdz_jul', [0:9], 'color', 'k');
    clabel(C, h, [0 1 2 3 4 5 6 7 8], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtmdz_zt', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_m_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtmdz_zt', [5 6 7 8 9], 'color', 'k');
    clabel(C, h, [5 6 7 8 9], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_diff/ga_m_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_ga_diff_mon_lat(type, par) % plot inversion strength
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/%s/ga_diff_si_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % mon x lat of diff
        figure(); clf; hold all;
        cmp = colCog(30);
        colormap(cmp);
        % imagesc([1:12], [lat(1) lat(end)], ga_diff.(land));
        contourf(mesh_lat, mesh_mon, ga_diff.(land), 'linecolor', 'none')
        caxis([-100 100]);
        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, %s', upper(type), land_text));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, %s', par.model, land_text));
        end
        cb = colorbar('limits', [-20 100], 'ytick', [-20:20:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_m(T_m) - \\Gamma)/\\Gamma_m(T_m)\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_diff/%s/ga_diff_mon_lat', par.plotdir, land), '-dpng', '-r300');
        close;

        % mon x lat of diff
        figure(); clf; hold all;
        cmp = colCog(30);
        colormap(cmp);
        % imagesc([1:12], [lat(1) lat(end)], ga_bl_diff.(land));
        contourf(mesh_lat, mesh_mon, ga_bl_diff.(land), 'linecolor', 'none')
        caxis([-100 100]);
        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, %s', upper(type), land_text));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, %s', par.model, land_text));
        end
        cb = colorbar('limits', [-20 100], 'ytick', [-20:20:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_diff/%s/ga_bl_diff_mon_lat', par.plotdir, land), '-dpng', '-r300');
        close;

    end
end
function plot_alb(type, par) % plot albedo
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/alb.mat', prefix)); % read clear sky albedo data
    load(sprintf('%s/albcs.mat', prefix)); % read clear sky albedo data

    alb = squeeze(nanmean(alb,1)); % zonal mean
    albcs = squeeze(nanmean(albcs,1)); % zonal mean

    % mon x lat of clear sky surface albedo
    figure(); clf; hold all;
    cmp = colCog(30);
    colormap(cmp);
    imagesc([1:12], [grid.dim2.lat(1) grid.dim2.lat(end)], alb);
    caxis([-1 1]);
    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 1], 'ytick', [0:0.1:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\alpha$ (unitless)'));
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/alb/alb_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

    % mon x lat of clear sky surface albedo
    figure(); clf; hold all;
    cmp = colCog(30);
    colormap(cmp);
    imagesc([1:12], [grid.dim2.lat(1) grid.dim2.lat(end)], albcs);
    caxis([-1 1]);
    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 1], 'ytick', [0:0.1:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$\\alpha$ (unitless)'));
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/alb/albcs_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_sn(type, par)
% plot snow depth
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sn.mat', prefix)); % read clear sky albedo data

    sn = squeeze(nanmean(sn,1)); % zonal mean

    % mon x lat of snow depth
    figure(); clf; hold all;
    cmp = colCog(30);
    colormap(flip(cmp));
    imagesc([1:12], [grid.dim2.lat(1) grid.dim2.lat(end)], sn);
    caxis([-1e-1 1e-1]);
    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 1e-1], 'ytick', 1e-1*[0:0.1:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Ice Depth (m)'));
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/sn/sn_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_inv_str(type, par) % plot inversion strength
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for si_eval = [0.8 0.85 0.9]; % evaluation sigma level
            diff = permute(tasi.(land), [3 1 2]);
            diff = squeeze(interp1(grid.dim3.si, diff, si_eval) - diff(1,:,:)); % evaluate difference between surface and si_eval
            % lat x lon of diff
            figure(); clf; hold all;
            cmp = colCog(30);
            colormap(cmp);
            imagesc([1:12], [lat(1) lat(end)], diff);
            caxis([-15 15]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-15 15], 'ytick', [-15:5:15], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$T_{\\sigma = %g} - T_s$ (K)', si_eval));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/inv_str/si_%g/%s/inv_str_lat_lon', par.plotdir, si_eval, land), '-dpng', '-r300');
            close;
        end
    end
end
function plot_inv_str_alt(type, par) % plot inversion strength (alt order of operations: calculate inv_str at every lat lon time)
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/inv_mon_lat.mat', type, par.lat_interp));
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/inv_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for isig = 1:length(par.si_eval); sig = par.si_eval(isig); % evaluation sigma level
            % lat x lon of inversion strength
            figure(); clf; hold all;
            cmp = colCog(30);
            colormap(cmp);
            imagesc([1:12], [lat(1) lat(end)], inv.(land)(:,:,isig));
            caxis([-15 15]);
            xlabel('Month'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-15 15], 'ytick', [-15:5:15], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$T_{\\sigma = %g} - T_s$ (K)', sig));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/inv_str/si_%g/%s/inv_str_mon_lat', par.plotdir, sig, land), '-dpng', '-r300');
            close;
        end
    end
end
function plot_flux_block(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        % % lat x mon w500
        % var_text = '$\omega500$';
        % figure(); clf; hold all;
        % cmp = colCog(20);
        % colormap(cmp);
        % imagesc([1 12], [lat(1) lat(end)], flux_z.(land).w500*864);
        % cb = colorbar('ticks', [-50:10:50], 'ticklabelinterpreter', 'latex');
        % ylabel(cb, '$\omega500$ (hPa/d)', 'interpreter', 'latex');
        % caxis([-50 50]);
        % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
        % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
        % xlabel('Month'); ylabel('Latitude (deg)');
        % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/flag/%s/0_w500_mon_lat', par.plotdir, land), '-dpng', '-r300');
        % close;

        % % lat x mon vas
        % var_text = '$v_{\,\mathrm{2\,m}}$';
        % figure(); clf; hold all;
        % cmp = colCog(20);
        % colormap(cmp);
        % imagesc([1 12], [lat(1) lat(end)], flux_z.(land).vas);
        % cb = colorbar('ticks', [-5:1:5], 'ticklabelinterpreter', 'latex');
        % ylabel(cb, '$v_{\,\mathrm{2\,m}}$ (m/s)', 'interpreter', 'latex');
        % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
        % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
        % caxis([-5 5]);
        % xlabel('Month'); ylabel('Latitude (deg)');
        % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/flag/%s/0_vas_mon_lat', par.plotdir, land), '-dpng', '-r300');
        % close;

        % % lat x mon P-E
        % var_text = '$P-E$';
        % figure(); clf; hold all;
        % cmp = colCog(20);
        % colormap(cmp);
        % imagesc([1 12], [lat(1) lat(end)], (flux_z.(land).pr-flux_z.(land).evspsbl)*86400);
        % cb = colorbar('ticks', [-5:1:5], 'ticklabelinterpreter', 'latex');
        % ylabel(cb, '$P-E$ (mm/d)', 'interpreter', 'latex');
        % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
        % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
        % caxis([-5 5]);
        % xlabel('Month'); ylabel('Latitude (deg)');
        % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/flag/%s/0_pe_mon_lat', par.plotdir, land), '-dpng', '-r300');
        % close;

        % % lat x mon convective vs large-scale precip
        % var_text = '$P_{\mathrm{ls}}/P_{\mathrm{c}}$';
        % figure(); clf; hold all;
        % cmp = colCog(20);
        % colormap(cmp);
        % imagesc([1 12], [lat(1) lat(end)], ((flux_z.(land).pr-flux_z.(land).prc)./flux_z.(land).prc));
        % cb = colorbar('ticks', [-20:4:20], 'ticklabelinterpreter', 'latex');
        % ylabel(cb, '$P_{\mathrm{ls}}/P_{\mathrm{c}}$ (unitless)', 'interpreter', 'latex');
        % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
        % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
        % caxis([-20 20]);
        % xlabel('Month'); ylabel('Latitude (deg)');
        % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/flag/%s/0_cp_mon_lat', par.plotdir, land), '-dpng', '-r300');
        % close;

        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            % % lat x lon w500
            % var_text = '$\omega500$';
            % figure(); clf; hold all;
            % cmp = colCog(20);
            % colormap(cmp);
            % imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).w500');
            % caxis([-nanmax(abs(flux_t.(land).(time).w500(:))) nanmax(abs(flux_t.(land).(time).w500(:)))]);
            % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
            % elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
            % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            % set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            % print(sprintf('%s/flag/%s/%s/w500_lat_lon', par.plotdir, land, time), '-dpng', '-r300');
            % close;

            % % lat x lon vas
            % var_text = '$v_{\,\mathrm{2\,m}}$';
            % figure(); clf; hold all;
            % cmp = colCog(20);
            % colormap(cmp);
            % imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).vas');
            % caxis([-nanmax(abs(flux_t.(land).(time).vas(:))) nanmax(abs(flux_t.(land).(time).vas(:)))]);
            % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
            % elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
            % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            % set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            % print(sprintf('%s/flag/%s/%s/vas_lat_lon', par.plotdir, land, time), '-dpng', '-r300');
            % close;

        end

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            % lat x mon dependence of RCE and RAE
            if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
            else var_text = '$\nabla \cdot F_m$'; end
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], flux_z.(land).res.(fw));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, var_text)); end;
            caxis([-150 150]);
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_div_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % R1 lat x mon dependence of RCE and RAE
            var_text = '$R_1$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(flipud(cmp));
            imagesc([1 12], [lat(1) lat(end)], flux_z.(land).r1.(fw));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-1 2]);
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_r1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            if strcmp(fw, 'mse2')
                % lat x mon dependence of RCE and RAE
                var_text = 'LW';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], flux_z.(land).lw);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
                caxis([-250 250]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_lw_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
                close;
            else
                % lat x mon dependence of RCE and RAE
                var_text = '$R_a$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], flux_z.(land).ra.(fw));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
                caxis([-150 150]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_ra_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
                close;
            end

            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                % lat x lon of RCE and RAE
                if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
                else var_text = '$\nabla \cdot F_m$'; end
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).res.(fw)');
                caxis([-150 150]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/div_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                var_text = '$R_1$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(flipud(cmp));
                imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).r1.(fw)');
                caxis([-2 2]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/r1_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                if strcmp(fw, 'mse2')
                    var_text = 'LW';
                    figure(); clf; hold all;
                    cmp = colCog(20);
                    colormap(cmp);
                    imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).lw');
                    caxis([-250 250]);
                    if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/lw_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                    close;
                else
                    var_text = '$R_a$';
                    figure(); clf; hold all;
                    cmp = colCog(20);
                    colormap(cmp);
                    imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], flux_t.(land).(time).ra.(fw)');
                    caxis([-150 150]);
                    if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text)); end;
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/ra_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                    close;
                end

            end % for time
        end % for mse dse
    end % for land
end
function plot_tend_comp(type, par)
    make_dirs(type, par)

    if strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    else
        error('This comparison only applies to ERA-Interim.');
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_zt.mat', prefix_proc, par.lat_interp));

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        figure(); clf; hold all;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        h1=plot(lat, flux_zt.lo.(time).TETEN);
        h2=plot(lat, flux_zt.lo.(time).tend);
        title('MSE tendency comparison')
        xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
        legend([h1 h2], 'DB13', 'ERA-I', 'location', 'eastoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/energy-flux/lo/%s/tend_comp', par.plotdir, time), '-dpng', '-r300');
        close;
    end
end
function plot_dr1(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            dr1 = flux_z.(land).r1.(fw) - repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);

            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, dr1, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, dr1, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_1$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % FRACTIONAL DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta R_1}{R_1}$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, dr1./flux_z.(land).r1.(fw), [-1e10 -1 -0.25:0.025:0.25 1 1e10], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, dr1./flux_z.(land).r1.(fw), [0 0], 'linecolor', 0.75*[ 1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta R_1}{R_1}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_frac_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 DEF COMP1 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\nabla\cdot F_m)}{R_a}$';
            delta_fm = flux_z.(land).res.(fw) - repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp1 = delta_fm./ann_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp1, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp1, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_m)}{R_a}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1c1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 DEF COMP2 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\nabla\cdot F_m}{R_a^2}\Delta R_a$';
            delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            ann_fm = repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            comp2 = -ann_fm./(ann_ra).^2.*delta_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp2, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp2, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\nabla\cdot F_m)}{R_a}-\frac{\nabla\cdot F_m}{R_a^2}\Delta R_a$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp1+comp2, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp1+comp2, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_m)}{R_a}-\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_c1+c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % REAL DELTA R1 - DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1 - \left(\frac{\Delta (\nabla\cdot F_m)}{R_a}-\frac{\nabla\cdot F_m}{R_a^2}\Delta R_a\right)$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, dr1 - (comp1+comp2), [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, dr1 - (comp1+comp2), [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_1 - \\left(\\frac{\\Delta (\\nabla\\cdot F_m)}{R_a}-\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta R_a\\right)$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_res_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            end % for time
        end % for mse dse
    end % for land
function plot_dr1_alt(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            dr1 = flux_z.(land).r1.(fw) - repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);

            % DELTA R1 ALT COMP1 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{(\mathrm{SH+LH})}{R_a^2}\Delta (\nabla\cdot F_m)$';
            delta_fm = flux_z.(land).res.(fw) - repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
            ann_stf = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp1 = -ann_stf./(ann_ra).^2.*delta_fm;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp1, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp1, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta (\\nabla\\cdot F_m)$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_c1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 ALT COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
            delta_stf = flux_z.(land).stf.(fw) - repmat(nanmean(flux_z.(land).stf.(fw),2),[1 12]);
            ann_fm = repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            comp2 = ann_fm./(ann_ra).^2.*delta_stf;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp2, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp2, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta (\\mathrm{SH+LH})$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\mathrm{SH+LH}}{R_a^2}\Delta (\nabla\cdot F_m)+\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, comp1+comp2, [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, comp1+comp2, [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta (\\nabla\\cdot F_m)+\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta (\\mathrm{SH+LH})$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_c1+c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % REAL DELTA 1 - DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1-\left(-\frac{\mathrm{SH+LH}}{R_a^2}\Delta (\nabla\cdot F_m)+\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})\right)$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, dr1-(comp1+comp2), [-1 -0.25:0.025:0.25 1], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, dr1-(comp1+comp2), [0 0], 'linecolor', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.25 0.25]);
            cb = colorbar('limits', [-0.25 0.25], 'ytick', [-0.25:0.05:0.25], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_1 \\left(-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta (\\nabla\\cdot F_m)+\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta (\\mathrm{SH+LH})\\right)$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_res_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

        end % for mse dse
    end % for land
end
function plot_dmse_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_eval_list = [-85 -45 0 45 85]; % Latitude to evaluate R1 seasonality

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
    %         folder = sprintf('%s/dmse/%s/0_lat_%g', par.plotdir, fw, lat_eval);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         ra = flux_z.lo.ra.(fw);
    %         ra_lat = interp1(grid.dim3.lat, ra, lat_eval);
    %         comp1 = sftlf*1e-2.*flux_z.l.ra.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.ra.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % RA lat x mon dependence of RCE and RAE
    %         var_text = '$R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], ra_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$R_a$', '$R_{a,\mathrm{\,L+O}}$', '$R_{a,\mathrm{\,L}}$', '$R_{a,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_ra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         ra_ann = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat_eval);
    %         dra = flux_z.lo.ra.(fw) - ra_ann;
    %         dra_lat = interp1(grid.dim3.lat, dra, lat_eval);
    %         ra_ann_l = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.ra.(fw) - ra_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         ra_ann_o = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.ra.(fw) - ra_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % DELTA RA lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], dra_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$\Delta(R_a)$', '$\Delta(R_{a,\mathrm{\,L+O}})$', '$\Delta(R_{a,\mathrm{\,L}})$', '$\Delta(R_{a,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res = flux_z.lo.res.(fw);
    %         res_lat = interp1(grid.dim3.lat, res, lat_eval);
    %         comp1 = sftlf*1e-2.*flux_z.l.res.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.res.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % RES lat x mon dependence of RCE and RAE
    %         var_text = '$\nabla\cdot F_m$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], res_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$\nabla\cdot F_m$', '$\nabla\cdot F_{m,\mathrm{\,L+O}}$', '$\nabla\cdot F_{m,\mathrm{\,L}}$', '$\nabla\cdot F_{m,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\nabla\\cdot F_m$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_res_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res_ann = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         res_ann_lat = interp1(grid.dim3.lat, res_ann, lat_eval);
    %         dres = flux_z.lo.res.(fw) - res_ann;
    %         dres_lat = interp1(grid.dim3.lat, dres, lat_eval);
    %         res_ann_l = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.res.(fw) - res_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         res_ann_o = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.res.(fw) - res_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % DELTA RES lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta(\nabla\cdot F_m)$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], dres_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\nabla\cdot F_m)$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L+O}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta(\\nabla\\cdot F_m)$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dres_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf = flux_z.lo.stf.(fw);
    %         stf_lat = interp1(grid.dim3.lat, stf, lat_eval);
    %         comp1 = sftlf*1e-2.*flux_z.l.stf.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.stf.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % STF lat x mon dependence of RCE and STFE
    %         var_text = '$\mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], stf_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
    %         if any(strcmp(type, {'estf5', 'estfi'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$\mathrm{LH+SH}$', '$(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$(\mathrm{LH+SH})_{\mathrm{L}}$', '$(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_stf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf_ann = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat_eval);
    %         dstf = flux_z.lo.stf.(fw) - stf_ann;
    %         dstf_lat = interp1(grid.dim3.lat, dstf, lat_eval);
    %         stf_ann_l = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.stf.(fw) - stf_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);
    %         stf_ann_o = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.stf.(fw) - stf_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);
    %         % DELTA STF lat x mon dependence of RCE and STFE
    %         var_text = '$\Delta \mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], dstf_lat, 'k');
    %         c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
    %         c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
    %         if any(strcmp(type, {'estf5', 'estfi'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\mathrm{LH+SH})$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\delta \\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dstf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %     end
    % end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
                folder = sprintf('%s/dmse/%s/%s/0_lat_%g', par.plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat_eval);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat_eval);
                stf = flux_z.(land).stf.(fw);
                stf_lat = interp1(grid.dim3.lat, stf, lat_eval);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12], ra_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                stf=plot([1:12], stf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat_eval);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat_eval);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                res_ann_lat = interp1(grid.dim3.lat, res_ann, lat_eval);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat_eval);
                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat_eval);
                dstf = flux_z.(land).stf.(fw) - stf_ann;
                dstf_lat = interp1(grid.dim3.lat, dstf, lat_eval);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12], dra_lat, 'color', 0.5*[1 1 1]);
                dres=plot([1:12], dres_lat, 'color', par.maroon);
                dstf=plot([1:12], dstf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse', folder), '-dpng', '-r300');
                close;

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat_eval);
                sw = flux_z.(land).sw;
                sw_lat = interp1(grid.dim3.lat, sw, lat_eval);
                lw = flux_z.(land).lw;
                lw_lat = interp1(grid.dim3.lat, lw, lat_eval);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12], ra_lat, 'k');
                sw=plot([1:12], sw_lat, 'color', par.blue);
                lw=plot([1:12], lw_lat, 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra sw lw], '$R_a$', '$\mathrm{Net SW}$', '$\mathrm{Net LW}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_ra', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat_eval);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat_eval);
                sw_ann = repmat(nanmean(flux_z.(land).sw, 2), [1 12]);
                sw_ann_lat = interp1(grid.dim3.lat, sw_ann, lat_eval);
                dsw = flux_z.(land).sw - sw_ann;
                dsw_lat = interp1(grid.dim3.lat, dsw, lat_eval);
                lw_ann = repmat(nanmean(flux_z.(land).lw, 2), [1 12]);
                lw_ann_lat = interp1(grid.dim3.lat, lw_ann, lat_eval);
                dlw = flux_z.(land).lw - lw_ann;
                dlw_lat = interp1(grid.dim3.lat, dlw, lat_eval);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12], dra_lat, 'color', 0.5*[1 1 1]);
                dsw=plot([1:12], dsw_lat, 'color', par.blue);
                dlw=plot([1:12], dlw_lat, 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dsw dlw], '$\Delta R_a$', '$\Delta(\mathrm{Net SW})$', '$\Delta (\mathrm{Net LW})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dra', folder), '-dpng', '-r300');
                close;

                if strcmp(type, 'gcm')
                    stf = flux_z.(land).stf.(fw);
                    stf_lat = interp1(grid.dim3.lat, stf, lat_eval);
                    hfls = flux_z.(land).hfls;
                    hfls_lat = interp1(grid.dim3.lat, hfls, lat_eval);
                    hfss = flux_z.(land).hfss;
                    hfss_lat = interp1(grid.dim3.lat, hfss, lat_eval);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], stf_lat, 'k');
                    lh=plot([1:12], hfls_lat, 'color', par.blue);
                    sh=plot([1:12], hfss_lat, 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\mathrm{LH+SH}$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_stf', folder), '-dpng', '-r300');
                    close;

                    stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                    stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat_eval);
                    dstf = flux_z.(land).stf.(fw) - stf_ann;
                    dstf_lat = interp1(grid.dim3.lat, dstf, lat_eval);
                    hfls_ann = repmat(nanmean(flux_z.(land).hfls, 2), [1 12]);
                    hfls_ann_lat = interp1(grid.dim3.lat, hfls_ann, lat_eval);
                    dhfls = flux_z.(land).hfls - hfls_ann;
                    dhfls_lat = interp1(grid.dim3.lat, dhfls, lat_eval);
                    hfss_ann = repmat(nanmean(flux_z.(land).hfss, 2), [1 12]);
                    hfss_ann_lat = interp1(grid.dim3.lat, hfss_ann, lat_eval);
                    dhfss = flux_z.(land).hfss - hfss_ann;
                    dhfss_lat = interp1(grid.dim3.lat, dhfss, lat_eval);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], dstf_lat, 'k');
                    lh=plot([1:12], dhfls_lat, 'color', par.blue);
                    sh=plot([1:12], dhfss_lat, 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, land_text, lat_eval));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\Delta(\mathrm{LH+SH})$', '$\Delta\mathrm{LH}$', '$\Delta\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dstf', folder), '-dpng', '-r300');
                    close;

                end

            end

        end % for mse dse
    end % for land

end % for function
function plot_dmse_equator_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [5 10 20 30];

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    for f = f_vec; fw = f{1};
        for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
            dlat = 0.25; % step size for standard lat grid
            lat = -lat_bound:dlat:lat_bound;
            clat = cosd(lat); % cosine of latitude for cosine weighting
            clat_mon = repmat(clat', [1 12]);

            folder = sprintf('%s/dmse/%s/0_equatorward_of_lat_%g', par.plotdir, fw, lat_bound);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            ra = flux_z.lo.ra.(fw);
            ra_lat = interp1(grid.dim3.lat, ra, lat);
            ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
            comp1 = sftlf*1e-2.*flux_z.l.ra.(fw);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            comp2 = (1-sftlf*1e-2).*flux_z.o.ra.(fw);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % RA lat x mon dependence of RCE and RAE
            var_text = '$R_a$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], ra_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$R_a$', '$R_{a,\mathrm{\,L+O}}$', '$R_{a,\mathrm{\,L}}$', '$R_{a,\mathrm{\,O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$R_a$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_ra_lo_decomp', folder), '-dpng', '-r300');
            close;

            ra_ann = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
            ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
            ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
            dra = flux_z.lo.ra.(fw) - ra_ann;
            dra_lat = interp1(grid.dim3.lat, dra, lat);
            dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
            ra_ann_l = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.ra.(fw) - ra_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            ra_ann_o = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.ra.(fw) - ra_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % DELTA RA lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_a$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dra_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\Delta(R_a)$', '$\Delta(R_{a,\mathrm{\,L+O}})$', '$\Delta(R_{a,\mathrm{\,L}})$', '$\Delta(R_{a,\mathrm{\,O}})$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\Delta R_a$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dra_lo_decomp', folder), '-dpng', '-r300');
            close;

            res = flux_z.lo.res.(fw);
            res_lat = interp1(grid.dim3.lat, res, lat);
            res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
            comp1 = sftlf*1e-2.*flux_z.l.res.(fw);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            comp2 = (1-sftlf*1e-2).*flux_z.o.res.(fw);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % RES lat x mon dependence of RCE and RAE
            var_text = '$\nabla\cdot F_m$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], res_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\nabla\cdot F_m$', '$\nabla\cdot F_{m,\mathrm{\,L+O}}$', '$\nabla\cdot F_{m,\mathrm{\,L}}$', '$\nabla\cdot F_{m,\mathrm{\,O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\nabla\\cdot F_m$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_res_lo_decomp', folder), '-dpng', '-r300');
            close;

            res_ann = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
            res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
            res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
            dres = flux_z.lo.res.(fw) - res_ann;
            dres_lat = interp1(grid.dim3.lat, dres, lat);
            dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
            res_ann_l = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.res.(fw) - res_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            res_ann_o = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.res.(fw) - res_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % DELTA RES lat x mon dependence of RCE and RAE
            var_text = '$\Delta(\nabla\cdot F_m)$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dres_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\Delta(\nabla\cdot F_m)$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L+O}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,O}})$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\Delta(\\nabla\\cdot F_m)$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dres_lo_decomp', folder), '-dpng', '-r300');
            close;

            stf = flux_z.lo.stf.(fw);
            stf_lat = interp1(grid.dim3.lat, stf, lat);
            stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
            comp1 = sftlf*1e-2.*flux_z.l.stf.(fw);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            comp2 = (1-sftlf*1e-2).*flux_z.o.stf.(fw);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % STF lat x mon dependence of RCE and STFE
            var_text = '$\mathrm{LH+SH}$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], stf_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':k', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\mathrm{LH+SH}$', '$(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$(\mathrm{LH+SH})_{\mathrm{L}}$', '$(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\mathrm{LH+SH}$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_stf_lo_decomp', folder), '-dpng', '-r300');
            close;

            stf_ann = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
            stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
            stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
            dstf = flux_z.lo.stf.(fw) - stf_ann;
            dstf_lat = interp1(grid.dim3.lat, dstf, lat);
            dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
            stf_ann_l = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.stf.(fw) - stf_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
            stf_ann_o = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.stf.(fw) - stf_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
            % DELTA STF lat x mon dependence of RCE and STFE
            var_text = '$\Delta \mathrm{LH+SH}$';
            figure(); clf; hold all; box on;
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dstf_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\Delta(\mathrm{LH+SH})$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\delta \\mathrm{LH+SH}$ (Wm$^{-2}$)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dstf_lo_decomp', folder), '-dpng', '-r300');
            close;

        end
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                lat = -lat_bound:dlat:lat_bound;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);


                folder = sprintf('%s/dmse/%s/%s/0_equatorward_of_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat);
                ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                stf = flux_z.(land).stf.(fw);
                stf_lat = interp1(grid.dim3.lat, stf, lat);
                stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12], ra_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                stf=plot([1:12], stf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
                ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat);
                dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
                res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat);
                dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
                stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
                dstf = flux_z.(land).stf.(fw) - stf_ann;
                dstf_lat = interp1(grid.dim3.lat, dstf, lat);
                dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12], dra_lat, 'color', 0.5*[1 1 1]);
                dres=plot([1:12], dres_lat, 'color', par.maroon);
                dstf=plot([1:12], dstf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse', folder), '-dpng', '-r300');
                close;

                ra = flux_z.(land).ra.(fw);
                ra_lat = interp1(grid.dim3.lat, ra, lat);
                ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
                sw = flux_z.(land).sw;
                sw_lat = interp1(grid.dim3.lat, sw, lat);
                sw_lat = nansum(sw_lat.*clat_mon)/nansum(clat);
                lw = flux_z.(land).lw;
                lw_lat = interp1(grid.dim3.lat, lw, lat);
                lw_lat = nansum(lw_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12], ra_lat, 'k');
                sw=plot([1:12], sw_lat, 'color', par.blue);
                lw=plot([1:12], lw_lat, 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra sw lw], '$R_a$', '$\mathrm{Net SW}$', '$\mathrm{Net LW}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_ra', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw), 2), [1 12]);
                ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
                ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
                dra = flux_z.(land).ra.(fw) - ra_ann;
                dra_lat = interp1(grid.dim3.lat, dra, lat);
                dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
                sw_ann = repmat(nanmean(flux_z.(land).sw, 2), [1 12]);
                sw_ann_lat = interp1(grid.dim3.lat, sw_ann, lat);
                sw_ann_lat = nansum(sw_ann_lat.*clat_mon)/nansum(clat);
                dsw = flux_z.(land).sw - sw_ann;
                dsw_lat = interp1(grid.dim3.lat, dsw, lat);
                dsw_lat = nansum(dsw_lat.*clat_mon)/nansum(clat);
                lw_ann = repmat(nanmean(flux_z.(land).lw, 2), [1 12]);
                lw_ann_lat = interp1(grid.dim3.lat, lw_ann, lat);
                lw_ann_lat = nansum(lw_ann_lat.*clat_mon)/nansum(clat);
                dlw = flux_z.(land).lw - lw_ann;
                dlw_lat = interp1(grid.dim3.lat, dlw, lat);
                dlw_lat = nansum(dlw_lat.*clat_mon)/nansum(clat);
                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12], dra_lat, 'color', 0.5*[1 1 1]);
                dsw=plot([1:12], dsw_lat, 'color', par.blue);
                dlw=plot([1:12], dlw_lat, 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dsw dlw], '$\Delta R_a$', '$\Delta(\mathrm{Net SW})$', '$\Delta (\mathrm{Net LW})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dra', folder), '-dpng', '-r300');
                close;

                if strcmp(type, 'gcm')
                    stf = flux_z.(land).stf.(fw);
                    stf_lat = interp1(grid.dim3.lat, stf, lat);
                    stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
                    hfls = flux_z.(land).hfls;
                    hfls_lat = interp1(grid.dim3.lat, hfls, lat);
                    hfls_lat = nansum(hfls_lat.*clat_mon)/nansum(clat);
                    hfss = flux_z.(land).hfss;
                    hfss_lat = interp1(grid.dim3.lat, hfss, lat);
                    hfss_lat = nansum(hfss_lat.*clat_mon)/nansum(clat);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], stf_lat, 'k');
                    lh=plot([1:12], hfls_lat, 'color', par.blue);
                    sh=plot([1:12], hfss_lat, 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\mathrm{LH+SH}$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_stf', folder), '-dpng', '-r300');
                    close;

                    stf_ann = repmat(nanmean(flux_z.(land).stf.(fw), 2), [1 12]);
                    stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
                    stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
                    dstf = flux_z.(land).stf.(fw) - stf_ann;
                    dstf_lat = interp1(grid.dim3.lat, dstf, lat);
                    dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
                    hfls_ann = repmat(nanmean(flux_z.(land).hfls, 2), [1 12]);
                    hfls_ann_lat = interp1(grid.dim3.lat, hfls_ann, lat);
                    hfls_ann_lat = nansum(hfls_ann_lat.*clat_mon)/nansum(clat);
                    dhfls = flux_z.(land).hfls - hfls_ann;
                    dhfls_lat = interp1(grid.dim3.lat, dhfls, lat);
                    dhfls_lat = nansum(dhfls_lat.*clat_mon)/nansum(clat);
                    hfss_ann = repmat(nanmean(flux_z.(land).hfss, 2), [1 12]);
                    hfss_ann_lat = interp1(grid.dim3.lat, hfss_ann, lat);
                    hfss_ann_lat = nansum(hfss_ann_lat.*clat_mon)/nansum(clat);
                    dhfss = flux_z.(land).hfss - hfss_ann;
                    dhfss_lat = interp1(grid.dim3.lat, dhfss, lat);
                    dhfss_lat = nansum(dhfss_lat.*clat_mon)/nansum(clat);
                    % ALL lat x mon dependence of RCE and RAE
                    figure(); clf; hold all; box on;
                    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                    stf=plot([1:12], dstf_lat, 'k');
                    lh=plot([1:12], dhfls_lat, 'color', par.blue);
                    sh=plot([1:12], dhfss_lat, 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound, lat_bound));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\Delta(\mathrm{LH+SH})$', '$\Delta\mathrm{LH}$', '$\Delta\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dstf', folder), '-dpng', '-r300');
                    close;

                end

            end

        end % for mse dse
    end % for land

end % for function
function plot_dmse_toasfc_midlatitude_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [5 10 15 20 -5 -10 -15 -20];
    lat_bound_list = [5 -5];

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center;
                else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dmse_toasfc/%s/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                rtoa = flux_z.(land).rtoa;
                rtoa_lat = interp1(grid.dim3.lat, rtoa, lat);
                rtoa_lat = nansum(rtoa_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                sfc = flux_z.(land).sfc.(fw);
                sfc_lat = interp1(grid.dim3.lat, sfc, lat);
                sfc_lat = nansum(sfc_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                sfc=plot([1:12], sfc_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([rtoa res sfc], '$F_{\mathrm{TOA}}$', '$\nabla\cdot F_m$', '$F_\mathrm{SFC}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_toa_sfc', folder), '-dpng', '-r300');
                close;

                % NO LEGEND ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'k');
                res=plot([1:12], res_lat, 'color', par.maroon);
                sfc=plot([1:12], sfc_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_toa_sfc_noleg', folder), '-dpng', '-r300');
                close;

                sw = flux_z.(land).sw;
                sw_lat = interp1(grid.dim3.lat, sw, lat);
                sw_lat = nansum(sw_lat.*clat_mon)/nansum(clat);
                olr = flux_z.(land).olr;
                olr_lat = interp1(grid.dim3.lat, olr, lat);
                olr_lat = nansum(olr_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                shf = flux_z.(land).shf.(fw);
                shf_lat = interp1(grid.dim3.lat, shf, lat);
                shf_lat = nansum(shf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                sw=plot([1:12], sw_lat, 'color', par.yellow);
                olr=plot([1:12], olr_lat, 'color', par.green);
                res=plot([1:12], res_lat, 'color', par.maroon);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([sw olr res shf], '$\mathrm{SWABS}$', '$\mathrm{OLR}$', '$\nabla\cdot F_m$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swabs_olr_shf', folder), '-dpng', '-r300');
                close;

                sw_ann = repmat(nanmean(flux_z.(land).sw, 2), [1 12]);
                olr_ann = repmat(nanmean(flux_z.(land).olr, 2), [1 12]);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                shf_ann = repmat(nanmean(flux_z.(land).shf.(fw), 2), [1 12]);

                dsw = flux_z.(land).sw - sw_ann;
                dsw_lat = interp1(grid.dim3.lat, dsw, lat);
                dsw_lat = nansum(dsw_lat.*clat_mon)/nansum(clat);
                dolr = flux_z.(land).olr - olr_ann;
                dolr_lat = interp1(grid.dim3.lat, dolr, lat);
                dolr_lat = nansum(dolr_lat.*clat_mon)/nansum(clat);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat);
                dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
                dshf = flux_z.(land).shf.(fw) - shf_ann;
                dshf_lat = interp1(grid.dim3.lat, dshf, lat);
                dshf_lat = nansum(dshf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dsw=plot([1:12],  dsw_lat, 'color', par.yellow);
                dolr=plot([1:12], dolr_lat, 'color', par.green);
                dres=plot([1:12], dres_lat, 'color', par.maroon);
                dshf=plot([1:12], dshf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dsw dolr dres dshf], '$\mathrm{\Delta SWABS}$', '$\mathrm{\Delta OLR}$', '$\Delta(\nabla\cdot F_m)$', '$\Delta F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dswabs_dolr_dshf', folder), '-dpng', '-r300');
                close;

                rtoa = flux_z.(land).rtoa;
                rtoa_lat = interp1(grid.dim3.lat, rtoa, lat);
                rtoa_lat = nansum(rtoa_lat.*clat_mon)/nansum(clat);
                swsfc = flux_z.(land).swsfc;
                swsfc_lat = interp1(grid.dim3.lat, swsfc, lat);
                swsfc_lat = nansum(swsfc_lat.*clat_mon)/nansum(clat);
                res = flux_z.(land).res.(fw);
                res_lat = interp1(grid.dim3.lat, res, lat);
                res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
                shf = flux_z.(land).shf.(fw);
                shf_lat = interp1(grid.dim3.lat, shf, lat);
                shf_lat = nansum(shf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                rtoa=plot([1:12], rtoa_lat, 'color', 'k');
                swsfc=plot([1:12], swsfc_lat, 'color', par.yellow);
                res=plot([1:12], res_lat, 'color', par.maroon);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([rtoa swsfc res shf], '$F_{\mathrm{TOA}}$', '$F_\mathrm{SW,\,SFC}$', '$\nabla\cdot F_m$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_rtoa_swsfc_shf', folder), '-dpng', '-r300');
                close;

                % just SWSFC and SHF lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', par.yellow);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([swsfc shf], '$F_\mathrm{SW,\,SFC}$', '$F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_shf', folder), '-dpng', '-r300');
                close;

                % NOLEGEND just SWSFC and SHF lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', par.yellow);
                shf=plot([1:12], shf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_shf_noleg', folder), '-dpng', '-r300');
                close;

                rtoa_ann = repmat(nanmean(flux_z.(land).rtoa, 2), [1 12]);
                swsfc_ann = repmat(nanmean(flux_z.(land).swsfc, 2), [1 12]);
                res_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);
                shf_ann = repmat(nanmean(flux_z.(land).shf.(fw), 2), [1 12]);

                drtoa = flux_z.(land).rtoa - rtoa_ann;
                drtoa_lat = interp1(grid.dim3.lat, drtoa, lat);
                drtoa_lat = nansum(drtoa_lat.*clat_mon)/nansum(clat);
                dswsfc = flux_z.(land).swsfc - swsfc_ann;
                dswsfc_lat = interp1(grid.dim3.lat, dswsfc, lat);
                dswsfc_lat = nansum(dswsfc_lat.*clat_mon)/nansum(clat);
                dres = flux_z.(land).res.(fw) - res_ann;
                dres_lat = interp1(grid.dim3.lat, dres, lat);
                dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
                dshf = flux_z.(land).shf.(fw) - shf_ann;
                dshf_lat = interp1(grid.dim3.lat, dshf, lat);
                dshf_lat = nansum(dshf_lat.*clat_mon)/nansum(clat);
                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                drtoa=plot([1:12],  drtoa_lat, 'color', 'k');
                dswsfc=plot([1:12], dswsfc_lat, 'color', par.yellow);
                dres=plot([1:12], dres_lat, 'color', par.maroon);
                dshf=plot([1:12], dshf_lat, 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([drtoa dswsfc dres dshf], '$\mathrm{\Delta F_\mathrm{TOA}}$', '$\mathrm{\Delta F_\mathrm{SW,\,SFC}}$', '$\Delta(\nabla\cdot F_m)$', '$\Delta F_\mathrm{SHF}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_drtoa_dswsfc_dshf', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
function plot_trop(type, par)
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
    end

    load(sprintf('%s/ztrop.mat', prefix)); % read tropopause in z
    load(sprintf('%s/ztrop_z.mat', prefix)); ztrop_zon = ztrop_z; clear ztrop_z; % read tropopause of latitudinally averaged temperature in z
    load(sprintf('%s/ptrop.mat', prefix)); % read tropopause in pressure
    load(sprintf('%s/ptrop_z.mat', prefix)); ptrop_zon = ptrop_z; clear ptrop_z; % read tropopause of latitudinally averaged temperature in pressure

    % zonal mean
    ztrop_z = squeeze(nanmean(ztrop,1));
    ptrop_z = squeeze(nanmean(ptrop,1));
    % annual mean
    ztrop_zt = squeeze(nanmean(ztrop_z,2));
    ztrop_zont = squeeze(nanmean(ztrop_zon,2));
    ptrop_zt = squeeze(nanmean(ptrop_z,2));
    ptrop_zont = squeeze(nanmean(ptrop_zon,2));

    % z tropopause, lat x height
    figure(); clf; hold all; box on;
    plot(grid.dim3.lat, ztrop_zt*1e-3, 'k');
    xlabel('latitude (deg)'); ylabel('z (km)');
    title(sprintf('Tropopause Height'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0 20], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/trop/ztrop_lat', par.plotdir), '-dpng', '-r300');
    close;

    % z tropopause zonally averaged first, lat x height
    figure(); clf; hold all; box on;
    plot(grid.dim3.lat, ztrop_zont*1e-3, 'k');
    xlabel('latitude (deg)'); ylabel('z (km)');
    title(sprintf('Tropopause Height'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0 20], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/trop/ztrop_zon_lat', par.plotdir), '-dpng', '-r300');
    close;

    % p tropopause, lat x height
    figure(); clf; hold all; box on;
    plot(grid.dim3.lat, ptrop_zt*1e-2, 'k');
    xlabel('latitude (deg)'); ylabel('p (hPa)');
    title(sprintf('Tropopause Height'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ydir', 'reverse', 'yscale', 'log', 'ylim', [50 1000], 'ytick', [50 100 200 300 400:200:1000], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/trop/ptrop_lat', par.plotdir), '-dpng', '-r300');
    close;

    % p tropopause zonally averaged first, lat x height
    figure(); clf; hold all; box on;
    plot(grid.dim3.lat, ptrop_zont*1e-2, 'k');
    xlabel('latitude (deg)'); ylabel('p (hPa)');
    title(sprintf('Tropopause Height'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ydir', 'reverse', 'yscale', 'log', 'ylim', [50 1000], 'ytick', [50 100 200 300 400:200:1000], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/trop/ptrop_zon_lat', par.plotdir), '-dpng', '-r300');
    close;

end

function plot_rad_lat(par)
    load('/project2/tas1/miyawaki/projects/002/data/proc/comp/comp_zt.mat');

    for fn = {'tsr', 'ttr', 'net', 'ssr', 'swabs', 'str', 'ra', 'surface_turbulent_plus_LW'}; fname = fn{1};
        if strcmp(fname, 'tsr'); ytext = 'Net SW TOA';
        elseif strcmp(fname, 'ttr'); ytext = 'OLR';
        elseif strcmp(fname, 'net'); ytext = 'Net TOA';
        elseif strcmp(fname, 'ssr'); ytext = 'Net SW SFC';
        elseif strcmp(fname, 'str'); ytext = 'Net LW SFC';
        elseif strcmp(fname, 'swabs'); ytext = 'SWABS';
        elseif strcmp(fname, 'ra'); ytext = '$R_a$';
        elseif strcmp(fname, 'surface_turbulent_plus_LW'); ytext = 'LH + SH + Net LW SFC';
        end

        figure(); clf; hold all;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        h.erai = plot(lat, erai_zt.rad.(fname), 'color', par.yellow);
        h.era5 = plot(lat, era5_zt.rad.(fname), 'color', par.orange);
        if ~strcmp(fname, 'surface_turbulent_plus_LW'); h.ceres = plot(lat, ceres_zt.rad.(fname), 'color', par.green); end;
        if ~any(strcmp(fname, {'str', 'ra'})); h.don = plot(lat, don_zt.(fname), 'color', par.blue); end;
        xlabel('latitude (deg)'); ylabel(sprintf('%s (Wm$^{-2}$)', ytext));
        if any(strcmp(fname, {'str', 'ra'})); legend([h.erai h.era5 h.ceres], 'ERA-I', 'ERA5', 'CERES4.1', 'location', 'eastoutside');
        elseif strcmp(fname, 'surface_turbulent_plus_LW'); legend([h.erai h.era5 h.don], 'ERA-I', 'ERA5', 'DB13', 'location', 'eastoutside');
        else legend([h.erai h.era5 h.ceres h.don], 'ERA-I', 'ERA5', 'CERES4.1', 'DB13', 'location', 'eastoutside'); end;
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s', fname), '-dpng', '-r300');
        close;

        if ~any(strcmp(fname, {'str', 'ra'}));
            figure(); clf; hold all;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            h.erai = plot(lat, erai_zt.rad.(fname) - don_zt.(fname), 'color', par.yellow);
            h.era5 = plot(lat, era5_zt.rad.(fname) - don_zt.(fname), 'color', par.orange);
            if ~strcmp(fname, 'surface_turbulent_plus_LW'); h.ceres = plot(lat, ceres_zt.rad.(fname) - don_zt.(fname), 'color', par.green); end;
            xlabel('latitude (deg)'); ylabel(sprintf('$\\Delta$ (%s) (Wm$^{-2}$)', ytext));
            if ~strcmp(fname, 'surface_turbulent_plus_LW'); legend([h.erai h.era5 h.ceres], 'ERA-I $-$ DB13', 'ERA5 $-$ DB13', 'CERES4.1 $-$ DB13', 'location', 'eastoutside');
            else; legend([h.erai h.era5], 'ERA-I $-$ DB13', 'ERA5 $-$ DB13', 'location', 'eastoutside'); end;
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s-diff-db13', fname), '-dpng', '-r300');
            close;
        else
            figure(); clf; hold all;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            h.erai = plot(lat, erai_zt.rad.(fname) - ceres_zt.rad.(fname), 'color', par.yellow);
            h.era5 = plot(lat, era5_zt.rad.(fname) - ceres_zt.rad.(fname), 'color', par.orange);
            xlabel('latitude (deg)'); ylabel(sprintf('$\\Delta$ (%s) (Wm$^{-2}$)', ytext));
            legend([h.erai h.era5], 'ERA-I $-$ CERES4.1', 'ERA5 $-$ CERES4.1', 'location', 'eastoutside')
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/%s-diff-ceres', fname), '-dpng', '-r300');
            close;
        end

    end

end
function plot_rad_lon_lat(par)
    load('/project2/tas1/miyawaki/projects/002/data/proc/comp/ceres_2001_2009.mat');
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    [mesh_lat, mesh_lon] = meshgrid(grid.ceres.dim2.lon, lat);

    figure(); clf; hold all;
    cmp = colCog(48);
    colormap(cmp);
    contourf(mesh_lat, mesh_lon, ceres_t.rad.ra', [-360 -240:60:-120 -90:30:-30 -15], 'linecolor', 'none');
    contour(par.landlon+180, par.landlat, par.land, [1 1], 'linecolor', 0.3*[1 1 1], 'linewidth', 0.5);
    caxis([-360 360]);
    xlabel('longitude (deg)'); ylabel('latitude (deg)');
    title('$R_a$, CERES4.1 2001--2009');
    cb = colorbar('limits', [-360 -15], 'ytick', [-360 -240:60:-120 -90:30:-30 -15], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('$R_a$ (Wm$^{-2}$)'));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [0:60:360], 'xminortick', 'on', 'ytick', [-90:30:90], 'yminortick', 'on')
    print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lon_lat/ra_ceres'), '-dpng', '-r300');
    close;

end
function plot_tediv_lat(par)
    par.model = 'MPI-ESM-LR';
    gcm=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/gcm/%s/%s/std/flux_zt.mat', par.model, par.gcm.clim));
    don=load('/project2/tas1/miyawaki/projects/002/data/raw/don/radiation_dynamics_climatology.mat');
    don_zt.TEDIV = squeeze(nanmean(nanmean(don.TEDIV, 1), 3));
    don_zt.TEDIV = interp1(don.lat, don_zt.TEDIV, gcm.lat);

    figure(); clf; hold all;
    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
    h_gcm = plot(gcm.lat, gcm.flux_zt.lo.ann.res.mse, 'color', par.maroon);
    h_don = plot(gcm.lat, don_zt.TEDIV, 'color', par.blue);
    xlabel('latitude (deg)'); ylabel(sprintf('$\\nabla \\cdot F_m$ (Wm$^{-2}$)'));
    legend([h_gcm h_don], par.model, 'DB13', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('/project2/tas1/miyawaki/projects/002/figures/comp/lat/tediv'), '-dpng', '-r300');
    close;

end

function plot_rcae(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end
                % lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_z.(land).(fw).(crit));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, crit_text, land_text)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_mon_lat', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_lat_lon', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_rcae_alt_overlay(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE_ALT data
    load(sprintf('%s/%s/ta_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ma_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_malr_diff_orig = ga_malr_diff; clear ga_diff;
    load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_dalr_bl_diff_orig = ga_dalr_bl_diff; clear ga_bl_diff;
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    % load(sprintf('%s/albedo.mat', prefix));
    % load(sprintf('%s/sn.mat', prefix));
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % albedo0 = squeeze(nanmean(albedo,1)); % zonal average
    % albedo = nan([size(albedo0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with igagesc)
    % albedo(:,2:13) = albedo0;
    % albedo(:,1) = 1/2*(albedo0(:,1)+albedo0(:,12));
    % albedo(:,14) = 1/2*(albedo0(:,1)+albedo0(:,12));
    % albedo = interp1(grid.dim3.lat, albedo, lat); % interpolate to standard lat

    % sn0 = squeeze(nanmean(sn,1)); % zonal average
    % sn = nan([size(sn0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with igagesc)
    % sn(:,2:13) = sn0;
    % sn(:,1) = 1/2*(sn0(:,1)+sn0(:,12));
    % sn(:,14) = 1/2*(sn0(:,1)+sn0(:,12));
    % sn = interp1(grid.dim3.lat, sn, lat); % interpolate to standard lat

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        ga_malr_diff0 = ga_malr_diff_orig.(land);
        ga_malr_diff = nan([size(ga_malr_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_malr_diff(:,2:13) = ga_malr_diff0;
        ga_malr_diff(:,1) = 1/2*(ga_malr_diff0(:,1)+ga_malr_diff0(:,12));
        ga_malr_diff(:,14) = 1/2*(ga_malr_diff0(:,1)+ga_malr_diff0(:,12));

        ga_dalr_bl_diff0 = ga_dalr_bl_diff_orig.(land);
        ga_dalr_bl_diff = nan([size(ga_dalr_bl_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_dalr_bl_diff(:,2:13) = ga_dalr_bl_diff0;
        ga_dalr_bl_diff(:,1) = 1/2*(ga_dalr_bl_diff0(:,1)+ga_dalr_bl_diff0(:,12));
        ga_dalr_bl_diff(:,14) = 1/2*(ga_dalr_bl_diff0(:,1)+ga_dalr_bl_diff0(:,12));

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_alt_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end
                [mesh_lat, mesh_mon] = meshgrid([0.5 1:12 12.5], lat);

                % GA MALR contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, ga_malr_diff, par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$', par.model, crit_text, land_text, par.ga_thresh));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh)); end
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_malr_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA_BL DALR contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), crit_text, land_text, par.ga_bl_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, crit_text, land_text, par.ga_bl_thresh));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), crit_text, land_text, par.ga_bl_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_dalr_bl_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA MALR and GA_BL DALR contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, ga_malr_diff, par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_malr_ga_dalr_bl_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_alt_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text));
                        elseif any(strcmp(type, {'echam'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text)); end
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_alt_lat_lon_overlay', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_rcae_alt_rc_overlay(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_rc_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_rc_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE_ALT data
    load(sprintf('%s/%s/inv_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_diff_mon_lat.mat', prefix_proc, par.lat_interp)); ga_diff_orig = ga_diff; clear ga_diff;
    load(sprintf('%s/dtdzz.mat', prefix));
    load(sprintf('%s/dtmdzz.mat', prefix));
    load(sprintf('%s/albedo.mat', prefix));
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % ga_diff0 = (dtmdzz - dtdzz)./dtmdzz * 1e2; % percentage difference of moist adiabatic vs GCM lapse rates
    % clear dtmdzz dtdzz;
    % ga_diff0 = squeeze(nanmean(ga_diff0,1)); % zonal average
    % ga_diff = nan([size(ga_diff0,1) size(ga_diff0,2) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with igagesc)
    % ga_diff(:,:,2:13) = ga_diff0;
    % ga_diff(:,:,1) = 1/2*(ga_diff0(:,:,1)+ga_diff0(:,:,12));
    % ga_diff(:,:,14) = 1/2*(ga_diff0(:,:,1)+ga_diff0(:,:,12));
    % ga_diff = interp1(grid.dim3.lat, ga_diff, lat); % interpolate to standard lat
    % ga_diff = permute(ga_diff, [2 1 3]); % bring height front
    % ga_diff = interp1(par.pa, ga_diff, 1e2*linspace(1000,200,100)); % prepare to average between 1000-200 hPa
    % ga_diff = squeeze(nanmean(ga_diff,1)); % take vertical average

    albedo0 = squeeze(nanmean(albedo,1)); % zonal average
    albedo = nan([size(albedo0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with igagesc)
    albedo(:,2:13) = albedo0;
    albedo(:,1) = 1/2*(albedo0(:,1)+albedo0(:,12));
    albedo(:,14) = 1/2*(albedo0(:,1)+albedo0(:,12));
    albedo = interp1(grid.dim3.lat, albedo, lat); % interpolate to standard lat

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_alt_rc_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end

                ga_diff0 = ga_diff_orig.(land);
                ga_diff = nan([size(ga_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
                ga_diff(:,2:13) = ga_diff0;
                ga_diff(:,1) = 1/2*(ga_diff0(:,1)+ga_diff0(:,12));
                ga_diff(:,14) = 1/2*(ga_diff0(:,1)+ga_diff0(:,12));

                inv_str0 = squeeze(inv.(land)(:,:,2));
                inv_str = nan([size(inv_str0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
                inv_str(:,2:13) = inv_str0;
                inv_str(:,1) = 1/2*(inv_str0(:,1)+inv_str0(:,12));
                inv_str(:,14) = 1/2*(inv_str0(:,1)+inv_str0(:,12));

                [mesh_lat, mesh_mon] = meshgrid([0.5 1:12 12.5], lat);

                % GA contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_rc_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $(\\Gamma_m - \\Gamma )/ \\Gamma_m < %g \\%%$', par.model, crit_text, land_text, par.ga_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_rc_mon_lat_ga_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % ALBEDO contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_rc_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, albedo, par.albedo_thresh*[1 1], 'k', 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\alpha=%g$', par.model, crit_text, land_text, par.albedo_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_rc_mon_lat_albedo_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_alt_rc_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_alt_rc_lat_lon_overlay', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_rcae_rc(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_rc_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_rc_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_rc_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end
                % lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_rc_z.(land).(fw).(crit));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, crit_text, land_text)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_rc_mon_lat', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_rc_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_rc_lat_lon', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_rcae_alt(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE_ALT data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_alt_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end
                % lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, crit_text, land_text)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_alt_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_alt_lat_lon', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function
function plot_rcae_alt_rc(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_rc_z.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/eps_%g_ga_%g/rcae_alt_rc_t.mat', prefix_proc, par.lat_interp, par.ep, par.ga)); % load lat x lon RCAE_ALT data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for fn = fieldnames(rcae_alt_rc_z.(land).(fw))'; crit = fn{1};
                if strcmp(crit, 'def'); crit_text = 'No flags';
                elseif strcmp(crit, 'jak');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 50$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 50$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak30');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 30$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 30$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'jak10');
                    if strcmp(fw, 'dse'); crit_text = '$|\nabla \cdot F_s| < 10$ Wm$^{-2}$';
                    else crit_text = '$|\nabla \cdot F_m| < 10$ Wm$^{-2}$'; end;
                elseif strcmp(crit, 'pe'); crit_text = '$P-E>0$ flag';
                elseif strcmp(crit, 'cp'); crit_text = '$P_{ls}/P_{c}<0.5$ flag';
                elseif strcmp(crit, 'w500'); crit_text = '$\omega500<0$ flag';
                elseif strcmp(crit, 'vas2'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/2$ flag';
                elseif strcmp(crit, 'vas4'); crit_text = '$|v_{\,\mathrm{2\,m}}| < \max(|v_{\,\mathrm{2\,m}}|)/4$ flag';
                elseif strcmp(crit, 'vh2');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/2$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/2$ flag'; end;
                elseif strcmp(crit, 'vh3');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/3$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/3$ flag'; end;
                elseif strcmp(crit, 'vh4');
                    if strcmp(fw, 'dse'); crit_text = '$|F_s| < \max(|F_s|)/4$ flag';
                    else crit_text = '$|F_m| < \max(|F_m|)/4$ flag'; end;
                end
                % lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_rc_z.(land).(fw).(crit));
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, crit_text, land_text)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_rc_mon_lat', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                    % lat x lon of RCE and RAE
                    if ~any(strcmp(crit, {'vh2', 'vh3', 'vh4'}))
                        figure(); clf; hold all;
                        cmp = colCog(10);
                        colormap(cmp);
                        imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_alt_rc_t.(land).(time).(fw).(crit)');
                        caxis([-2 2]);
                        if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), crit_text, upper(time), land_text));
                        elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, crit_text, upper(time), land_text)); end;
                        xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                        set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                        print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/rcae_alt_rc_lat_lon', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                        close;
                    end % if vh2 vh3 vh4
                end % for time
            end % for crit
        end % for mse dse
    end % for land
end % function

function plot_dr1_equator_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [5 10 20 30];

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    for f = f_vec; fw = f{1};
        for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
            dlat = 0.25; % step size for standard lat grid
            lat = -lat_bound:dlat:lat_bound;
            clat = cosd(lat); % cosine of latitude for cosine weighting
            clat_mon = repmat(clat', [1 12]);

            folder = sprintf('%s/dr1/%s/0_equatorward_of_lat_%g', par.plotdir, fw, lat_bound);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
            r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

            dr1 = flux_z.lo.r1.(fw) - r1_ann;
            dr1_lat = interp1(grid.dim3.lat, dr1, lat);
            dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

            r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

            r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
            ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
            figure(); clf; hold all; box on;
                rcemax = par.ep-r1_ann_lat(1);
                vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                raemin = par.ga-r1_ann_lat(1);
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dr1_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound, lat_bound));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, -lat_bound, lat_bound)); end;
            legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
            xlabel('Month');
            ylabel(sprintf('$\\Delta R_1$ (unitless)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
            close;
        end
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                lat = -lat_bound:dlat:lat_bound;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_equatorward_of_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1 = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1_lat = interp1(grid.dim3.lat, comp1, lat);
                comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

                comp2_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2 = fm_ann./(ra_ann).^2.*delta_stf;
                comp2_lat = interp1(grid.dim3.lat, comp2, lat);
                comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound, lat_bound)); end;
                xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                c12=plot([1:12], comp1_lat+comp2_lat, '-.k');
                c1=plot([1:12], comp1_lat, '--k');
                c2=plot([1:12], comp2_lat, ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound, lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound, lat_bound)); end;
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta (LH + SH)$', 'location', 'eastoutside');
                xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land
end % for function
function plot_dr1_polar_line_old(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'merra2')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    sftlf = nanmean(sftlf, 1); % zonal average
    sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [-80 80];

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    for f = f_vec; fw = f{1};
        for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
            dlat = 0.25; % step size for standard lat grid
            if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; monlabel=par.monlabel; shiftby=0;
            else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; monlabel=par.monlabelsh; shiftby=6; end;
            clat = cosd(lat); % cosine of latitude for cosine weighting
            clat_mon = repmat(clat', [1 12]);

            folder = sprintf('%s/dr1/%s/0_poleward_of_lat_%g', par.plotdir, fw, lat_bound);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
            r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

            dr1 = flux_z.lo.r1.(fw) - r1_ann;
            dr1_lat = interp1(grid.dim3.lat, dr1, lat);
            dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

            r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat);
            comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

            r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat);
            comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
            ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
            figure(); clf; hold all; box on;
                rcemax = par.ep-r1_ann_lat(1);
                vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                raemin = par.ga-r1_ann_lat(1);
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
            c12=plot([1:12], circshift(comp1_lat+comp2_lat, shiftby, 2), '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], circshift(comp1_lat, shiftby, 2), '--', 'color', par.maroon);
            c2=plot([1:12], circshift(comp2_lat, shiftby, 2), ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
            elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
            legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
            % xlabel('Month');
            ylabel(sprintf('$\\Delta R_1$ (unitless)'));
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
            close;
        end
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'L+O';
        elseif strcmp(land, 'l'); land_text = 'L';
        elseif strcmp(land, 'o'); land_text = 'O';
        end
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
        elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; monlabel=par.monlabel; shiftby=0;
                else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; monlabel=par.monlabelsh; shiftby=6; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_poleward_of_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                r1_lat = nansum(r1_lat.*clat_mon)/nansum(clat);

                % R1 lat x mon dependence of RCE and RAE
                var_text = '$R_1$';
                figure(); clf; hold all; box on;
                ylim_lo = -0.5;
                ylim_up = 1.5;
                    rcemax = par.ep;
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    raemin = par.ga;
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([1 12], [1 1], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_r1', folder), '-dpng', '-r300');
                close;

                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1a = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1a_lat = interp1(grid.dim3.lat, comp1a, lat);
                comp1a_lat = nansum(comp1a_lat.*clat_mon)/nansum(clat);

                comp2a_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2a = fm_ann./(ra_ann).^2.*delta_stf;
                comp2a_lat = interp1(grid.dim3.lat, comp2a, lat);
                comp2a_lat = nansum(comp2a_lat.*clat_mon)/nansum(clat);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k','k'})
                yyaxis left
                ylim_lo = min(r1_lat); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max(r1_lat); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                r1=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k','k'})
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                r1=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                c12=plot([1:12], circshift(comp1a_lat+comp2a_lat, shiftby, 2), '-.k');
                c1=plot([1:12],  circshift(comp1a_lat, shiftby, 2), '--k');
                c2=plot([1:12],  circshift(comp2a_lat, shiftby, 2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta (LH + SH)$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_alt', folder), '-dpng', '-r300');
                close;

                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                comp1s = delta_fm./ra_ann;
                comp1s_lat = interp1(grid.dim3.lat, comp1s, lat);
                comp1s_lat = nansum(comp1s_lat.*clat_mon)/nansum(clat);

                delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                comp2s = -fm_ann./(ra_ann).^2.*delta_ra;
                comp2s_lat = interp1(grid.dim3.lat, comp2s, lat);
                comp2s_lat = nansum(comp2s_lat.*clat_mon)/nansum(clat);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k','k'})
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                r1=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                c12=plot([1:12], circshift(comp1s_lat+comp2s_lat, shiftby, 2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat, shiftby, 2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat, shiftby, 2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp', folder), '-dpng', '-r300');
                close;

                % NOLEG DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k','k'})
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                r1=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                rcemax = par.ep-r1_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                c12=plot([1:12], circshift(comp1s_lat+comp2s_lat, shiftby, 2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat, shiftby, 2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat, shiftby, 2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_noleg', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land
end % for function
function plot_dr2(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            % DELTA R2 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_2$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], flux_z.(land).r2.(fw) - repmat(nanmean(flux_z.(land).r2.(fw), 2), [1 12]));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_2$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr2/%s/%s/0_dr2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % FRACTIONAL DELTA R2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta R_2}{R_2}$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], (flux_z.(land).r2.(fw) - repmat(nanmean(flux_z.(land).r2.(fw), 2), [1 12]))./flux_z.(land).r2.(fw));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta R_2}{R_2}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr2/%s/%s/0_dr2_frac_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R2 DEF COMP1 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\mathrm{SH+LH})}{R_a}$';
            delta_fm = flux_z.(land).stf.(fw) - repmat(nanmean(flux_z.(land).stf.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp1 = delta_fm./ann_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\mathrm{SH+LH})}{R_a}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr2/%s/%s/0_dr2c1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R2 DEF COMP2 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\mathrm{SH+LH}}{R_a^2}\Delta R_a$';
            delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            ann_fm = repmat(nanmean(flux_z.(land).stf.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            comp2 = -ann_fm./(ann_ra).^2.*delta_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr2/%s/%s/0_dr2c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\mathrm{SH+LH})}{R_a}-\frac{\mathrm{SH+LH}}{R_a^2}\Delta R_a$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1+comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\mathrm{SH+LH})}{R_a}-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr2/%s/%s/0_dr1_c1+c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

        end % for time
    end % for mse dse
end % for land
function plot_dr1_decomp4(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            % DELTA R1 R3 R4 COMP1 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\nabla\cdot F_{TOA})}{R_a}$';
            delta_ftoa = flux_z.(land).ftoa.(fw) - repmat(nanmean(flux_z.(land).ftoa.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp1 = delta_ftoa./ann_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_{TOA})}{R_a}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\nabla\cdot F_{SFC})}{R_a}$';
            delta_fsfc = flux_z.(land).fsfc.(fw) - repmat(nanmean(flux_z.(land).fsfc.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp2 = delta_fsfc./ann_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_{SFC})}{R_a}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 COMP3 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\nabla\cdot F_{TOA}}{R_a^2}\Delta R_a$';
            delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            ann_ftoa = repmat(nanmean(flux_z.(land).ftoa.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            comp3 = -ann_ftoa./(ann_ra).^2.*delta_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp3);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\nabla\\cdot F_{TOA}}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c3_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 COMP4 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\nabla\cdot F_{SFC}}{R_a^2}\Delta R_a$';
            delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            ann_fsfc = repmat(nanmean(flux_z.(land).fsfc.(fw),2),[1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
            comp4 = -ann_fsfc./(ann_ra).^2.*delta_ra;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp4);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\nabla\\cdot F_{SFC}}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c4_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 SUM OF COMP1 and COMP3 (R3) lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_3$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1+comp3);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_3$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c1and3_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 SUM OF COMP2 and COMP4 (R4) lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_4$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp2+comp4);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_4$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c2and4_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 R3 R4 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_3+\Delta R_4$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1+comp2+comp3+comp4);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_1+\\Delta R_3$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1r3r4c1to4_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            end % for time
        end % for mse dse
    end % for land
function plot_dmse_line_comp(type, par)
    make_dirs(type, par)

    par.plotdir = sprintf('./figures/%s/%s', 'erai', par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'erai');
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', 'erai');
    era_grid = load(sprintf('%s/grid.mat', prefix)); % read grid data
    era_flux_z = load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    era_flux_t = load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    gcm_grid = load(sprintf('%s/grid.mat', prefix)); % read grid data
    gcm_flux_z = load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    gcm_flux_t = load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data

    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    lat_eval_list = [-85 -45 0 45 85]; % Latitude to evaluate R1 seasonality

    for l = {'lo'}; land = l{1};
        land_text = 'Land + Ocean';

        f_vec = par.era.fw;
        for f = f_vec; fw = f{1};
            for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
                folder = sprintf('%s/dmse_comp/%s/%s/0_lat_%g', par.plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                era_ra = era_flux_z.flux_z.(land).ra.(fw);
                era_ra_lat = interp1(era_grid.grid.dim3.lat, era_ra, lat_eval);
                era_res = era_flux_z.flux_z.(land).res.(fw);
                era_res_lat = interp1(era_grid.grid.dim3.lat, era_res, lat_eval);
                era_stf = era_flux_z.flux_z.(land).stf.(fw);
                era_stf_lat = interp1(era_grid.grid.dim3.lat, era_stf, lat_eval);

                gcm_ra = gcm_flux_z.flux_z.(land).ra.mse;
                gcm_ra_lat = interp1(gcm_grid.grid.dim3.lat, gcm_ra, lat_eval);
                gcm_res = gcm_flux_z.flux_z.(land).res.mse;
                gcm_res_lat = interp1(gcm_grid.grid.dim3.lat, gcm_res, lat_eval);
                gcm_stf = gcm_flux_z.flux_z.(land).stf.mse;
                gcm_stf_lat = interp1(gcm_grid.grid.dim3.lat, gcm_stf, lat_eval);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                era_ra=plot([1:12], era_ra_lat, 'k');
                era_res=plot([1:12], era_res_lat, 'color', par.maroon);
                era_stf=plot([1:12], era_stf_lat, 'color', par.blue);
                gcm_ra=plot([1:12], gcm_ra_lat, '--k');
                gcm_res=plot([1:12], gcm_res_lat, '--', 'color', par.maroon);
                gcm_stf=plot([1:12], gcm_stf_lat, '--', 'color', par.blue);
                title(sprintf('ERA-I (solid), ECHAM (dashed), $\\phi=%g^\\circ$', lat_eval));
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([era_ra era_res era_stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                era_ra_ann = repmat(nanmean(era_flux_z.flux_z.(land).ra.(fw), 2), [1 12]);
                era_ra_ann_lat = interp1(era_grid.grid.dim3.lat, era_ra_ann, lat_eval);
                era_dra = era_flux_z.flux_z.(land).ra.(fw) - era_ra_ann;
                era_dra_lat = interp1(era_grid.grid.dim3.lat, era_dra, lat_eval);
                era_res_ann = repmat(nanmean(era_flux_z.flux_z.(land).res.(fw), 2), [1 12]);
                era_res_ann_lat = interp1(era_grid.grid.dim3.lat, era_res_ann, lat_eval);
                era_dres = era_flux_z.flux_z.(land).res.(fw) - era_res_ann;
                era_dres_lat = interp1(era_grid.grid.dim3.lat, era_dres, lat_eval);
                era_stf_ann = repmat(nanmean(era_flux_z.flux_z.(land).stf.(fw), 2), [1 12]);
                era_stf_ann_lat = interp1(era_grid.grid.dim3.lat, era_stf_ann, lat_eval);
                era_dstf = era_flux_z.flux_z.(land).stf.(fw) - era_stf_ann;
                era_dstf_lat = interp1(era_grid.grid.dim3.lat, era_dstf, lat_eval);

                gcm_ra_ann = repmat(nanmean(gcm_flux_z.flux_z.(land).ra.mse, 2), [1 12]);
                gcm_ra_ann_lat = interp1(gcm_grid.grid.dim3.lat, gcm_ra_ann, lat_eval);
                gcm_dra = gcm_flux_z.flux_z.(land).ra.mse - gcm_ra_ann;
                gcm_dra_lat = interp1(gcm_grid.grid.dim3.lat, gcm_dra, lat_eval);
                gcm_res_ann = repmat(nanmean(gcm_flux_z.flux_z.(land).res.mse, 2), [1 12]);
                gcm_res_ann_lat = interp1(gcm_grid.grid.dim3.lat, gcm_res_ann, lat_eval);
                gcm_dres = gcm_flux_z.flux_z.(land).res.mse - gcm_res_ann;
                gcm_dres_lat = interp1(gcm_grid.grid.dim3.lat, gcm_dres, lat_eval);
                gcm_stf_ann = repmat(nanmean(gcm_flux_z.flux_z.(land).stf.mse, 2), [1 12]);
                gcm_stf_ann_lat = interp1(gcm_grid.grid.dim3.lat, gcm_stf_ann, lat_eval);
                gcm_dstf = gcm_flux_z.flux_z.(land).stf.mse - gcm_stf_ann;
                gcm_dstf_lat = interp1(gcm_grid.grid.dim3.lat, gcm_dstf, lat_eval);

                % DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                era_dra=plot([1:12],  era_dra_lat, 'color', 0.5*[1 1 1]);
                era_dres=plot([1:12], era_dres_lat, 'color', par.maroon);
                era_dstf=plot([1:12], era_dstf_lat, 'color', par.blue);
                gcm_dra=plot([1:12],  gcm_dra_lat,  '--', 'color', 0.5*[1 1 1]);
                gcm_dres=plot([1:12], gcm_dres_lat, '--', 'color', par.maroon);
                gcm_dstf=plot([1:12], gcm_dstf_lat, '--', 'color', par.blue);
                title(sprintf('ERA-I (solid), ECHAM (dashed), $\\phi=%g^\\circ$', lat_eval));
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([era_dra era_dres era_dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
function plot_temp_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    lat=length(grid_pl.grid.dim3.lat);
    % lat=1;
    lon=1;
    mon=1;

    folder = sprintf('%s/temp/mon_%g/lat_%g/lon_%g', par.plotdir, mon, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    grid_pl.grid.dim3.lon(lon)
    ps_pl(lon,lat,mon)

    ta_pl0=squeeze(ta_pl(lon,lat,:,mon));
    si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,mon)).*(grid_pl.grid.dim3.plev<squeeze(ps_pl(lon,lat,mon)));
    ta_pl0(find(~si_pl0)) = [];
    si_pl0(find(~si_pl0)) = [];

    ta_ml0=squeeze(ta_ml(lon,lat,:,mon));
    si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,mon))+grid_ml.grid.dim3.b;

    figure(); clf; hold all; box on;
    si_pl0(1)
    si_ml0(end)
    h_pl = plot(ta_pl0, si_pl0, '-*');
    h_ml = plot(ta_ml0, si_ml0, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp', folder), '-dpng', '-r300');
    close;

    % interpolation comparisons
    % add surface data
    ta_pl0(end+1) = squeeze(tas_pl(lon,lat,mon));
    si_pl0(end+1) = 1;
    ta_pl0 = circshift(ta_pl0,1);
    si_pl0 = circshift(si_pl0,1);

    ta_ml0(end+1) = squeeze(tas_ml(lon,lat,mon));
    si_ml0(end+1) = 1;

    % LINEAR
    ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
    ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/lin_comp', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/lin_diff', folder), '-dpng', '-r300');
    close;

    % CUBIC
    ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
    ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/cub_comp', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/cub_diff', folder), '-dpng', '-r300');
    close;

    % SPLINE
    ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
    ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/spl_comp', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/spl_diff', folder), '-dpng', '-r300');
    close;

    % MODIFIED AKIMA
    ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
    ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/mak_comp', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/mak_diff', folder), '-dpng', '-r300');
    close;

    % ALL PL
    figure(); clf; hold all; box on;
    plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
    hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
    close;

    % ALL ML
    figure(); clf; hold all; box on;
    hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
    hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
    close;

    % ALL DIFF
    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
    close;

end
function plot_temp_pl2_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl2/ATM_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*.ymonmean.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    lat=length(grid_pl.grid.dim3.lat);
    % lat=1;
    lon=1;
    mon=1;

    folder = sprintf('%s/temp/mon_%g/lat_%g/lon_%g', par.plotdir, mon, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    grid_pl.grid.dim3.lon(lon)
    ps_pl(lon,lat,mon)

    ta_ml0=squeeze(ta_ml(lon,lat,:,mon));
    si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,mon))+grid_ml.grid.dim3.b;

    ta_pl0=squeeze(ta_pl(lon,lat,:,mon));
    si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,mon)).*(grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,mon))<si_ml0(end));
    ta_pl0(find(~si_pl0)) = [];
    si_pl0(find(~si_pl0)) = [];

    figure(); clf; hold all; box on;
    si_pl0(1)
    si_ml0(end)
    h_pl = plot(ta_pl0, si_pl0, '-*');
    h_ml = plot(ta_ml0, si_ml0, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp2', folder), '-dpng', '-r300');
    close;

    % interpolation comparisons
    % add surface data
    ta_pl0(end+1) = squeeze(tas_pl(lon,lat,mon));
    si_pl0(end+1) = 1;
    ta_pl0 = circshift(ta_pl0,1);
    si_pl0 = circshift(si_pl0,1);

    ta_ml0(end+1) = squeeze(tas_ml(lon,lat,mon));
    si_ml0(end+1) = 1;

    % LINEAR
    ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
    ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/lin_comp2', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/lin_diff2', folder), '-dpng', '-r300');
    close;

    % CUBIC
    ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
    ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/cub_comp2', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/cub_diff2', folder), '-dpng', '-r300');
    close;

    % SPLINE
    ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
    ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/spl_comp2', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/spl_diff2', folder), '-dpng', '-r300');
    close;

    % MODIFIED AKIMA
    ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
    ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/mak_comp2', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/mak_diff2', folder), '-dpng', '-r300');
    close;

    % ALL PL
    figure(); clf; hold all; box on;
    plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
    hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_pl2', folder), '-dpng', '-r300');
    close;

    % ALL ML
    figure(); clf; hold all; box on;
    hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
    hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_ml2', folder), '-dpng', '-r300');
    close;

    % ALL DIFF
    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_diff2', folder), '-dpng', '-r300');
    close;

end
function plot_temp_daily_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    lat=length(grid_pl.grid.dim3.lat);
    % lat=1;
    day=1;
    mon=1;

    folder = sprintf('%s/temp/mon_%g/day_%g/lat_%g', par.plotdir, mon, day, grid_pl.grid.dim3.lat(lat));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    ta_pl0=squeeze(ta_pl(:,lat,:,day));
    pa_pl = repmat(grid_pl.grid.dim3.plev', [length(grid_pl.grid.dim3.lon) 1]);
    ps_vert_pl = repmat(squeeze(ps_pl(:,lat,day)), [1 length(grid_pl.grid.dim3.plev)]);
    si_pl0=pa_pl./ps_vert_pl.*(pa_pl<ps_vert_pl);
    for lon = 1:length(grid_pl.grid.dim3.lon)
        ta_tmp = ta_pl0(lon, :);
        si_tmp = si_pl0(lon, :);
        ta_tmp(find(~si_tmp)) = [];
        si_tmp(find(~si_tmp)) = [];
        ta_pli(lon, :) = interp1(si_tmp, ta_tmp, grid_pl.grid.dim3.si);
    end

    ta_ml0=squeeze(ta_ml(:,lat,:,day));
    a_ml = repmat(grid_ml.grid.dim3.a', [length(grid_pl.grid.dim3.lon) 1]);
    b_ml = repmat(grid_ml.grid.dim3.b', [length(grid_pl.grid.dim3.lon) 1]);
    ps_vert_ml = repmat(squeeze(ps_pl(:,lat,day)), [1 length(grid_ml.grid.dim3.a)]);
    si_ml0=a_ml./ps_vert_ml+b_ml;
    for lon = 1:length(grid_ml.grid.dim3.lon)
        ta_tmp = ta_ml0(lon, :);
        si_tmp = si_ml0(lon, :);
        ta_tmp(find(~si_tmp)) = [];
        si_tmp(find(~si_tmp)) = [];
        ta_mli(lon, :) = interp1(si_tmp, ta_tmp, grid_ml.grid.dim3.si);
    end

    % zonal mean
    ta_pli = squeeze(mean(ta_pli, 1));
    ta_mli = squeeze(mean(ta_mli, 1));

    % nanmean zonal mean
    ta_pli_nanmean = squeeze(nanmean(ta_pli, 1));
    ta_mli_nanmean = squeeze(nanmean(ta_mli, 1));

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pli, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_mli, grid_ml.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pli_nanmean, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_mli_nanmean, grid_ml.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp_nanmean', folder), '-dpng', '-r300');
    close;

    for day=1:size(ta_pl,4)

        lat=length(grid_pl.grid.dim3.lat);
        % lat=1;
        lon=1;
        mon=1;
        % day=1;

        folder = sprintf('%s/temp/mon_%g/day_%g/lat_%g/lon_%g', par.plotdir, mon, day, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
        if ~exist(folder, 'dir'); mkdir(folder); end;

        grid_pl.grid.dim3.lon(lon)
        ps_pl(lon,lat,day)

        ta_pl0=squeeze(ta_pl(lon,lat,:,day));
        si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day)).*(grid_pl.grid.dim3.plev<squeeze(ps_pl(lon,lat,day)));
        ta_pl0(find(~si_pl0)) = [];
        si_pl0(find(~si_pl0)) = [];

        ta_ml0=squeeze(ta_ml(lon,lat,:,day));
        si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,day))+grid_ml.grid.dim3.b;

        figure(); clf; hold all; box on;
        si_pl0(1)
        si_ml0(end)
        h_pl = plot(ta_pl0, si_pl0, '-*');
        h_ml = plot(ta_ml0, si_ml0, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/comp', folder), '-dpng', '-r300');
        close;

        % interpolation comparisons
        % add surface data
        ta_pl0(end+1) = squeeze(tas_pl(lon,lat,day));
        si_pl0(end+1) = 1;
        ta_pl0 = circshift(ta_pl0,1);
        si_pl0 = circshift(si_pl0,1);

        ta_ml0(end+1) = squeeze(tas_ml(lon,lat,day));
        si_ml0(end+1) = 1;

        % LINEAR
        ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
        ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_comp', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_diff', folder), '-dpng', '-r300');
        close;

        % CUBIC
        ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
        ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_comp', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_diff', folder), '-dpng', '-r300');
        close;

        % SPLINE
        ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
        ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_comp', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_diff', folder), '-dpng', '-r300');
        close;

        % MODIFIED AKIMA
        ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
        ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_comp', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_diff', folder), '-dpng', '-r300');
        close;

        % ALL PL
        figure(); clf; hold all; box on;
        plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
        hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
        close;

        % ALL ML
        figure(); clf; hold all; box on;
        hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
        hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
        close;

        % ALL DIFF
        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
        close;
    end

end
function plot_temp_daily_avg_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    lat=length(grid_pl.grid.dim3.lat);
    % lat=1;
    lon=1;
    mon=1;

    ta_ml0=squeeze(ta_ml(lon,lat,:,:));
    a_time = repmat(grid_ml.grid.dim3.a, [1 size(ta_ml0, 2)]);
    b_time = repmat(grid_ml.grid.dim3.b, [1 size(ta_ml0, 2)]);
    ps_ml = repmat(permute(squeeze(ps_ml(lon,lat,:)), [2 1]), [size(a_time,1) 1]);
    si_ml0=a_time./ps_ml+grid_ml.grid.dim3.b;
    for day = 1:size(ta_ml0, 4)
        ta_tmp = ta_ml0(:, day);
        si_tmp = si_ml0(:, day);
        ta_tmp(find(~si_tmp)) = [];
        si_tmp(find(~si_tmp)) = [];
        ta_mli(:, day) = interp1(si_tmp, ta_tmp, grid_ml.grid.dim3.si);
    end

    pa = repmat(grid_pl.grid.dim3.plev, [1 size(ta_pl,4)]);
    ps = repmat(squeeze(ps_pl(lon,lat,:))', [size(pa, 1) 1]);
    ta_pl0=squeeze(ta_pl(lon,lat,:,:));
    si_pl0_mask=pa./ps.*((pa./ps)<repmat(si_ml0(end,:), [size(pa, 1) 1]));
    si_pl0=pa./ps;
    % ta_pl0(find(~si_pl0)) = nan;
    % si_pl0(find(~si_pl0)) = nan;
    % si_pl0=pa_pl./ps_vert_pl.*(pa_pl<ps_vert_pl);
    for day = 1:size(ta_pl0, 4)
        ta_tmp = ta_pl0(:, day);
        si_tmp = si_pl0(:, day);
        si_mask_tmp = si_pl0_mask(:, day);
        ta_tmp(find(~si_mask_tmp)) = nan;
        ta_pli(:, day) = interp1(si_tmp, ta_tmp, grid_pl.grid.dim3.si, 'linear', nan)
    end

    ta_ml_mean = mean(ta_mli, 2);
    ta_pl_mean = mean(ta_pli, 2);
    ta_ml_nanmean = nanmean(ta_mli, 2);
    ta_pl_nanmean = nanmean(ta_pli, 2);

    folder = sprintf('%s/temp/mon_%g/lat_%g/lon_%g', par.plotdir, mon, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_mean, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_mean, grid_ml.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp_mean', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_nanmean, grid_pl.grid.dim3.si, '-');
    h_ml = plot(ta_ml_nanmean, grid_ml.grid.dim3.si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/comp_nanmean', folder), '-dpng', '-r300');
    close;

end
function plot_temp_daily_mod_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    for day=1:size(ta_pl,4)

        lat=length(grid_pl.grid.dim3.lat);
        % lat=1;
        lon=1;
        mon=1;
        % day=1;

        folder = sprintf('%s/temp/mon_%g/day_%g/lat_%g/lon_%g', par.plotdir, mon, day, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
        if ~exist(folder, 'dir'); mkdir(folder); end;

        grid_pl.grid.dim3.lon(lon)
        ps_pl(lon,lat,day)

        ta_ml0=squeeze(ta_ml(lon,lat,:,day));
        si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,day))+grid_ml.grid.dim3.b;

        ta_pl0=squeeze(ta_pl(lon,lat,:,day));
        si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day)).*(grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day))<si_ml0(end));
        ta_pl0(find(~si_pl0)) = [];
        si_pl0(find(~si_pl0)) = [];

        figure(); clf; hold all; box on;
        si_pl0(1)
        si_ml0(end)
        h_pl = plot(ta_pl0, si_pl0, '-*');
        h_ml = plot(ta_ml0, si_ml0, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/comp_mod', folder), '-dpng', '-r300');
        close;

        % interpolation comparisons
        % add surface data
        ta_pl0(end+1) = squeeze(tas_pl(lon,lat,day));
        si_pl0(end+1) = 1;
        ta_pl0 = circshift(ta_pl0,1);
        si_pl0 = circshift(si_pl0,1);

        ta_ml0(end+1) = squeeze(tas_ml(lon,lat,day));
        si_ml0(end+1) = 1;

        % LINEAR
        ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
        ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_comp_mod', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_diff_mod', folder), '-dpng', '-r300');
        close;

        % CUBIC
        ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
        ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_comp_mod', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_diff_mod', folder), '-dpng', '-r300');
        close;

        % SPLINE
        ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
        ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_comp_mod', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_diff_mod', folder), '-dpng', '-r300');
        close;

        % MODIFIED AKIMA
        ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
        ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_comp_mod', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_diff_mod', folder), '-dpng', '-r300');
        close;

        % ALL PL
        figure(); clf; hold all; box on;
        plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
        hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_pl_mod', folder), '-dpng', '-r300');
        close;

        % ALL ML
        figure(); clf; hold all; box on;
        hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
        hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_ml_mod', folder), '-dpng', '-r300');
        close;

        % ALL DIFF_MOD
        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_diff_mod', folder), '-dpng', '-r300');
        close;
    end

end
function plot_temp_daily_pl2_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl2/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    for day=1:size(ta_pl,4)

        lat=length(grid_pl.grid.dim3.lat);
        % lat=1;
        lon=1;
        mon=1;
        % day=1;

        folder = sprintf('%s/temp/mon_%g/day_%g/lat_%g/lon_%g', par.plotdir, mon, day, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
        if ~exist(folder, 'dir'); mkdir(folder); end;

        grid_pl.grid.dim3.lon(lon)
        ps_pl(lon,lat,day)

        ta_pl0=squeeze(ta_pl(lon,lat,:,day));
        si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day)).*(grid_pl.grid.dim3.plev<squeeze(ps_pl(lon,lat,day)));
        ta_pl0(find(~si_pl0)) = [];
        si_pl0(find(~si_pl0)) = [];

        ta_ml0=squeeze(ta_ml(lon,lat,:,day));
        si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,day))+grid_ml.grid.dim3.b;

        figure(); clf; hold all; box on;
        si_pl0(1)
        si_ml0(end)
        h_pl = plot(ta_pl0, si_pl0, '-*');
        h_ml = plot(ta_ml0, si_ml0, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/comp2', folder), '-dpng', '-r300');
        close;

        % interpolation comp2arisons
        % add surface data
        ta_pl0(end+1) = squeeze(tas_pl(lon,lat,day));
        si_pl0(end+1) = 1;
        ta_pl0 = circshift(ta_pl0,1);
        si_pl0 = circshift(si_pl0,1);

        ta_ml0(end+1) = squeeze(tas_ml(lon,lat,day));
        si_ml0(end+1) = 1;

        % LINEAR
        ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
        ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_comp2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_diff2', folder), '-dpng', '-r300');
        close;

        % CUBIC
        ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
        ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_comp2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_diff2', folder), '-dpng', '-r300');
        close;

        % SPLINE
        ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
        ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_comp2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_diff2', folder), '-dpng', '-r300');
        close;

        % MODIFIED AKIMA
        ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
        ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_comp2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_diff2', folder), '-dpng', '-r300');
        close;

        % ALL PL
        figure(); clf; hold all; box on;
        plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
        hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_pl2', folder), '-dpng', '-r300');
        close;

        % ALL ML
        figure(); clf; hold all; box on;
        hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
        hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_ml2', folder), '-dpng', '-r300');
        close;

        % ALL DIFF2
        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_diff2', folder), '-dpng', '-r300');
        close;
    end

end
function plot_temp_daily_mod_pl2_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_ml = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_ml = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_ml = double(ncread(fullpath, var));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    var = 't';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl2/ATM_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ta_pl = double(ncread(fullpath, var));
    var = 'aps';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    ps_pl = double(ncread(fullpath, var));
    var = 'temp2';
    file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_pl/BOT_dm_*1949.jan.nc'));
    fullpath=sprintf('%s/%s', file.folder, file.name);
    tas_pl = double(ncread(fullpath, var));

    for day=1:size(ta_pl,4)

        lat=length(grid_pl.grid.dim3.lat);
        % lat=1;
        lon=1;
        mon=1;
        % day=1;

        folder = sprintf('%s/temp/mon_%g/day_%g/lat_%g/lon_%g', par.plotdir, mon, day, grid_pl.grid.dim3.lat(lat), grid_pl.grid.dim3.lon(lon));
        if ~exist(folder, 'dir'); mkdir(folder); end;

        grid_pl.grid.dim3.lon(lon)
        ps_pl(lon,lat,day)

        ta_ml0=squeeze(ta_ml(lon,lat,:,day));
        si_ml0=grid_ml.grid.dim3.a/squeeze(ps_ml(lon,lat,day))+grid_ml.grid.dim3.b;

        ta_pl0=squeeze(ta_pl(lon,lat,:,day));
        si_pl0=grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day)).*(grid_pl.grid.dim3.plev/squeeze(ps_pl(lon,lat,day))<si_ml0(end));
        ta_pl0(find(~si_pl0)) = [];
        si_pl0(find(~si_pl0)) = [];

        figure(); clf; hold all; box on;
        si_pl0(1)
        si_ml0(end)
        h_pl = plot(ta_pl0, si_pl0, '-*');
        h_ml = plot(ta_ml0, si_ml0, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/comp_mod_pl2', folder), '-dpng', '-r300');
        close;

        % interpolation comparisons
        % add surface data
        ta_pl0(end+1) = squeeze(tas_pl(lon,lat,day));
        si_pl0(end+1) = 1;
        ta_pl0 = circshift(ta_pl0,1);
        si_pl0 = circshift(si_pl0,1);

        ta_ml0(end+1) = squeeze(tas_ml(lon,lat,day));
        si_ml0(end+1) = 1;

        % LINEAR
        ta_pl_lin = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'linear', nan);
        ta_ml_lin = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'linear', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_lin, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_lin, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_comp_mod_pl2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/lin_diff_mod_pl2', folder), '-dpng', '-r300');
        close;

        % CUBIC
        ta_pl_cub = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'pchip', nan);
        ta_ml_cub = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'pchip', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_cub, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_cub, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_comp_mod_pl2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/cub_diff_mod_pl2', folder), '-dpng', '-r300');
        close;

        % SPLINE
        ta_pl_spl = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'spline', nan);
        ta_ml_spl = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'spline', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_spl, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_spl, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_comp_mod_pl2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/spl_diff_mod_pl2', folder), '-dpng', '-r300');
        close;

        % MOD_PL2IFIED AKIMA
        ta_pl_mak = interp1(si_pl0, ta_pl0, grid_pl.grid.dim3.si, 'makima', nan);
        ta_ml_mak = interp1(si_ml0, ta_ml0, grid_pl.grid.dim3.si, 'makima', nan);

        figure(); clf; hold all; box on;
        h_pl = plot(ta_pl_mak, grid_pl.grid.dim3.si, '-');
        h_ml = plot(ta_ml_mak, grid_pl.grid.dim3.si, '--');
        xlabel('T (K)');
        ylabel('$\sigma$ (unitless)');
        legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_comp_mod_pl2', folder), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, 'k')
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/mak_diff_mod_pl2', folder), '-dpng', '-r300');
        close;

        % ALL PL
        figure(); clf; hold all; box on;
        plot(ta_ml_lin, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
        hlin=plot(ta_pl_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Mod_Pl2ified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_pl_mod_pl2', folder), '-dpng', '-r300');
        close;

        % ALL ML
        figure(); clf; hold all; box on;
        hlin=plot(ta_ml_lin, grid_ml.grid.dim3.si, '-k');
        hcub=plot(ta_ml_cub, grid_ml.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_ml_spl, grid_ml.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_ml_mak, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Mod_Pl2ified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_ml_mod_pl2', folder), '-dpng', '-r300');
        close;

        % ALL DIFF_MOD_PL2
        figure(); clf; hold all; box on;
        line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
        hlin=plot(ta_pl_lin - ta_ml_lin, grid_pl.grid.dim3.si, '-k');
        hcub=plot(ta_pl_cub - ta_ml_cub, grid_pl.grid.dim3.si, '--', 'color', par.blue);
        hspl=plot(ta_pl_spl - ta_ml_spl, grid_pl.grid.dim3.si, ':', 'color', par.orange);
        hmak=plot(ta_pl_mak - ta_ml_mak, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
        xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
        ylabel('$\sigma$ (unitless)');
        legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Mod_Pl2ified Akima', 'location', 'eastoutside');
        title(sprintf('ECHAM, Month=%g, Day=%g, Lat=$%g^\\circ$, Lon=$%g^\\circ$', mon, day, grid_ml.grid.dim3.lat(lat), grid_ml.grid.dim3.lon(lon)));
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
        set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
        print(sprintf('%s/all_diff_mod_pl2', folder), '-dpng', '-r300');
        close;
    end

end
function plot_temp_zon_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', 'echam_pl', par.lat_interp));

    lat = grid_ml.grid.dim3.lat;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            ta_ml.tasi_z.(land).(time) = squeeze(mean(ta_ml.tasi_t.(land).(time), 1)); % zonal average
            ta_pl.tasi_z.(land).(time) = squeeze(mean(ta_pl.tasi_t.(land).(time), 1)); % zonal average

            idx_f = min(find(~isnan(ta_ml.tasi_z.(land).(time)(:,floor(size(ta_ml.tasi_z.(land).(time),2)/2))))); % first non-nan latitude index
            [~, idx_fi] = min(abs(lat - ceil(lat(idx_f)))); % first non-nan integer latitude index
            step = 1; % step size of adjacent latitudes
            nums = 9; % number of latitudes to plot at once

            if par.do_surf; v_vec = {'si'};
            else v_vec = {'si'}; end
            for v = v_vec; vert = v{1};
                for h = {'sh', 'nh', 'eq', 'nmid', 'smid', 'nhall', 'shall'}; hemi = h{1};
                    figure(); clf; hold all; box on;
                    if strcmp(h, 'sh'); indices = idx_fi:step:step*(nums+1);
                    elseif strcmp(h, 'nh'); indices = length(lat)-idx_fi+1:-step:length(lat)-step*nums;
                    elseif strcmp(h, 'eq');
                        if mod(length(lat),2)==0; indices = length(lat)/2+step*floor(nums/2):-step:length(lat)/2-step*floor(nums/2);
                        elseif mod(length(lat),2)==1; indices = floor(length(lat)/2)+1+step*floor(nums/2):-step:floor(length(lat)/2)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nmid');
                        if mod(length(lat),4)==0; indices = 3*length(lat)/4+step*floor(nums/2):-step:3*length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(3*length(lat)/4)+1+step*floor(nums/2):-step:floor(3*length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'smid');
                        if mod(length(lat),4)==0; indices = length(lat)/4+step*floor(nums/2):-step:length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(length(lat)/4)+1+step*floor(nums/2):-step:floor(length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nhall');
                        if mod(length(lat),2)==0; indices = length(lat)-idx_fi+1:-step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = length(lat)-idx_fi+1:-step:floor(length(lat)/2); end;
                    elseif strcmp(h, 'shall');
                        if mod(length(lat),2)==0; indices = idx_fi:step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = idx_fi:step:floor(length(lat)/2); end;
                    end
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, '--');
                        else
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, 'color', (c/(1.3*nums))*[1 1 1]);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, '--', 'color', (c/(1.3*nums))*[1 1 1]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, '.');
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, 'x');
                        else
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, '.', 'color', (c/(1.3*nums))*[1 1 1]);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, 'x', 'color', (c/(1.3*nums))*[1 1 1]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_points_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(squeeze(ta_ml.tasi_t.(land).(time)(2,i,:)), grid_ml.grid.dim3.si);
                            h{c} = plot(squeeze(ta_pl.tasi_t.(land).(time)(2,i,:)), grid_pl.grid.dim3.si, '--');
                        else
                            h{c} = plot(squeeze(ta_ml.tasi_t.(land).(time)(2,i,:)), grid_ml.grid.dim3.si, 'color', (c/(2.3*nums))*[2 2 2]);
                            h{c} = plot(squeeze(ta_pl.tasi_t.(land).(time)(2,i,:)), grid_pl.grid.dim3.si, '--', 'color', (c/(2.3*nums))*[2 2 2]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_lon_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;
                    clear h h_label;

                end % hemi
            end % vert

        end
    end
end
function plot_temp_zonnanmean_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', 'echam_pl', par.lat_interp));

    lat = grid_ml.grid.dim3.lat;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            ta_ml.tasi_z.(land).(time) = squeeze(nanmean(ta_ml.tasi_t.(land).(time), 1)); % zonal average
            ta_pl.tasi_z.(land).(time) = squeeze(nanmean(ta_pl.tasi_t.(land).(time), 1)); % zonal average

            idx_f = min(find(~isnan(ta_ml.tasi_z.(land).(time)(:,floor(size(ta_ml.tasi_z.(land).(time),2)/2))))); % first non-nan latitude index
            [~, idx_fi] = min(abs(lat - ceil(lat(idx_f)))); % first non-nan integer latitude index
            step = 1; % step size of adjacent latitudes
            nums = 9; % number of latitudes to plot at once

            if par.do_surf; v_vec = {'si'};
            else v_vec = {'si'}; end
            for v = v_vec; vert = v{1};
                for h = {'sh', 'nh', 'eq', 'nmid', 'smid', 'nhall', 'shall'}; hemi = h{1};
                    figure(); clf; hold all; box on;
                    if strcmp(h, 'sh'); indices = idx_fi:step:step*(nums+1);
                    elseif strcmp(h, 'nh'); indices = length(lat)-idx_fi+1:-step:length(lat)-step*nums;
                    elseif strcmp(h, 'eq');
                        if mod(length(lat),2)==0; indices = length(lat)/2+step*floor(nums/2):-step:length(lat)/2-step*floor(nums/2);
                        elseif mod(length(lat),2)==1; indices = floor(length(lat)/2)+1+step*floor(nums/2):-step:floor(length(lat)/2)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nmid');
                        if mod(length(lat),4)==0; indices = 3*length(lat)/4+step*floor(nums/2):-step:3*length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(3*length(lat)/4)+1+step*floor(nums/2):-step:floor(3*length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'smid');
                        if mod(length(lat),4)==0; indices = length(lat)/4+step*floor(nums/2):-step:length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(length(lat)/4)+1+step*floor(nums/2):-step:floor(length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nhall');
                        if mod(length(lat),2)==0; indices = length(lat)-idx_fi+1:-step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = length(lat)-idx_fi+1:-step:floor(length(lat)/2); end;
                    elseif strcmp(h, 'shall');
                        if mod(length(lat),2)==0; indices = idx_fi:step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = idx_fi:step:floor(length(lat)/2); end;
                    end
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, '--');
                        else
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, 'color', (c/(1.3*nums))*[1 1 1]);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, '--', 'color', (c/(1.3*nums))*[1 1 1]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, '.');
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, 'x');
                        else
                            h{c} = plot(ta_ml.tasi_z.(land).(time)(i,:), grid_ml.grid.dim3.si, '.', 'color', (c/(1.3*nums))*[1 1 1]);
                            h{c} = plot(ta_pl.tasi_z.(land).(time)(i,:), grid_pl.grid.dim3.si, 'x', 'color', (c/(1.3*nums))*[1 1 1]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_points_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            h{c} = plot(squeeze(ta_ml.tasi_t.(land).(time)(2,i,:)), grid_ml.grid.dim3.si);
                            h{c} = plot(squeeze(ta_pl.tasi_t.(land).(time)(2,i,:)), grid_pl.grid.dim3.si, '--');
                        else
                            h{c} = plot(squeeze(ta_ml.tasi_t.(land).(time)(2,i,:)), grid_ml.grid.dim3.si, 'color', (c/(2.3*nums))*[2 2 2]);
                            h{c} = plot(squeeze(ta_pl.tasi_t.(land).(time)(2,i,:)), grid_pl.grid.dim3.si, '--', 'color', (c/(2.3*nums))*[2 2 2]);
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/comp_lon_%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;
                    clear h h_label;

                end % hemi
            end % vert

        end
    end
end
function plot_temp_zon_interp_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', type, par.lat_interp));
    ta_ml_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', type, par.lat_interp));
    ta_ml_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', type, par.lat_interp));
    ta_ml_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mean_mon_lat_lin.mat', 'echam_pl', par.lat_interp));
    ta_pl_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mean_mon_lat_cub.mat', 'echam_pl', par.lat_interp));
    ta_pl_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mean_mon_lat_spl.mat', 'echam_pl', par.lat_interp));
    ta_pl_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mean_mon_lat_mak.mat', 'echam_pl', par.lat_interp));

    si = grid_pl.grid.dim3.si;

    lat = length(grid_pl.grid.dim3.lat);
    % lat = 1;
    mon = 1;

    folder = sprintf('%s/temp/mon_%g/lat_%g', par.plotdir, mon, grid_pl.grid.dim3.lat(lat));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    ta_pl_lin0 = squeeze(ta_pl_lin.tasi.lo(lat,mon,:));
    ta_ml_lin0 = squeeze(ta_ml_lin.tasi.lo(lat,mon,:));
    ta_pl_cub0 = squeeze(ta_pl_cub.tasi.lo(lat,mon,:));
    ta_ml_cub0 = squeeze(ta_ml_cub.tasi.lo(lat,mon,:));
    ta_pl_spl0 = squeeze(ta_pl_spl.tasi.lo(lat,mon,:));
    ta_ml_spl0 = squeeze(ta_ml_spl.tasi.lo(lat,mon,:));
    ta_pl_mak0 = squeeze(ta_pl_mak.tasi.lo(lat,mon,:));
    ta_ml_mak0 = squeeze(ta_ml_mak.tasi.lo(lat,mon,:));

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_lin0, si, '-');
    h_ml = plot(ta_ml_lin0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Linear', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zon_comp_lin', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_cub0, si, '-');
    h_ml = plot(ta_ml_cub0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Piecewise Cubic Hermite', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zon_comp_cub', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_spl0, si, '-');
    h_ml = plot(ta_ml_spl0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Cubic Spline', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zon_comp_spl', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_mak0, si, '-');
    h_ml = plot(ta_ml_mak0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Modified Akima', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zon_comp_mak', folder), '-dpng', '-r300');
    close;

    % ALL PL
    figure(); clf; hold all; box on;
    plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
    hlin=plot(ta_pl_lin0, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
    close;

    % ALL ML
    figure(); clf; hold all; box on;
    hlin=plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-k');
    hcub=plot(ta_ml_cub0, grid_ml.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_ml_spl0, grid_ml.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_ml_mak0, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
    close;

    % ALL DIFF
    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    hlin=plot(ta_pl_lin0 - ta_ml_lin0, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub0 - ta_ml_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl0 - ta_ml_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak0 - ta_ml_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
    close;

end
function plot_temp_zonnanmean_interp_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', type, par.lat_interp));
    ta_ml_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', type, par.lat_interp));
    ta_ml_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', type, par.lat_interp));
    ta_ml_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', 'echam_pl', par.lat_interp));
    ta_pl_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', 'echam_pl', par.lat_interp));
    ta_pl_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', 'echam_pl', par.lat_interp));
    ta_pl_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', 'echam_pl', par.lat_interp));

    si = grid_pl.grid.dim3.si;

    lat = length(grid_pl.grid.dim3.lat);
    % lat = 1;
    mon = 1;

    folder = sprintf('%s/temp/mon_%g/lat_%g', par.plotdir, mon, grid_pl.grid.dim3.lat(lat));
    if ~exist(folder, 'dir'); mkdir(folder); end;

    ta_pl_lin0 = squeeze(ta_pl_lin.tasi.lo(lat,mon,:));
    ta_ml_lin0 = squeeze(ta_ml_lin.tasi.lo(lat,mon,:));
    ta_pl_cub0 = squeeze(ta_pl_cub.tasi.lo(lat,mon,:));
    ta_ml_cub0 = squeeze(ta_ml_cub.tasi.lo(lat,mon,:));
    ta_pl_spl0 = squeeze(ta_pl_spl.tasi.lo(lat,mon,:));
    ta_ml_spl0 = squeeze(ta_ml_spl.tasi.lo(lat,mon,:));
    ta_pl_mak0 = squeeze(ta_pl_mak.tasi.lo(lat,mon,:));
    ta_ml_mak0 = squeeze(ta_ml_mak.tasi.lo(lat,mon,:));

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_lin0, si, '-');
    h_ml = plot(ta_ml_lin0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Linear', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zonnanmean_comp_lin', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_cub0, si, '-');
    h_ml = plot(ta_ml_cub0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Piecewise Cubic Hermite', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zonnanmean_comp_cub', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_spl0, si, '-');
    h_ml = plot(ta_ml_spl0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Cubic Spline', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zonnanmean_comp_spl', folder), '-dpng', '-r300');
    close;

    figure(); clf; hold all; box on;
    h_pl = plot(ta_pl_mak0, si, '-');
    h_ml = plot(ta_ml_mak0, si, '--');
    xlabel('T (K)');
    ylabel('$\sigma$ (unitless)');
    legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$, Modified Akima', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/zonnanmean_comp_mak', folder), '-dpng', '-r300');
    close;

    % ALL PL
    figure(); clf; hold all; box on;
    plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
    hlin=plot(ta_pl_lin0, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
    close;

    % ALL ML
    figure(); clf; hold all; box on;
    hlin=plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-k');
    hcub=plot(ta_ml_cub0, grid_ml.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_ml_spl0, grid_ml.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_ml_mak0, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
    close;

    % ALL DIFF
    figure(); clf; hold all; box on;
    line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
    hlin=plot(ta_pl_lin0 - ta_ml_lin0, grid_pl.grid.dim3.si, '-k');
    hcub=plot(ta_pl_cub0 - ta_ml_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
    hspl=plot(ta_pl_spl0 - ta_ml_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
    hmak=plot(ta_pl_mak0 - ta_ml_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
    xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
    ylabel('$\sigma$ (unitless)');
    legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
    title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$', mon, grid_ml.grid.dim3.lat(lat)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
    close;

end
function plot_temp_polar_interp_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', type, par.lat_interp));
    ta_ml_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', type, par.lat_interp));
    ta_ml_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', type, par.lat_interp));
    ta_ml_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', 'echam_pl', par.lat_interp));
    ta_pl_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', 'echam_pl', par.lat_interp));
    ta_pl_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', 'echam_pl', par.lat_interp));
    ta_pl_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', 'echam_pl', par.lat_interp));

    si = grid_pl.grid.dim3.si;

    lat_bound_list = [-85 -80 -70 70 80 85];
    mon_list = 1:12;

    for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
        for m = 1:length(mon_list); mon = mon_list(m)
            dlat = 0.25; % step size for standard lat grid
            if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole;
            else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; end;
            clat = cosd(lat); % cosine of latitude for cosine weighting

            folder = sprintf('%s/temp/mon_%g/poleward_of_lat_%g', par.plotdir, mon, lat_bound);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            ta_pl_lin0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_lin.tasi.lo(:,mon,:), lat));
            ta_ml_lin0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_lin.tasi.lo(:,mon,:), lat));
            ta_pl_cub0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_cub.tasi.lo(:,mon,:), lat));
            ta_ml_cub0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_cub.tasi.lo(:,mon,:), lat));
            ta_pl_spl0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_spl.tasi.lo(:,mon,:), lat));
            ta_ml_spl0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_spl.tasi.lo(:,mon,:), lat));
            ta_pl_mak0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_mak.tasi.lo(:,mon,:), lat));
            ta_ml_mak0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_mak.tasi.lo(:,mon,:), lat));

            clat_vert = repmat(clat', [1 size(ta_pl_lin0, 2)]);

            ta_pl_lin0 = nansum(ta_pl_lin0.*clat_vert, 1)/nansum(clat);
            ta_ml_lin0 = nansum(ta_ml_lin0.*clat_vert, 1)/nansum(clat);
            ta_pl_cub0 = nansum(ta_pl_cub0.*clat_vert, 1)/nansum(clat);
            ta_ml_cub0 = nansum(ta_ml_cub0.*clat_vert, 1)/nansum(clat);
            ta_pl_spl0 = nansum(ta_pl_spl0.*clat_vert, 1)/nansum(clat);
            ta_ml_spl0 = nansum(ta_ml_spl0.*clat_vert, 1)/nansum(clat);
            ta_pl_mak0 = nansum(ta_pl_mak0.*clat_vert, 1)/nansum(clat);
            ta_ml_mak0 = nansum(ta_ml_mak0.*clat_vert, 1)/nansum(clat);

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_lin0, si, '-');
            h_ml = plot(ta_ml_lin0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Linear', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_lin', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_cub0, si, '-');
            h_ml = plot(ta_ml_cub0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Piecewise Cubic Hermite', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_cub', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_spl0, si, '-');
            h_ml = plot(ta_ml_spl0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Cubic Spline', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_spl', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_mak0, si, '-');
            h_ml = plot(ta_ml_mak0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Modified Akima', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_mak', folder), '-dpng', '-r300');
            close;

            % ALL PL
            figure(); clf; hold all; box on;
            plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
            hlin=plot(ta_pl_lin0, grid_pl.grid.dim3.si, '-k');
            hcub=plot(ta_pl_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_pl_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_pl_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
            close;

            % ALL ML
            figure(); clf; hold all; box on;
            hlin=plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-k');
            hcub=plot(ta_ml_cub0, grid_ml.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_ml_spl0, grid_ml.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_ml_mak0, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
            close;

            % ALL DIFF
            figure(); clf; hold all; box on;
            line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
            hlin=plot(ta_pl_lin0 - ta_ml_lin0, grid_pl.grid.dim3.si, '-k');
            hcub=plot(ta_pl_cub0 - ta_ml_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_pl_spl0 - ta_ml_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_pl_mak0 - ta_ml_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
            close;

        end
    end
end
function plot_temp_polar_interp_mean_comp(type, par)
    make_dirs(type, par)

    % load data
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    grid_ml = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    ta_ml_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', type, par.lat_interp));
    ta_ml_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', type, par.lat_interp));
    ta_ml_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', type, par.lat_interp));
    ta_ml_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', type, par.lat_interp));

    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', 'echam_pl');
    grid_pl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'echam_pl'));
    ta_pl_lin = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_lin.mat', 'echam_pl', par.lat_interp));
    ta_pl_cub = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_cub.mat', 'echam_pl', par.lat_interp));
    ta_pl_spl = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_spl.mat', 'echam_pl', par.lat_interp));
    ta_pl_mak = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat_mak.mat', 'echam_pl', par.lat_interp));

    si = grid_pl.grid.dim3.si;

    lat_bound_list = [-85 -80 -70 70 80 85];
    mon_list = 1:12;

    for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
        for m = 1:length(mon_list); mon = mon_list(m)
            dlat = 0.25; % step size for standard lat grid
            if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole;
            else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; end;
            clat = cosd(lat); % cosine of latitude for cosine weighting

            folder = sprintf('%s/temp/mon_%g/poleward_of_lat_%g', par.plotdir, mon, lat_bound);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            ta_pl_lin0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_lin.tasi.lo(:,mon,:), lat));
            ta_ml_lin0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_lin.tasi.lo(:,mon,:), lat));
            ta_pl_cub0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_cub.tasi.lo(:,mon,:), lat));
            ta_ml_cub0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_cub.tasi.lo(:,mon,:), lat));
            ta_pl_spl0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_spl.tasi.lo(:,mon,:), lat));
            ta_ml_spl0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_spl.tasi.lo(:,mon,:), lat));
            ta_pl_mak0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_pl_mak.tasi.lo(:,mon,:), lat));
            ta_ml_mak0 = squeeze(interp1(grid_pl.grid.dim3.lat, ta_ml_mak.tasi.lo(:,mon,:), lat));

            clat_vert = repmat(clat', [1 size(ta_pl_lin0, 2)]);

            ta_pl_lin0 = nansum(ta_pl_lin0.*clat_vert, 1)/nansum(clat);
            ta_ml_lin0 = nansum(ta_ml_lin0.*clat_vert, 1)/nansum(clat);
            ta_pl_cub0 = nansum(ta_pl_cub0.*clat_vert, 1)/nansum(clat);
            ta_ml_cub0 = nansum(ta_ml_cub0.*clat_vert, 1)/nansum(clat);
            ta_pl_spl0 = nansum(ta_pl_spl0.*clat_vert, 1)/nansum(clat);
            ta_ml_spl0 = nansum(ta_ml_spl0.*clat_vert, 1)/nansum(clat);
            ta_pl_mak0 = nansum(ta_pl_mak0.*clat_vert, 1)/nansum(clat);
            ta_ml_mak0 = nansum(ta_ml_mak0.*clat_vert, 1)/nansum(clat);

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_lin0, si, '-');
            h_ml = plot(ta_ml_lin0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Linear', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_lin', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_cub0, si, '-');
            h_ml = plot(ta_ml_cub0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Piecewise Cubic Hermite', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_cub', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_spl0, si, '-');
            h_ml = plot(ta_ml_spl0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Cubic Spline', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_spl', folder), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_pl = plot(ta_pl_mak0, si, '-');
            h_ml = plot(ta_ml_mak0, si, '--');
            xlabel('T (K)');
            ylabel('$\sigma$ (unitless)');
            legend([h_pl h_ml], '$T_{\mathrm{pl}(\sigma)}$', '$T_{\mathrm{ml}}(\sigma)$', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$, Modified Akima', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/zon_comp_mak', folder), '-dpng', '-r300');
            close;

            % ALL PL
            figure(); clf; hold all; box on;
            plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-', 'linewidth', 1.5, 'color', 0.5*[1 1 1]);
            hlin=plot(ta_pl_lin0, grid_pl.grid.dim3.si, '-k');
            hcub=plot(ta_pl_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_pl_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_pl_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{pl}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_pl', folder), '-dpng', '-r300');
            close;

            % ALL ML
            figure(); clf; hold all; box on;
            hlin=plot(ta_ml_lin0, grid_ml.grid.dim3.si, '-k');
            hcub=plot(ta_ml_cub0, grid_ml.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_ml_spl0, grid_ml.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_ml_mak0, grid_ml.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{ml}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_ml', folder), '-dpng', '-r300');
            close;

            % ALL DIFF
            figure(); clf; hold all; box on;
            line([0 0], [0.1 1], 'color', 'k', 'linewidth', 0.5);
            hlin=plot(ta_pl_lin0 - ta_ml_lin0, grid_pl.grid.dim3.si, '-k');
            hcub=plot(ta_pl_cub0 - ta_ml_cub0, grid_pl.grid.dim3.si, '--', 'color', par.blue);
            hspl=plot(ta_pl_spl0 - ta_ml_spl0, grid_pl.grid.dim3.si, ':', 'color', par.orange);
            hmak=plot(ta_pl_mak0 - ta_ml_mak0, grid_pl.grid.dim3.si, '-.', 'color', par.maroon);
            xlabel('$T_{\mathrm{pl}}(\sigma)-T_{\mathrm{ml}}(\sigma)$ (K)');
            ylabel('$\sigma$ (unitless)');
            legend([hlin hcub hspl hmak], 'Linear', 'Piecewise Cubic Hermite', 'Cubic Spline', 'Modified Akima', 'location', 'eastoutside');
            title(sprintf('ECHAM, Month=%g, Lat=$%g^\\circ$ to $%g^\\circ$', mon, lat_bound, lat_pole));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
            print(sprintf('%s/all_diff', folder), '-dpng', '-r300');
            close;

        end
    end
end
