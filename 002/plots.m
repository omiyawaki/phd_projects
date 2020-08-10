clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.ep_swp = 0.5; %[0.25 0.3 0.35]; % threshold value for determining RCE
par.ga_swp = 0.9; % threshold for determining RAE
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget (option for ERA-Interim)
% if == stf, use SH and LH data from ERA-Interim to close energy budget (option for ERA-Interim)
% if == era5, use ERA5 radiative cooling and surface turbulent fluxes to close energy budget (only option for ERA5)
par.closure = 'era5';
par.do_surf = 0; % whether or not to calculate temperature profile in pressure grids including ts and ps data
par.cpd = 1005.7; par.Rd = 287; par.Rv = 461; par.g = 9.81; par.L = 2.501e6; par.a = 6357e3; par.eps = 0.622; % common constants, all in SI units for brevity
% set default figure parameters
if 1
    par.ppos = [0 0 10/3 7/3];
    par.ppos_sq = [0 0 10/3 10/3];
    par.ppos_wide = [0 0 13/3 7/3];
    par.ppos_verywide = [0 0 16/3 7/3];
    par.fs = 10;
    set(0, 'DefaultLineLineWidth', 1.1);
    set(0, 'DefaultFigureUnits', 'inches', 'DefaultFigurePosition', par.ppos_wide);
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    par.navy = 0.2*[0, 0.447, 0.741];
    par.darkblue = 0.5*[0, 0.447, 0.741];
    par.blue = [0, 0.447, 0.741];
    par.orange = [0.85, 0.325, 0.098];
    par.yellow = [0.929, 0.694, 0.125];
    par.purple = [0.494, 0.184, 0.556];
    par.green = [0.466, 0.674, 0.188];
    par.cyan = [0.301, 0.745, 0.933];
    par.maroon = [0.635, 0.078, 0.184];
    par.brown = 0.5*[0.635, 0.078, 0.184];
    par.darkbrown = 0.2*[0.635, 0.078, 0.184];
    par.purple = 0.5*[0.4940, 0.1840, 0.5560];
    par.gray = 0.5*[1 1 1];
end
gcm_info

%% call functions
% plot_rad_lat(par)
% plot_rad_lon_lat(par)
% plot_tediv_lat(par)

type = 'era5';
% choose_plots(type, par);
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % choose_plots(type, par);
end

% sweep through various threshold values
for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
    type = 'era5';
    % choose_plots_ep(type, par)
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        choose_plots_ep(type, par)
    end
end

%% define functions
function choose_plots(type, par)
    % plot_va(type, par) % plot lat x height structure of meridional velocity
    % plot_temp_zon(type, par) % plot temperature profiles at specific latitudes
    % plot_ma_diff(type, par) % plot difference of temperature profile from moist adiabat
    plot_flux(type, par) % plot various energy fluxes in mon x lat and lon x lat space
end % select which functions to run at a time
function plot_va(type, par)
    make_dirs(type, par)

    % load data
    if strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        lat = grid.dim3.lat;
        plev = grid.dim3.plev/100;
        var = 'va';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/gcm/%s/%s_Amon_%s_piControl_r1i1p1_*.ymonmean.nc', par.model, var, par.model));
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
function plot_temp_zon(type, par)
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_lon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.lat_interp));
        plev = grid.dim3.plev/100;
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            ta_z.(land).(time) = squeeze(nanmean(ta_t.(land).(time), 1)); % zonal average
            tasi_z.(land).(time) = squeeze(nanmean(tasi_t.(land).(time), 1)); % zonal average
            if par.do_surf; tai_z.(land).(time) = squeeze(nanmean(tai_t.(land).(time), 1)); end % zonal average
            tsurf_z.(land).(time) = squeeze(nanmean(tsurf_t.(land).(time), 1)); % zonal average
            psurf_z.(land).(time) = squeeze(nanmean(psurf_t.(land).(time), 1)); % zonal average
            zsurf_z.(land).(time) = squeeze(nanmean(zsurf_t.(land).(time), 1)); % zonal average
            ta_sp = ta_z.(land).(time)(7,:); % sounding at -88.5 S
            ta_np = ta_z.(land).(time)(end-6,:); % sounding at 88.5 N
            ta_eq = interp1(lat, ta_z.(land).(time), 0); % sounding at equator

            z_int=[0:500:20e3]';
            taz_z.(land).(time) = squeeze(nanmean(taz_t.(land).(time), 1)); % zonal average
            taz_sp = taz_z.(land).(time)(7,:); % sounding at -88.5 S
            taz_np = taz_z.(land).(time)(end-6,:); % sounding at 88.5 N
            taz_eq = interp1(lat, taz_z.(land).(time), 0); % sounding at equator

            idx_f = min(find(~isnan(ta_z.(land).(time)(:,floor(size(ta_z.(land).(time),2)/2))))); % first non-nan latitude index
            [~, idx_fi] = min(abs(lat - ceil(lat(idx_f)))); % first non-nan integer latitude index
            step = 8; % step size of adjacent latitudes
            nums = 9; % number of latitudes to plot at once

            if par.do_surf; v_vec = {'p', 'z', 'si', 'pi'};
            else v_vec = {'p', 'z', 'si'}; end
            for v = v_vec; vert = v{1};
                for h = {'sh', 'nh', 'eq', 'nmid', 'smid'}; hemi = h{1};
                    figure(); clf; hold all;
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
                    end
                    c = 1;
                    for i = indices;
                        if strcmp(vert, 'p'); h{c} = plot(ta_z.(land).(time)(i,:), plev, 'color', (c/(1.3*nums))*[1 1 1]);
                        elseif strcmp(vert, 'z'); h{c} = plot(taz_z.(land).(time)(i,:), z_int*10^-3, 'color', (c/(1.3*nums))*[1 1 1]);
                        elseif strcmp(vert, 'si'); h{c} = plot(tasi_z.(land).(time)(i,:), grid.dim3.si, 'color', (c/(1.3*nums))*[1 1 1]);
                        elseif strcmp(vert, 'pi'); h{c} = plot(tai_z.(land).(time)(i,:), par.pa/1e2, 'color', (c/(1.3*nums))*[1 1 1]); end;
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;
                end % hemi
            end % vert

            figure(); clf; hold all;
            h_eq = plot(ta_eq, plev, 'color', par.orange);
            h_np = plot(ta_np, plev, 'color', par.blue);
            h_sp = plot(ta_sp, plev, '--', 'color', par.blue);
            xlabel('T (K)'); ylabel('p (hPa)');
            legend([h_eq h_np h_sp], 'Equator', '88.5', '$-88.5$', 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon/%s/%s/p/all', par.plotdir, land, time), '-dpng', '-r300');
            close;

            figure(); clf; hold all;
            h_eq = plot(taz_eq, z_int*10^-3, 'color', par.orange);
            h_np = plot(taz_np, z_int*10^-3, 'color', par.blue);
            h_sp = plot(taz_sp, z_int*10^-3, '--', 'color', par.blue);
            xlabel('T (K)'); ylabel('z (km)');
            legend([h_eq h_np h_sp], 'Equator', '88.5', '-$88.5$', 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on')
            print(sprintf('%s/temp_zon/%s/%s/z/all', par.plotdir, land, time), '-dpng', '-r300');
            close;
        end
    end
end
function plot_ma_diff(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ta_mon_lat.mat', type, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ta_mon_lat.mat', type, par.model, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ma_mon_lat.mat', type, par.model, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100;
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        for plev_eval = [300:100:500] % evaluation plev in hPa
            diff = permute(ta.(land) - ma.(land).ta, [3 1 2]); % bring plev to front
            diff = squeeze(interp1(plev, diff, plev_eval)); % evaluate difference at plev_eval
            % lat x lon of RCE and RAE
            figure(); clf; hold all;
            cmp = colCog(40);
            colormap(cmp);
            imagesc([1:12], [lat(1) lat(end)], diff);
            caxis([-20 20]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-20 20], 'ytick', [-20:5:20], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$(T - T_m)_{%g \\,\\mathrm{hPa}}$ (K)', plev_eval));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ma_diff/plev_%g/%s/ma_diff_lat_lon', par.plotdir, plev_eval, land), '-dpng', '-r300');
            close;
        end
    end
end
function plot_flux(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
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

        if any(strcmp(type, {'era5', 'erai'})); f_vec = {'mse', 'dse', 'db13', 'db13s'};
        elseif strcmp(type, 'gcm'); f_vec = {'mse', 'dse', 'mse2'}; end
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

            % lat x mon dependence of RCE and RAE
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
end % function

function choose_plots_ep(type, par)
    plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_rcae(type, par) % plot RCE/RAE regimes
    % plot_rcae_rc(type, par) % plot RCE/RAE regimes with R1 recomputed at the very end (order of operations test)
    % plot_temp(type, par) % plot temperature profiles
end % select which ep-functions to run at a time
function plot_energy_lat(type, par) % latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, par] = load_flux(type, par); % load data
    make_dirs(type, par)

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        if any(strcmp(type, {'erai'})); f_vec = {'mse', 'dse', 'db13', 'db13s'};
        elseif any(strcmp(type, {'era5'})); f_vec = {'mse', 'dse', 'db13', 'db13s', 'div'};
        elseif strcmp(type, 'gcm'); f_vec = {'mse', 'dse', 'mse2'}; end
        for f = f_vec; fw = f{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

                if strcmp(fw, 'mse2')
                    figure(); clf; hold all;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux.(land).(time).lw, 'color', par.gray); text(0, 0.85*interp1(lat,flux.(land).(time).lw,0), '\boldmath{$\mathrm{LW}$}', 'color', par.gray);
                    plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        if contains(fw, 'db')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif strcmp(fw, 'div')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue); end
                            plot(lat, -flux.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'gcm')
                        if contains(fw, 'mse'); plot(lat, flux.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    if any(strcmp(type, {'era5', 'erai'}));
                        if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                        else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end;
                    elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text)); end
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf 150]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', par.plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                else
                    figure(); clf; hold all;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).res.dse, '--', 'color', par.maroon); text(-30, 2*interp1(lat,flux.(land).(time).res.dse,-30), '\boldmath{$\nabla\cdot F_s$}', 'color', par.maroon);
                    else; plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m$}', 'color', par.maroon); end;
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        if contains(fw, 'db')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif strcmp(fw, 'div')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue);
                            end
                            plot(lat, -flux.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'gcm')
                        if strcmp(fw, 'mse'); plot(lat, flux.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    if any(strcmp(type, {'era5', 'erai'}));
                        if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                        else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end;
                    elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                    end
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf 150]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', par.plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                end

                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                if any(strcmp(fw, {'mse', 'mse2', 'db13', 'db13s'})); plot(lat,flux.(land).(time).r1.(fw), '-k');
                elseif strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                end
                if any(strcmp(type, {'era5', 'erai'}));
                    if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                    else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                end
                xlabel('latitude (deg)');
                if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                else ylabel('$R_1$ (unitless)'); end
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                print(sprintf('%s/energy-flux/%s/%s/%s-r1', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;
            % northward M/DSE transport
                figure(); clf; hold all;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([0 0], [min(vh.(land).(time).(fw)) max(vh.(land).(time).(fw))]*10^-15, 'linewidth', 0.5, 'color', 'k');
                if any(strcmp(fw, {'mse', 'db13', 'db13s'})); plot(lat, vh.(land).(time).(fw)*10^-15, 'color', par.maroon);
                elseif strcmp(fw, 'dse'); plot(lat, vh.(land).(time).(fw)*10^-15, '--', 'color', par.maroon);
                end
                xlabel('latitude (deg)'); ylabel('PW')
                if strcmp(fw, 'db13s'); title(sprintf('Northward %s Transport, %s', 'DB13*', upper(time)));
                else title(sprintf('Northward %s Transport, %s', upper(fw), upper(time))); end;
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                if strcmp(fw, 'dse'); set(gca, 'ytick', [-5:5]); end;
                print(sprintf('%s/transport/%s/%s/%s', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;
                if ~contains(fw, {'db'})
                % MSE/DSE transport plotted together
                    figure(); clf; hold all;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    line([0 0], [min([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse]) max([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse])]*10^-15, 'linewidth', 0.5, 'color', 'k');
                    h_mse = plot(lat, vh.(land).(time).mse*10^-15, 'color', par.maroon);
                    h_dse = plot(lat, vh.(land).(time).dse*10^-15, '--', 'color', par.maroon);
                    h_lh = plot(lat, (vh.(land).(time).mse-vh.(land).(time).dse)*10^-15, ':', 'color', par.maroon);
                    legend([h_mse h_dse h_lh], '$F_m$', '$F_s$', '$F_{m}-F_{s}$', 'location', 'eastoutside');
                    xlabel('latitude (deg)'); ylabel('PW')
                    title(sprintf('Northward Energy Transport, %s', upper(time)));
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    print(sprintf('%s/transport/%s/%s/all', par.plotdir, land, time), '-dpng', '-r300');
                    close;
                end
            end % land
        end % end mse/dse loop

    end % time

    for f = f_vec; fw = f{1};
        for l = {'lo', 'l', 'o'}; land = l{1};
            % northward M/DSE transport, mon x lat
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
            figure(); clf; hold all;
            cmp = colCog(30);
            colormap(cmp);
            [C, h] = contour(mesh_lat, mesh_mon, vh_mon.(land).(fw)*10^-15, -5:1:5);
            clabel(C, h, [-4:2:4], 'fontsize', 6, 'interpreter', 'latex');
            caxis([-5 5]);
            xlabel('Month'); ylabel('latitude (deg)');
            title(sprintf('Northward %s Transport (PW)', upper(fw)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'fontsize', par.fs, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on')
            print(sprintf('%s/transport/%s/all/%s', par.plotdir, land, fw), '-dpng', '-r300');
            close;
        end
    end

end
function plot_rcae(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
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

        if any(strcmp(type, {'erai'})); f_vec = {'mse', 'dse', 'mse2', 'db13', 'db13s'};
        elseif any(strcmp(type, {'era5'})); f_vec = {'mse', 'dse', 'mse2', 'db13', 'db13s', 'div'};
        elseif strcmp(type, 'gcm'); f_vec = {'mse', 'dse', 'mse2'}; end
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
function plot_rcae_rc(type, par)
    [~, ~, ~, lat, par] = load_flux(type, par); % load data
    make_dirs_ep(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
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

        if any(strcmp(type, {'erai'})); f_vec = {'mse', 'dse', 'db13', 'db13s'};
        elseif any(strcmp(type, {'era5'})); f_vec = {'mse', 'dse', 'db13', 'db13s', 'div'};
        elseif strcmp(type, 'gcm'); f_vec = {'mse', 'dse', 'mse2'}; end
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
function plot_temp(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ta.mat', type, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ma.mat', type, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ta.mat', type, par.model, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ma.mat', type, par.model, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    end
    for f = {'mse', 'dse'}; fw = f{1};
        for c = fieldnames(ta.rce.tp.(fw))'; crit = c{1};
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                for l = {'lo', 'l', 'o'}; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

                % RCE and RAE separated into NH and SH
                    figure(); clf; hold all;
                    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
                    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
                    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
                    h_rae_nh = plot(ta.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
                    h_rae_sh = plot(ta.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('%s, %s', upper(time), land_text));
                    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE', 'NH ML RCE', 'SH ML RCE', 'NH RAE', 'SH RAE', 'location', 'northeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    line([0 0], [100 1000], 'linewidth', 0.5, 'color', 'k');
                    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time)-ma.rce.tp.(fw).(crit).(land).(time).ta, plev, 'color', par.maroon);
                    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time)-ma.rce.nh.(fw).(crit).(land).(time).ta, plev, 'color', par.orange);
                    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time)-ma.rce.sh.(fw).(crit).(land).(time).ta, plev, '--', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [200 1000], 'xtick', [-5:5:40], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_diff', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % Tropical RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
                    h_rce_ma = plot(ma.rce.tp.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.maroon);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_tp', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % NH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, 'color', par.orange);
                    h_rce_ma = plot(ma.rce.nh.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('NH ML RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;
                % SH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, 'color', par.orange);
                    h_rce_ma = plot(ma.rce.sh.(fw).(crit).(land).(time).ta, plev, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('p (hPa)');
                    title(sprintf('SH ML RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;
                end % land
                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    h_rce_l = plot(ta.rce.nh.(fw).(crit).l.(time)-ma.rce.nh.(fw).(crit).l.(time).ta, plev, 'color', par.orange);
                    h_rce_o = plot(ta.rce.nh.(fw).(crit).o.(time)-ma.rce.nh.(fw).(crit).o.(time).ta, plev, ':', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('p (hPa)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    legend([h_rce_l, h_rce_o], 'Land', 'Ocean', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [100 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_lo_diff', par.plotdir, par.ep, par.ga, fw, crit, 'lo', time), '-dpng', '-r300');
                    close;
            end % time avg
        end % RCE/RAE definition
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
    figure(); clf; hold all;
    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
    h_rce_nh = plot(ta.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
    h_rce_sh = plot(ta.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
    h_rae_nh = plot(ta.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
    h_rae_sh = plot(ta.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
    xlabel('T (K)'); ylabel('p (hPa)');
    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE ($\pm 30^\circ$)', 'NH ML RCE ($>+30^\circ$)', 'SH ML RCE ($<-30^\circ$)', 'NH RAE', 'SH RAE', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp', par.plotdir), '-dpng', '-r300');
    close;
    % Legend for moist adiabat comparisons
    figure(); clf; hold all;
    h_rce_tp = plot(ta.rce.tp.(fw).(crit).(land).(time), plev, '-k');
    h_rce_tp_ma = plot(ma.rce.tp.(fw).(crit).(land).(time).ta, plev, ':k');
    xlabel('T (K)'); ylabel('p (hPa)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        legend([h_rce_tp h_rce_tp_ma], upper(type), 'Moist adiabat', 'location', 'eastoutside');
    elseif strcmp(type, 'gcm')
        legend([h_rce_tp h_rce_tp_ma], par.model, 'Moist adiabat', 'location', 'eastoutside');
    end
    title(upper(sprintf('%s', time)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp_ma', par.plotdir), '-dpng', '-r300');
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
    gcm=load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/gcm/%s/std/flux_zt.mat', par.model));
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

function [flux_zt, vh, vh_mon, lat, par] = load_flux(type, par)
    % load processed data/proc
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh_mon.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/flux_zt.mat', type, par.model, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh.mat', type, par.model, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh_mon.mat', type, par.model, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
    end
end
function make_dirs(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        for plev_eval = [300:100:500]
            if ~exist(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land), 'dir')
                mkdir(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land));
            end
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time));
            end
            if ~exist(sprintf('%s/transport/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/transport/%s/%s', par.plotdir, land, time));
            end
            if par.do_surf; v_vec = {'p', 'z', 'si', 'pi'};
            else v_vec = {'p', 'z', 'si'}; end
            for v = v_vec; vert = v{1};
                if ~exist(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert));
                end
            end
            if strcmp(type, 'erai')
                for f = {'mse', 'dse', 'mse2', 'db13', 'db13s'}; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            elseif strcmp(type, 'era5')
                for f = {'mse', 'dse', 'mse2', 'db13', 'db13s', 'div'}; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            elseif strcmp(type, 'gcm')
                for f = {'mse', 'dse', 'mse2'}; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            end
        end
    end

    if ~exist(sprintf('%s/va', par.plotdir), 'dir')
        mkdir(sprintf('%s/va', par.plotdir));
    end

    if ~exist(sprintf('%s/legends', par.plotdir), 'dir')
        mkdir(sprintf('%s/legends', par.plotdir));
    end
end
function make_dirs_ep(type, par)
    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        if any(strcmp(type, {'erai'})); f_vec = {'mse', 'dse', 'mse2', 'db13', 'db13s'};
        elseif any(strcmp(type, {'era5'})); f_vec = {'mse', 'dse', 'mse2', 'db13', 'db13s', 'div'};
        elseif strcmp(type, 'gcm'); f_vec = {'mse', 'dse', 'mse2'}; end
        for f = f_vec; fw = f{1};
            for j = {'def', 'jak', 'jak30', 'jak10', 'pe', 'cp', 'w500', 'vh2', 'vh3', 'vh4'}; crit = j{1};
                for l = {'lo', 'l', 'o'}; land = l{1};
                    for k = {'ann', 'djf', 'jja', 'mam', 'son'}; time = k{1};
                        if ~exist(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp', par.plotdir, par.ep, par.ga, fw, crit, land, time), 'dir')
                            mkdir(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp', par.plotdir, par.ep, par.ga, fw, crit, land, time));
                        end
                    end
                end
            end
        end
    end
end
