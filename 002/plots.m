clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.gcm.clim = 'piControl'; % choose either piControl or abrupt4xCO2
par.echam.clim = '20170908'; % choose from 20170908 (snowball) or 20170915_2 (modern)
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold value for determining RCE
par.ga_swp = 0.9; % threshold for determining RAE
par.si_eval = [0.8 0.85 0.9]; % sigma level for calculating inversion strength
par.pa_eval = 500e2; % pressure level for calculating ma_diff
par.si_bl = 0.85; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.si_up = 0.4; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
par.ta_thresh = 6.5; % criteria for temperature difference in K of plotting contour for closeness to moist adiabat
par.ga_thresh = 10; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.ga_bl_thresh = 90; % criteria for lapse rate difference in % of plotting contour for closeness to moist adiabat
par.inv_thresh = -4; % criteria for temperature difference in K of plotting contour for inversion strength
par.albedo_thresh = 0.6; % criteria for surface albedo of plotting contour for inversion strength
par.sn_thresh = 1e-1; % criteria for surface albedo of plotting contour for inversion strength
par.alb_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
par.albcs_thresh = 0.8; % criteria for surface albedo of plotting contour for inversion strength
% par.era.fw = {'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt', 'div79'};
% par.era.fw = {'div79', 'mse', 'dse', 'db13', 'db13s', 'db13t', 'div', 'divt'};
par.era.fw = {'div00'};
par.gcm.fw = {'mse', 'dse'};
par.echam.fw = {'mse', 'dse'};
par.pa = linspace(1000,10,100)*1e2; % high resolution vertical grid to interpolate to
par.z = [0:500:40e3]';
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
    par.ppos_larger = [0 0 16/3 13/3];
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

type = 'echam_ml';
% choose_plots(type, par);
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % choose_plots(type, par);
end

% sweep through various threshold values
for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
    type = 'erai';
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
    % plot_ga_diff(type, par) % plot difference of temperature profile from moist adiabat
    % plot_ga_si_diff(type, par) % plot difference of temperature profile from moist adiabat in sigma
    % plot_ga_malr_si_diff(type, par) % plot difference of temperature profile from moist adiabat in sigma
    % plot_ga_z_diff(type, par) % plot difference of temperature profile from moist adiabat in z coord
    % plot_ga_diff_mon_lat(type, par) % plot ga diff
    plot_ga_malr_diff_mon_lat(type, par) % plot ga diff
    % plot_ga_malr_diff_lon_lat(type, par) % plot ga diff
    % plot_inv_str(type, par) % plot strength of inversion
    % plot_inv_str_alt(type, par) % plot strength of inversion
    % plot_flux(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_flux_block(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_dr1(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dr1_decomp4(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dr1_alt(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dr2(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_tend_comp(type, par) % plot MSE tendency of DB13 with that calculated from ERA-Interim
    % plot_trop(type, par) % plot WMO tropopause
    % plot_alb(type, par) % plot albedo
    % plot_sn(type, par) % plot snow depth
end % select which functions to run at a time
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
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_lon_lat.mat', type, par.echam.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        plev = 1:47;
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
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

            z_int=par.z;
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
                for h = {'sh', 'nh', 'eq', 'nmid', 'smid', 'nhall', 'shall'}; hemi = h{1};
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
                    elseif strcmp(h, 'nhall');
                        if mod(length(lat),2)==0; indices = length(lat)-idx_fi+1:-16:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = length(lat)-idx_fi+1:-16:floor(length(lat)/2); end;
                    elseif strcmp(h, 'shall');
                        if mod(length(lat),2)==0; indices = idx_fi:16:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = idx_fi:16:floor(length(lat)/2); end;
                    end
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            if strcmp(vert, 'p'); h{c} = plot(ta_z.(land).(time)(i,:), plev);
                            elseif strcmp(vert, 'z'); h{c} = plot(taz_z.(land).(time)(i,:), z_int*10^-3);
                            elseif strcmp(vert, 'si'); h{c} = plot(tasi_z.(land).(time)(i,:), grid.dim3.si);
                            elseif strcmp(vert, 'pi'); h{c} = plot(tai_z.(land).(time)(i,:), par.pa/1e2); end;
                        else
                            if strcmp(vert, 'p'); h{c} = plot(ta_z.(land).(time)(i,:), plev, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'z'); h{c} = plot(taz_z.(land).(time)(i,:), z_int*10^-3, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'si'); h{c} = plot(tasi_z.(land).(time)(i,:), grid.dim3.si, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'pi'); h{c} = plot(tai_z.(land).(time)(i,:), par.pa/1e2, 'color', (c/(1.3*nums))*[1 1 1]); end;
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
    make_dirs(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
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
            if strcmp(type, 'era5') | strcmp(type, 'erai')
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
function plot_ga_malr_si_diff(type, par) % sigma coordinates
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
    % load(sprintf('%s/dtdzzsi.mat', prefix)); % read model lapse rate
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read model lapse rate
    % load(sprintf('%s/malrzsi.mat', prefix)); % read moist adiabatic lapse rate
    load(sprintf('%s/malrsi.mat', prefix)); dtmdzzsi = dtmdzsi; clear dtmdzsi; % read moist adiabatic lapse rate

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
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JANUARY
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jan', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JULY
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jul', -10:1:10);
    % clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'log', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_jul', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc_jan', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_sc_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtmdz_zt', -10:1:10, 'w');
    clabel(C, h, [-50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
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
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma', par.plotdir), '-dpng', '-r300');
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
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
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
function plot_ga_malr_diff_mon_lat(type, par) % plot inversion strength
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
    load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_dalr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp));

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % mon x lat of diff
        figure(); clf; hold all;
        cmp = colCog(12);
        colormap(cmp);
        % imagesc([1:12], [lat(1) lat(end)], ga_malr_diff.(land));
        contourf(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_malr_diff.(land), -100:10:100, '-k');
        clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex');
        caxis([-60 60]);
        xlabel('Month'); ylabel('Latitude (deg)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, %s', upper(type), land_text));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, %s', par.model, land_text));
        end
        cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/%s/ga_malr_diff_mon_lat', par.plotdir, land), '-dpng', '-r300');
        close;

        % mon x lat of diff
        figure(); clf; hold all;
        cmp = colCog(12);
        colormap(cmp);
        % imagesc([1:12], [lat(1) lat(end)], ga_dalr_diff.(land));
        contourf(mesh_lat, mesh_mon, ga_dalr_diff.(land), -100:10:100, 'linecolor', 'none')
        caxis([-60 60]);
        xlabel('Month'); ylabel('Latitude (deg)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, %s', upper(type), land_text));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, %s', par.model, land_text));
        end
        cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/%s/ga_dalr_diff_mon_lat', par.plotdir, land), '-dpng', '-r300');
        close;

        % DALR BL mon x lat of diff
        figure(); clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        % imagesc([1:12], [lat(1) lat(end)], ga_dalr_bl_diff.(land));
        contourf(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, 'linecolor', 'none');
        [C,h]=contour(mesh_lat, mesh_mon, ga_dalr_bl_diff.(land), -100:10:100, '-w');
        clabel(C, h, -100:20:100, 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
        caxis([-100 100]);
        xlabel('Month'); ylabel('Latitude (deg)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, %s', upper(type), land_text));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, %s', par.model, land_text));
        end
        cb = colorbar('limits', [0 100], 'ytick', [0:10:100], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, sprintf('$\\left[(\\Gamma_d - \\Gamma)/\\Gamma_d\\right]_{1.0-%g}$ (\\%%)', par.si_bl));
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/ga_malr_diff/%s/ga_dalr_bl_diff_mon_lat', par.plotdir, land), '-dpng', '-r300');
        close;

    end
end
function plot_ga_malr_diff_lon_lat(type, par) % plot inversion strength
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
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/%s/ga_malr_diff_si_lon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_dalr_diff_si_lon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_dalr_bl_diff_si_lon_lat.mat', prefix_proc, par.lat_interp));

    [mesh_lat, mesh_lon] = meshgrid(grid.dim3.lon, lat);

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
            [C,h]=contour(mesh_lat, mesh_lon, ga_malr_diff_t.(land).(time)', -100:10:100, '-k');
            clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex');
            caxis([-60 60]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-40 60], 'ytick', [-40:10:60], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{%g-%g}$ (\\%%)', par.si_bl, par.si_up));
            set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ga_malr_diff/%s/%s/ga_malr_diff_lon_lat', par.plotdir, land, time), '-dpng', '-r300');
            close;

            % DALR BL lon x lat of diff
            figure(); clf; hold all;
            cmp = colCog(12);
            colormap(cmp);
            contourf(mesh_lat, mesh_lon, ga_dalr_bl_diff_t.(land).(time)', -100:10:100, 'linecolor', 'none');
            [C,h]=contour(mesh_lat, mesh_lon, ga_dalr_bl_diff_t.(land).(time)', -100:10:100, '-k');
            clabel(C, h, -100:10:100, 'fontsize', 6, 'interpreter', 'latex');
            caxis([-100 100]);
            xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
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
function plot_flux(type, par)
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
    load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_malr_diff_orig = ga_malr_diff; clear ga_diff;
    load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_dalr_bl_diff_orig = ga_dalr_bl_diff; clear ga_bl_diff;
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        ga_malr_diff = ga_malr_diff_orig.(land);
        ga_dalr_bl_diff = ga_dalr_bl_diff_orig.(land);

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
        for f = f_vec; fw = f{1};
            % lat x mon dependence of RCE and RAE
            if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
            else var_text = '$\nabla \cdot F_m$'; end
            figure(); clf; hold all;
            cmp = colCog(12);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw), -150:25:150, 'linecolor', 'none');
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, var_text));
            elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, var_text)); end;
            caxis([-150 150]);
            cb = colorbar('limits', [-150 150], 'ytick', [-150:25:150], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_div_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % R1 lat x mon dependence of RCE and RAE
            var_text = '$R_1$';
            figure(); clf; hold all;
            cmp = colCog(10);
            colormap(flipud(cmp));
            contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'none');
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, var_text));
            elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, var_text)); end;
            caxis([-1 1]);
            cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (unitless)', var_text));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            close;

            % GA MALR and GA_BL DALR contour lat x mon dependence of RCE and RAE
            figure(); clf; hold all;
            cmp = colCog(10);
            colormap(flipud(cmp));
            contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'none');
            contour(mesh_lat, mesh_mon, abs(ga_malr_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
            contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
            elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
            caxis([-1 1]);
            cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (unitless)', var_text));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_r1_mon_lat_ga_malr_ga_dalr_bl_overlay', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            if strcmp(fw, 'mse2')
                % lat x mon dependence of RCE and RAE
                var_text = 'LW';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, flux_z.(land).lw, 'linecolor', 'none');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, var_text));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, var_text)); end;
                caxis([-250 250]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_lw_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
                close;
            else
                % lat x mon dependence of RCE and RAE
                var_text = '$R_a$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                contourf(mesh_lat, mesh_mon, flux_z.(land).ra.(fw), 'linecolor', 'none');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, var_text));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, var_text)); end;
                caxis([-150 150]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/0_ra_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
                close;
            end

            [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim3.lon, lat);
            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                % lat x lon of RCE and RAE
                if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
                else var_text = '$\nabla \cdot F_m$'; end
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).res.(fw)', 'linecolor', 'none');
                caxis([-150 150]);
                contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/div_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                var_text = '$R_1$';
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(flipud(cmp));
                contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).r1.(fw)', [-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'none');
                contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                caxis([-1 1]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s',upper(type), var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/r1_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                if strcmp(fw, 'mse2')
                    var_text = 'LW';
                    figure(); clf; hold all;
                    cmp = colCog(20);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).lw', 'linecolor', 'none');
                    contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    caxis([-250 250]);
                    if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s',upper(type), var_text, upper(time), land_text)); end;
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/lw_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                    close;
                else
                    var_text = '$R_a$';
                    figure(); clf; hold all;
                    cmp = colCog(20);
                    colormap(cmp);
                    contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).ra.(fw)', 'linecolor', 'none');
                    contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                    caxis([-150 150]);
                    if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                    elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s',upper(type), var_text, upper(time), land_text)); end;
                    xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                    set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/flux/%s/%s/%s/ra_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                    close;
                end

            end % for time
        end % for mse dse
    end % for land
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
            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], flux_z.(land).r1.(fw) - repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\Delta R_1$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % FRACTIONAL DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta R_1}{R_1}$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], (flux_z.(land).r1.(fw) - repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]))./flux_z.(land).r1.(fw));
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta R_1}{R_1}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
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
            imagesc([1 12], [lat(1) lat(end)], comp1);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_m)}{R_a}$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
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
            imagesc([1 12], [lat(1) lat(end)], comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$\frac{\Delta (\nabla\cdot F_m)}{R_a}-\frac{\nabla\cdot F_m}{R_a^2}\Delta R_a$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1+comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\Delta (\\nabla\\cdot F_m)}{R_a}-\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta R_a$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_c1+c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
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
            % DELTA R1 ALT COMP1 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{(\mathrm{SH+LH})}{R_a^2}\Delta (\nabla\cdot F_m)$';
            delta_fm = flux_z.(land).res.(fw) - repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
            ann_stf = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
            ann_ra = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
            comp1 = -ann_stf./(ann_ra).^2.*delta_fm;
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta (\\nabla\\cdot F_m)$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
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
            imagesc([1 12], [lat(1) lat(end)], comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta (\\mathrm{SH+LH})$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DELTA R1 SUM OF COMP1 and COMP2 lat x mon dependence of RCE and RAE
            var_text = '$-\frac{\mathrm{SH+LH}}{R_a^2}\Delta (\nabla\cdot F_m)+\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
            figure(); clf; hold all;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], comp1+comp2);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text)); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$-\\frac{\\mathrm{SH+LH}}{R_a^2}\\Delta (\\nabla\\cdot F_m)+\\frac{\\nabla\\cdot F_m}{R_a^2}\\Delta (\\mathrm{SH+LH})$ (unitless)'));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/dr1/%s/%s/0_dr1_alt_c1+c2_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

        end % for mse dse
    end % for land
end
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
function plot_sn(type, par) % plot snow depth
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

function choose_plots_ep(type, par)
    % plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_rcae(type, par) % plot RCE/RAE regimes
    % plot_rcae_rc(type, par) % plot RCE/RAE regimes with R1 recomputed at the very end (order of operations test)
    % plot_rcae_alt(type, par) % plot RCE/RAE regimes
    % plot_rcae_alt_rc(type, par) % plot RCE/RAE regimes with R1 recomputed at the very end (order of operations test)
    plot_rcae_alt_overlay(type, par) % plot RCE/RAE regimes
    % plot_rcae_alt_rc_overlay(type, par) % plot RCE/RAE regimes with R1 recomputed at the very end (order of operations test)
    % plot_temp(type, par) % plot temperature profiles
end % select which ep-functions to run at a time
function plot_energy_lat(type, par) % latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, par] = load_flux(type, par); % load data
    make_dirs(type, par)

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
        for f = f_vec; fw = f{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

                if strcmp(fw, 'mse2')
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux.(land).(time).lw, 'color', par.gray); text(0, 0.85*interp1(lat,flux.(land).(time).lw,0), '\boldmath{$\mathrm{LW}$}', 'color', par.gray);
                    if contains(fw, 'db')
                        plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    else
                        plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    end
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
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    if any(strcmp(type, {'era5', 'erai'}));
                        if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                        else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end;
                    elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                    elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', par.plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                else
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).res.dse, '--', 'color', par.maroon); text(-30, 2*interp1(lat,flux.(land).(time).res.dse,-30), '\boldmath{$\nabla\cdot F_s$}', 'color', par.maroon);
                    else; plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m$}', 'color', par.maroon); end;
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        if contains(fw, 'db')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif contains(fw, 'div')
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
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    if any(strcmp(type, {'era5', 'erai'}));
                        if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                        else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end;
                    elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                    elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', par.plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                end

                figure(); clf; hold all; box on;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                else; plot(lat,flux.(land).(time).r1.(fw), '-k'); end
                if any(strcmp(type, {'era5', 'erai'}));
                    if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                    else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
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
                figure(); clf; hold all; box on;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([0 0], [min(vh.(land).(time).(fw)) max(vh.(land).(time).(fw))]*10^-15, 'linewidth', 0.5, 'color', 'k');
                if strcmp(fw, 'dse'); plot(lat, vh.(land).(time).(fw)*10^-15, '--', 'color', par.maroon);
                else; plot(lat, vh.(land).(time).(fw)*10^-15, 'color', par.maroon);
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
                % if ~contains(fw, {'db'})
                % % MSE/DSE transport plotted together
                %     figure(); clf; hold all;
                %     line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                %     line([0 0], [min([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse]) max([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse])]*10^-15, 'linewidth', 0.5, 'color', 'k');
                %     h_mse = plot(lat, vh.(land).(time).mse*10^-15, 'color', par.maroon);
                %     h_dse = plot(lat, vh.(land).(time).dse*10^-15, '--', 'color', par.maroon);
                %     h_lh = plot(lat, (vh.(land).(time).mse-vh.(land).(time).dse)*10^-15, ':', 'color', par.maroon);
                %     legend([h_mse h_dse h_lh], '$F_m$', '$F_s$', '$F_{m}-F_{s}$', 'location', 'eastoutside');
                %     xlabel('latitude (deg)'); ylabel('PW')
                %     title(sprintf('Northward Energy Transport, %s', upper(time)));
                %     axis('tight');
                %     set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                %     set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                %     print(sprintf('%s/transport/%s/%s/all', par.plotdir, land, time), '-dpng', '-r300');
                %     close;
                % end
            end % land
        end % end mse/dse loop

    end % time

    for f = f_vec; fw = f{1};
        for l = {'lo', 'l', 'o'}; land = l{1};
            % northward M/DSE transport, mon x lat
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
            figure(); clf; hold all; box on;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], vh_mon.(land).(fw)*1e-15);
            cb = colorbar('ticks', [-5:1:5], 'ticklabelinterpreter', 'latex');
            ylabel(cb, '$F_m$ (PW)', 'interpreter', 'latex');
            % [C, h] = contour(mesh_lat, mesh_mon, vh_mon.(land).(fw)*10^-15, -5:1:5);
            % clabel(C, h, [-4:2:4], 'fontsize', 6, 'interpreter', 'latex');
            caxis([-5 5]);
            xlabel('Month'); ylabel('latitude (deg)');
            title(sprintf('Northward %s Transport (PW)', upper(fw)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
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
    load(sprintf('%s/%s/inv_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ga_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_diff_orig = ga_diff; clear ga_diff;
    load(sprintf('%s/%s/ga_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_bl_diff_orig = ga_bl_diff; clear ga_bl_diff;
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

        % % moist adiabat differences
        % ma_diff0 = ta.(land) - ma.(land).ta;
        % ma_diff = nan([size(ma_diff0,1) 14 size(ma_diff0,3)]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        % ma_diff(:,2:13,:) = ma_diff0;
        % ma_diff(:,1,:) = 1/2*(ma_diff0(:,1,:)+ma_diff0(:,12,:));
        % ma_diff(:,14,:) = 1/2*(ma_diff0(:,1,:)+ma_diff0(:,12,:));
        % ma_diff = permute(ma_diff, [3 1 2]); % bring plev to front
        % ma_diff = squeeze(interp1(grid.dim3.plev, ma_diff, par.pa_eval)); % evaluate ma_difference at plev_eval

        ga_diff0 = ga_diff_orig.(land);
        ga_diff = nan([size(ga_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_diff(:,2:13) = ga_diff0;
        ga_diff(:,1) = 1/2*(ga_diff0(:,1)+ga_diff0(:,12));
        ga_diff(:,14) = 1/2*(ga_diff0(:,1)+ga_diff0(:,12));

        ga_malr_diff0 = ga_malr_diff_orig.(land);
        ga_malr_diff = nan([size(ga_malr_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_malr_diff(:,2:13) = ga_malr_diff0;
        ga_malr_diff(:,1) = 1/2*(ga_malr_diff0(:,1)+ga_malr_diff0(:,12));
        ga_malr_diff(:,14) = 1/2*(ga_malr_diff0(:,1)+ga_malr_diff0(:,12));

        ga_bl_diff0 = ga_bl_diff_orig.(land);
        ga_bl_diff = nan([size(ga_bl_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_bl_diff(:,2:13) = ga_bl_diff0;
        ga_bl_diff(:,1) = 1/2*(ga_bl_diff0(:,1)+ga_bl_diff0(:,12));
        ga_bl_diff(:,14) = 1/2*(ga_bl_diff0(:,1)+ga_bl_diff0(:,12));

        ga_dalr_bl_diff0 = ga_dalr_bl_diff_orig.(land);
        ga_dalr_bl_diff = nan([size(ga_dalr_bl_diff0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        ga_dalr_bl_diff(:,2:13) = ga_dalr_bl_diff0;
        ga_dalr_bl_diff(:,1) = 1/2*(ga_dalr_bl_diff0(:,1)+ga_dalr_bl_diff0(:,12));
        ga_dalr_bl_diff(:,14) = 1/2*(ga_dalr_bl_diff0(:,1)+ga_dalr_bl_diff0(:,12));

        inv_str0 = squeeze(inv.(land)(:,:,2));
        inv_str = nan([size(inv_str0,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
        inv_str(:,2:13) = inv_str0;
        inv_str(:,1) = 1/2*(inv_str0(:,1)+inv_str0(:,12));
        inv_str(:,14) = 1/2*(inv_str0(:,1)+inv_str0(:,12));

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

                % MA contour lat x mon dependence of RCE and RAE
                % figure(); clf; hold all;
                % cmp = colCog(10);
                % colormap(cmp);
                % imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                % contour(mesh_lat, mesh_mon, abs(ma_diff), par.ta_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $|T - T_m|_{500\\,\\mathrm{hPa}} < %g$ K', par.model, crit_text, land_text, par.ta_thresh)); end;
                % caxis([-2 2]);
                % xlabel('Month'); ylabel('Latitude (deg)');
                % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ma_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                % close;

                % GA contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', par.model, crit_text, land_text, par.ga_thresh));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh)); end
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA MALR contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_malr_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', par.model, crit_text, land_text, par.ga_thresh));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$', upper(type), crit_text, land_text, par.ga_thresh)); end
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_malr_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA_BL contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, ga_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', upper(type), crit_text, land_text, par.ga_bl_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', par.model, crit_text, land_text, par.ga_bl_thresh));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', upper(type), crit_text, land_text, par.ga_bl_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_bl_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
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

                % INV contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, inv_str, par.inv_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $T_{\\sigma = 0.85} - T_s>%g$ K', upper(type), crit_text, land_text, par.inv_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $T_{\\sigma = 0.85} - T_s>%g$ K', par.model, crit_text, land_text, par.inv_thresh));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, %s, %s, $T_{\\sigma = 0.85} - T_s>%g$ K', upper(type), crit_text, land_text, par.inv_thresh)); end
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_inv_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA and GA_BL contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, ga_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_m - \\Gamma )/ \\Gamma_m > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_ga_bl_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA MALR and GA_BL DALR contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_malr_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_malr_ga_dalr_bl_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % GA and INV contour lat x mon dependence of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                contour(mesh_lat, mesh_mon, abs(ga_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, inv_str, par.inv_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $T_{\\sigma = 0.85} - T_s>%g$ K', upper(type), par.ga_thresh, par.inv_thresh));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $T_{\\sigma = 0.85} - T_s>%g$ K', par.model, par.ga_thresh, par.inv_thresh));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $|\\Gamma_m - \\Gamma |/ \\Gamma_m < %g \\%%$, $T_{\\sigma = 0.85} - T_s>%g$ K', upper(type), par.ga_thresh, par.inv_thresh)); end
                caxis([-2 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_ga_inv_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % R1 GA and INV contour lat x mon dependence of RCE and RAE
                r10 = flux_z.(land).r1.(fw);
                r1 = nan([size(r10,1) 14]); % prepare to add months 0.5 and 12.5 (cleaner plot when overlaying contour with imagesc)
                r1(:,2:13) = r10;
                r1(:,1) = 1/2*(r10(:,1)+r10(:,12));
                r1(:,14) = 1/2*(r10(:,1)+r10(:,12));
                var_text = '$R_1$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(flipud(cmp));
                contour(mesh_lat, mesh_mon, r1, par.ep*[1 1], '--', 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, r1, par.ga*[1 1], '--', 'color', par.blue, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, abs(ga_diff), par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
                contour(mesh_lat, mesh_mon, inv_str, par.inv_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text));
                elseif any(strcmp(type, {'echam'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text)); end
                caxis([-1 2]);
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_r1_mon_lat_ga_inv_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                close;

                % % ALBEDO contour lat x mon dependence of RCE and RAE
                % figure(); clf; hold all;
                % cmp = colCog(10);
                % colormap(cmp);
                % imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                % contour(mesh_lat, mesh_mon, albedo, par.albedo_thresh*[1 1], 'k', 'linewidth', 2);
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\alpha=%g$', par.model, crit_text, land_text, par.albedo_thresh)); end;
                % caxis([-2 2]);
                % xlabel('Month'); ylabel('Latitude (deg)');
                % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_albedo_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                % close;

                % % SN contour lat x mon dependence of RCE and RAE
                % figure(); clf; hold all;
                % cmp = colCog(10);
                % colormap(cmp);
                % imagesc([1 12], [lat(1) lat(end)], rcae_alt_z.(land).(fw).(crit));
                % contour(mesh_lat, mesh_mon, sn, par.sn_thresh*[1 1], 'k', 'linewidth', 2);
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), crit_text, land_text));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $d=%g$ m', par.model, crit_text, land_text, par.sn_thresh)); end;
                % caxis([-2 2]);
                % xlabel('Month'); ylabel('Latitude (deg)');
                % set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/0_rcae_alt_mon_lat_sn_overlay', par.plotdir, par.ep, par.ga, fw, crit, land), '-dpng', '-r300');
                % close;

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
function plot_temp(type, par)
    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ta.mat', type, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ma.mat', type, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ta.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ma.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
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

function [flux_zt, vh, vh_mon, lat, par] = load_flux(type, par)
    % load processed data/proc
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh_mon.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/flux_zt.mat', type, par.model, par.gcm.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/vh.mat', type, par.model, par.gcm.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/vh_mon.mat', type, par.model, par.gcm.clim, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/flux_zt.mat', type, par.echam.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh.mat', type, par.echam.clim, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/vh_mon.mat', type, par.echam.clim, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    elseif strcmp(type, 'echam_ml')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh_mon.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    end
end
function make_dirs(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        for plev_eval = [300:100:500]
            if ~exist(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land), 'dir')
                mkdir(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land));
            end
        end
        for si_eval = [0.8 0.85 0.9]
            if ~exist(sprintf('%s/inv_str/si_%g/%s', par.plotdir, si_eval, land), 'dir')
                mkdir(sprintf('%s/inv_str/si_%g/%s', par.plotdir, si_eval, land));
            end
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/energy-flux/%s/%s', par.plotdir, land, time));
            end
            if ~exist(sprintf('%s/transport/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/transport/%s/%s', par.plotdir, land, time));
            end
            if ~exist(sprintf('%s/ga_diff/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/ga_diff/%s/%s', par.plotdir, land, time));
            end
            if ~exist(sprintf('%s/ga_malr_diff/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/ga_malr_diff/%s/%s', par.plotdir, land, time));
            end
            if par.do_surf; v_vec = {'p', 'z', 'si', 'pi'};
            else v_vec = {'p', 'z', 'si'}; end
            for v = v_vec; vert = v{1};
                if ~exist(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert));
                end
            end
            if any(strcmp(type, {'era5', 'erai'}))
                for f = par.era.fw; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            elseif strcmp(type, 'gcm')
                for f = par.gcm.fw; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            elseif strcmp(type, 'echam')
                for f = par.echam.fw; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            elseif strcmp(type, 'echam_ml')
                for f = par.echam.fw; fw = f{1};
                    if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr1/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                    if ~exist(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                        mkdir(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time));
                    end
                end
            end
            if ~exist(sprintf('%s/flag/%s/%s', par.plotdir, land, time), 'dir')
                mkdir(sprintf('%s/flag/%s/%s', par.plotdir, land, time));
            end
        end
    end

    if ~exist(sprintf('%s/va', par.plotdir), 'dir')
        mkdir(sprintf('%s/va', par.plotdir));
    end
    if ~exist(sprintf('%s/trop', par.plotdir), 'dir')
        mkdir(sprintf('%s/trop', par.plotdir));
    end
    if ~exist(sprintf('%s/alb', par.plotdir), 'dir')
        mkdir(sprintf('%s/alb', par.plotdir));
    end
    if ~exist(sprintf('%s/sn', par.plotdir), 'dir')
        mkdir(sprintf('%s/sn', par.plotdir));
    end
    if ~exist(sprintf('%s/legends', par.plotdir), 'dir')
        mkdir(sprintf('%s/legends', par.plotdir));
    end
end
function make_dirs_ep(type, par)
    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw;
        elseif strcmp(type, 'echam_ml'); f_vec = par.echam.fw; end
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
