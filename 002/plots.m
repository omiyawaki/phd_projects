clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

gcm_info
echam_info

%% set parameters
% lat grid type
if 1
par.lat_interp = 'native'; % native: native model grid, don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.gcm.clim = 'piControl'; % choose either piControl or abrupt4xCO2
par.echam_clims = par.echam.sel; % par.echam.all_mld; % choose from 20170908 (snowball), 20170915_2 (modern), or rp000*** (various mixed layer depth and with/without sea ice)
par.ep_swp = 0.1; %[0.25 0.3 0.35]; % threshold value for determining RCE
par.ga_swp = 0.9; % threshold for determining RAE
par.si_eval = [0.8 0.85 0.9]; % sigma level for calculating inversion strength
par.ma_init = 0.95; % 'surf' for initializing with 2 m data, otherwise enter starting sigma level for moist adiabat
par.pa_eval = 500e2; % pressure level for calculating ma_diff
par.si_bl_swp = [0.85 0.9 0.95]; % sigma level to separate vertical average for close to moist adiabatic and stable surface stratifications
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
par.era.fw = {'mse', 'dse'};
par.merra2.fw = {'mse', 'dse'};
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
end

% set default figure parameters
if 1
    par.monlabel = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};
    par.monlabelsh = {'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J'};
    par.ppos = [0 0 10/3 7/3];
    par.ppos_larger = [0 0 16/3 13/3];
    par.ppos_vert = [0 0 7/3 10/3];
    par.ppos_sq = [0 0 10/3 10/3];
    par.ppos_wide = [0 0 13/3 7/3];
    par.ppos_verywide = [0 0 16/3 7/3];
    par.ppos_superwide = [0 0 21/3 7/3];
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

%% call functions
% plot_rad_lat(par)
% plot_rad_lon_lat(par)
% plot_tediv_lat(par)

type = 'era5';
% choose_plots(type, par);
for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
    type='echam';
    % disp(par.echam.clim)
    % choose_plots(type, par);
end
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    disp(par.model)
    choose_plots(type, par);
end

% % sweep through various boundary layer heights
% for i = 1:length(par.si_bl_swp); par.si_bl = par.si_bl_swp(i);
%     type = 'merra2';
%     % choose_plots_si_bl(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam';
%         % disp(par.echam.clim)
%         % choose_plots_si_bl(type, par);
%     end
%     for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
%         type = 'gcm';
%         disp(par.model)
%         choose_plots_si_bl(type, par)
%     end
% end

% % sweep through various threshold values
% for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i); par.ga = par.ga_swp(i);
%     type = 'era5';
%     % choose_plots_ep(type, par)
%     for k=1:length(par.echam_clims); par.echam.clim=par.echam_clims{k};
%         type='echam';
%         disp(par.echam.clim)
%         choose_plots_ep(type, par);
%     end
%     for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
%         type = 'gcm';
%         % disp(par.model)
%         % choose_plots_ep(type, par)
%     end
% end

function choose_plots(type, par)
    % plot_temp_zon(type, par) % plot temperature profiles at specific latitudes
    % plot_temp_zon_select(type, par) % plot temperature profiles at specific latitudes
    % plot_thetaeq_zon_select(type, par) % plot eq pot temp profiles at specific latitudes
    % plot_dtdz_zon_select(type, par) % plot temperature profiles at specific latitudes
    % plot_ma_diff(type, par) % plot difference of temperature profile from moist adiabat
    % plot_ga_diff(type, par) % plot difference of temperature profile from moist adiabat
    % plot_ga_malr_si_diff(type, par) % plot difference of temperature profile from moist adiabat in sigma
    % plot_dmse_polar_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_dmse_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    plot_dmse_toasfc_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space
    % plot_siced(type, par) % sea ice depth
    % plot_sftlf(type, par) % land fraction

    % plot_dmse_toasfc_midlatitude_line(type, par) % plot decomposition of R1 in mon x lat and lon x lat space

end % select which functions to run at a time
function choose_plots_si_bl(type, par)
    plot_ga_malr_diff_mon_lat(type, par) % plot ga diff
    % plot_ga_malr_diff_lon_lat(type, par) % plot ga diff
end
function choose_plots_ep(type, par)
    % plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
    % plot_r1z_lat(type, par); % compare r1 line plot with ERA5
    % plot_flux(type, par) % plot various energy fluxes in mon x lat and lon x lat space
    % plot_temp(type, par) % plot temperature profiles
    % plot_temp_ann(type, par) % plot temperature profiles
    plot_dr1_midlatitude_line(type, par) % plot decomposition of R1 in mon at specific latitudes
    % plot_dr1_polar_line(type, par) % plot decomposition of R1 in mon at specific latitudes
end % select which ep-functions to run at a time

%% define functions
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
    elseif strcmp(type, 'echam_pl')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
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
                    clear h h_label;
                end % hemi
            end % vert

            figure(); clf; hold all; box on;
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

            figure(); clf; hold all; box on;
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
function plot_mlev(type, par)
    make_dirs(type, par)

    if strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    else
        error('This code only works for data output in the model vertical grid.')
    end
    load(sprintf('%s/grid.mat', prefix));

    % compute sigma from a and b
    ps_vert = repmat(ps_orig, [1 1 1 length(grid.dim3.a)]); % dims (lon x lat x time x plev)
    ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
    b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
    pa = a + b.*ps_vert;

    % time and zonal mean
    pa = squeeze(nanmean(nanmean(pa,1),4));

    % standard pl grid
    pl = 1e2*[1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1 0.5 0.2 0.1];

    % pressure x lat
    figure(); clf; hold all; box on;
    for lev = 1:length(grid.dim3.a)
        h1=plot(grid.dim3.lat, 1e-2*pa(:,lev), '-k', 'linewidth', 0.5);
    end
    for lev = 1:length(pl); plev = pl(lev);
        h2=line([-90 90], 1e-2*plev*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 0.5);
    end
    legend([h1 h2], 'Model levels', 'Standard pressure levels', 'location', 'southoutside')
    xlabel('latitude (deg)'); ylabel('p (hPa)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on')
    print(sprintf('%s/mlev/mlev', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_temp_zon_select(type, par)
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
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 6]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July'; end;

            tasi_mon.(land) = squeeze(tasi.(land)(:,month,:));
            masi_mon.(land) = squeeze(masi.(land)(:,month,:));

            % remove moist adiabat data below initialization level
            if ~strcmp(par.ma_init, 'surf')
                masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
            end

            tasi_sp(:,m) = interp1(lat, tasi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tasi_np(:,m) = interp1(lat, tasi_mon.(land), lat_pole); % sounding at lat_pole N
            tasi_smid(:,m) = interp1(lat, tasi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tasi_nmid(:,m) = interp1(lat, tasi_mon.(land), lat_mid); % sounding at lat_mid N
            tasi_eq(:,m) = interp1(lat, tasi_mon.(land), 0); % sounding at equator

            masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
            masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
            masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
            masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
            masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator

            % ALL
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, '--', 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, '--', 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_np h_sp h_nmid h_smid], sprintf('%g N', lat_pole), sprintf('%g S', lat_pole), sprintf('%g N', lat_mid), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/all', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
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
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
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
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nh_only_vert', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sh_only_vert', par.plotdir, land, month), '-dpng', '-r300');
            close;

        end

        % ALL NH
        figure(); clf; hold all; box on;
        h_nmid_jan = plot(tasi_nmid(:,1), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_nmid_ma_jan = plot(masi_nmid(:,1), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_nmid_jun = plot(tasi_nmid(:,6), grid.dim3.si, 'color', par.orange);
        h_nmid_ma_jun = plot(masi_nmid(:,6), grid.dim3.si, ':', 'color', par.orange);
        h_np_jan = plot(tasi_np(:,1), grid.dim3.si, 'color', par.blue);
        h_np_jun = plot(tasi_np(:,6), grid.dim3.si, 'color', 0.25*[1 1 1]);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, NH', upper(type)));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, NH', par.model));
        elseif strcmp(type, 'echam')
            title(sprintf('%s, NH', upper(type)));
        end
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [200 290], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/temp_zon_sel/%s/nh_all', par.plotdir, land), '-dpng', '-r300');
        close;

        % ALL SH
        figure(); clf; hold all; box on;
        h_smid_jan = plot(tasi_smid(:,1), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_smid_ma_jan = plot(masi_smid(:,1), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_smid_jun = plot(tasi_smid(:,6), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_smid_ma_jun = plot(masi_smid(:,6), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_sp_jan = plot(tasi_sp(:,1), grid.dim3.si, 'color', par.blue);
        h_sp_jun = plot(tasi_sp(:,6), grid.dim3.si, 'color', par.blue);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        if strcmp(type, 'era5') | strcmp(type, 'erai')
            title(sprintf('%s, SH', upper(type)));
        elseif strcmp(type, 'gcm')
            title(sprintf('%s, SH', par.model));
        elseif strcmp(type, 'echam')
            title(sprintf('%s, SH', upper(type)));
        end
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [200 290], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/temp_zon_sel/%s/sh_all', par.plotdir, land), '-dpng', '-r300');
        close;

    end
end
function plot_thetaeq_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/thetaeq_mon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/thetaeq_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/thetaeq_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/thetaeq_mon_lat.mat', type, par.lat_interp));
        plev = 1:47;
    end

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 7]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==7; mon_str = 'July'; end;

            thetaeqsi_mon.(land) = squeeze(thetaeqsi.(land)(:,month,:));

            thetaeqsi_sp = interp1(lat, thetaeqsi_mon.(land), -lat_pole); % sounding at -lat_pole S
            thetaeqsi_np = interp1(lat, thetaeqsi_mon.(land), lat_pole); % sounding at lat_pole N
            thetaeqsi_smid = interp1(lat, thetaeqsi_mon.(land), -lat_mid); % sounding at -lat_mid S
            thetaeqsi_nmid = interp1(lat, thetaeqsi_mon.(land), lat_mid); % sounding at lat_mid N
            thetaeqsi_eq = interp1(lat, thetaeqsi_mon.(land), 0); % sounding at equator

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(thetaeqsi_np, grid.dim3.si, 'color', par.blue);
                h_nmid = plot(thetaeqsi_nmid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            elseif m == 7
                h_np = plot(thetaeqsi_np, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(thetaeqsi_nmid, grid.dim3.si, 'color', par.orange);
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
            set(gca, 'fontsize', par.fs, 'xlim', [250 330], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/nh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(thetaeqsi_sp, grid.dim3.si, 'color', par.blue);
            h_smid = plot(thetaeqsi_smid, grid.dim3.si, 'color', 0.25*[1 1 1]);
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
            set(gca, 'fontsize', par.fs, 'xlim', [250 330], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/sh_only', par.plotdir, land, month), '-dpng', '-r300');
            close;

        end
    end
end
function plot_dtdz_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
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
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
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
            if strcmp(type, 'era5') | strcmp(type, 'erai')
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
function plot_ga_malr_si_diff(type, par)
% sigma coordinates
    make_dirs(type, par)

    % load data
    % [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/dtdzsi.mat', prefix)); dtdzzsi = dtdzsi; clear dtdzsi; % read model lapse rate
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
    [C, h] = contourf(mesh_si, mesh_lat, diff_zt', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_jul', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, July', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, July', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_jul', par.plotdir), '-dpng', '-r300');
    close;

    % lat x plev of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, diff_jan', [-200 -100 -50 -20 -10 0 10 20 50 100 200]);
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-100 100]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, January', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, January', par.model));
    end
    cb = colorbar('limits', [-100 100], 'ytick', [-100:20:100], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$(\Gamma_m - \Gamma)/\Gamma_m$ (\%)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_malr_diff_lat_sigma_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of lapse rate percentage diff
    figure(); clf; hold all;
    cmp = colCog(40);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, diff_zt', [-200 -100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
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
    [C, h] = contour(mesh_si, mesh_lat, diff_jan', [-200 -100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
    clabel(C, h, [-100 -50 -20 -10 0 10 20 50 100], 'fontsize', 6, 'interpreter', 'latex');
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
    [C, h] = contour(mesh_si, mesh_lat, diff_jul', [-100 -50 -20 -10 0 10 20 50 100 200], 'color', 'k');
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
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_zt', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JANUARY
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jan', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
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
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_lat_sigma_jan', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of model lapse rate JULY
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtdz_jul', [-20 -10:2:10 20], 'linecolor', 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
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
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
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
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contourf(mesh_si, mesh_lat, dtmdz_zt', -10:2:10, 'w');
    clabel(C, h, [-10:2:10], 'fontsize', 6, 'interpreter', 'latex', 'color', 'w');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    cb = colorbar('limits', [0 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gca, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yscale', 'linear', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma', par.plotdir), '-dpng', '-r300');
    close;

    % lat x sigma of moist adiabatic lapse rate 10 S to 75 N
    figure(); clf; hold all; box on;
    cmp = colCog(10);
    colormap(cmp);
    [C, h] = contour(mesh_si, mesh_lat, dtmdz_zt', [5 6 7 8 9], 'color', 'k');
    clabel(C, h, [5 6 7 8 9], 'fontsize', 6, 'interpreter', 'latex');
    caxis([-10 10]);
    xlabel('Latitude (deg)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s, Annual', upper(type)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s, Annual', par.model));
    end
    % cb = colorbar('limits', [-10 10], 'ytick', [-10:2:10], 'location', 'eastoutside');
    % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    % ylabel(cb, '$\Gamma_m$ (K/km)');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'xlim', [-10 75], 'xtick', [-10:10:70 75], 'ylim', [0.2 1], 'ytick', [0.2:0.1:1], 'ydir', 'reverse', 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/ga_malr_diff/ga_m_lat_sigma_sc', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_dmse_polar_line(type, par)
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
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [-85 -80 -70 70 80 85];
    lat_bound_list = [-80 80];

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'L+O';
        elseif strcmp(land, 'l'); land_text = 'L';
        elseif strcmp(land, 'o'); land_text = 'O';
        end
        if strcmp(type, 'echam'); land_text = par.echam.(par.echam.clim); end;
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
        elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby=0; monlabel=par.monlabel;
                else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby=6; monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g', par.plotdir, fw, land, lat_bound);
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
                if any(strcmp(type,{'era5', 'erai'}))
                    lh = -flux_z.(land).slhf;
                    sh = -flux_z.(land).sshf;
                elseif any(strcmp(type,{'merra2'}))
                    lh = flux_z.(land).EFLUX;
                    sh = flux_z.(land).HFLUX;
                elseif strcmp(type, 'gcm')
                    lh = flux_z.(land).hfls;
                    sh = flux_z.(land).hfss;
                elseif strcmp(type, 'echam')
                    lh = -flux_z.(land).ahfl;
                    sh = -flux_z.(land).ahfs;
                end
                lh_lat = interp1(grid.dim3.lat, lh, lat);
                lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                sh_lat = interp1(grid.dim3.lat, sh, lat);
                sh_lat = nansum(sh_lat.*clat_mon)/nansum(clat);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res lhf shf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                % NOLEG ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                % legend([ra res stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'ylim', [-170 30], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse_noleg', folder), '-dpng', '-r300');
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
                dra=plot([1:12],  circshift(dra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(dres_lat,shiftby,2), 'color', par.maroon);
                dstf=plot([1:12], circshift(dstf_lat,shiftby,2), 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                ra=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
                sw=plot([1:12], circshift(sw_lat,shiftby,2), 'color', par.blue);
                lw=plot([1:12], circshift(lw_lat,shiftby,2), 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra sw lw], '$R_a$', '$\mathrm{Net SW}$', '$\mathrm{Net LW}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                dra=plot([1:12], circshift(dra_lat,shiftby,2), 'color', 0.5*[1 1 1]);
                dsw=plot([1:12], circshift(dsw_lat,shiftby,2), 'color', par.blue);
                dlw=plot([1:12], circshift(dlw_lat,shiftby,2), 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dsw dlw], '$\Delta R_a$', '$\Delta(\mathrm{Net SW})$', '$\Delta (\mathrm{Net LW})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                    stf=plot([1:12], circshift(stf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(hfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(hfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole));
                    xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\mathrm{LH+SH}$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                    stf=plot([1:12], circshift(dstf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(dhfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(dhfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_bound, lat_pole));
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\Delta(\mathrm{LH+SH})$', '$\Delta\mathrm{LH}$', '$\Delta\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dstf', folder), '-dpng', '-r300');
                    close;

                end

            end

        end % for mse dse
    end % for land

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    % elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
    %         dlat = 0.25; % step size for standard lat grid
    %         if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; shiftby=0; monlabel=par.monlabel;
    %         else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; shiftby=6; monlabel=par.monlabelsh; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dmse/%s/0_poleward_of_lat_%g', par.plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         ra = flux_z.lo.ra.(fw);
    %         ra_lat = interp1(grid.dim3.lat, ra, lat);
    %         ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.ra.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.ra.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RA lat x mon dependence of RCE and RAE
    %         var_text = '$R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$R_a$', '$R_{a,\mathrm{\,L+O}}$', '$R_{a,\mathrm{\,L}}$', '$R_{a,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_ra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         ra_ann = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
    %         ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
    %         dra = flux_z.lo.ra.(fw) - ra_ann;
    %         dra_lat = interp1(grid.dim3.lat, dra, lat);
    %         dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
    %         ra_ann_l = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.ra.(fw) - ra_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         ra_ann_o = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.ra.(fw) - ra_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RA lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(R_a)$', '$\Delta(R_{a,\mathrm{\,L+O}})$', '$\Delta(R_{a,\mathrm{\,L}})$', '$\Delta(R_{a,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res = flux_z.lo.res.(fw);
    %         res_lat = interp1(grid.dim3.lat, res, lat);
    %         res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.res.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.res.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RES lat x mon dependence of RCE and RAE
    %         var_text = '$\nabla\cdot F_m$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(res_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\nabla\cdot F_m$', '$\nabla\cdot F_{m,\mathrm{\,L+O}}$', '$\nabla\cdot F_{m,\mathrm{\,L}}$', '$\nabla\cdot F_{m,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\nabla\\cdot F_m$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_res_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res_ann = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
    %         res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
    %         dres = flux_z.lo.res.(fw) - res_ann;
    %         dres_lat = interp1(grid.dim3.lat, dres, lat);
    %         dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
    %         res_ann_l = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.res.(fw) - res_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         res_ann_o = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.res.(fw) - res_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RES lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta(\nabla\cdot F_m)$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dres_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\nabla\cdot F_m)$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L+O}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta(\\nabla\\cdot F_m)$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dres_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf = flux_z.lo.stf.(fw);
    %         stf_lat = interp1(grid.dim3.lat, stf, lat);
    %         stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.stf.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.stf.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % STF lat x mon dependence of RCE and STFE
    %         var_text = '$\mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(stf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\mathrm{LH+SH}$', '$(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$(\mathrm{LH+SH})_{\mathrm{L}}$', '$(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_stf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf_ann = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
    %         stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
    %         dstf = flux_z.lo.stf.(fw) - stf_ann;
    %         dstf_lat = interp1(grid.dim3.lat, dstf, lat);
    %         dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
    %         stf_ann_l = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.stf.(fw) - stf_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         stf_ann_o = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.stf.(fw) - stf_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA STF lat x mon dependence of RCE and STFE
    %         var_text = '$\Delta \mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dstf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\mathrm{LH+SH})$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\delta \\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dstf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %     end
    % end

end % for function
function plot_dmse_midlatitude_line(type, par)
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
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [5 -5];

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'L+O';
        elseif strcmp(land, 'l'); land_text = 'L';
        elseif strcmp(land, 'o'); land_text = 'O';
        end
        if strcmp(type, 'echam'); land_text = par.echam.(par.echam.clim); end;
        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
        elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
        for f = f_vec; fw = f{1};
            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
                else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dmse/%s/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, land, lat_center-lat_bound);
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
                if any(strcmp(type,{'era5', 'erai'}))
                    lh = -flux_z.(land).slhf;
                    sh = -flux_z.(land).sshf;
                elseif any(strcmp(type,{'merra2'}))
                    lh = flux_z.(land).EFLUX;
                    sh = flux_z.(land).HFLUX;
                elseif strcmp(type, 'gcm')
                    lh = flux_z.(land).hfls;
                    sh = flux_z.(land).hfss;
                elseif strcmp(type, 'echam')
                    lh = -flux_z.(land).ahfl;
                    sh = -flux_z.(land).ahfs;
                end
                lh_lat = interp1(grid.dim3.lat, lh, lat);
                lh_lat = nansum(lh_lat.*clat_mon)/nansum(clat);
                sh_lat = interp1(grid.dim3.lat, sh, lat);
                sh_lat = nansum(sh_lat.*clat_mon)/nansum(clat);

                % ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra res lhf, shf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside', 'numcolumns', 2);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'ylim', [-150 100], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse', folder), '-dpng', '-r300');
                close;

                % NOLEG ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                ra=plot([1:12],  circshift(ra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                res=plot([1:12], circshift(res_lat,shiftby,2), 'color', par.maroon);
                lhf=plot([1:12], circshift(lh_lat,shiftby,2), 'color', par.blue);
                shf=plot([1:12], circshift(sh_lat,shiftby,2), 'color', par.orange);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                % legend([ra res stf], '$R_a$', '$\nabla\cdot F_m$', '$\mathrm{LH+SH}$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'ylim', [-150 100], 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_mse_noleg', folder), '-dpng', '-r300');
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
                dra=plot([1:12],  circshift(dra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(dres_lat,shiftby,2), 'color', par.maroon);
                dstf=plot([1:12], circshift(dstf_lat,shiftby,2), 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse', folder), '-dpng', '-r300');
                close;

                % NOLEG DELTA ALL lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                dra=plot([1:12],  circshift(dra_lat ,shiftby,2), 'color', 0.5*[1 1 1]);
                dres=plot([1:12], circshift(dres_lat,shiftby,2), 'color', par.maroon);
                dstf=plot([1:12], circshift(dstf_lat,shiftby,2), 'color', par.blue);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                % legend([dra dres dstf], '$\Delta R_a$', '$\Delta(\nabla\cdot F_m)$', '$\Delta (\mathrm{LH+SH})$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dmse_noleg', folder), '-dpng', '-r300');
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
                ra=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
                sw=plot([1:12], circshift(sw_lat,shiftby,2), 'color', par.blue);
                lw=plot([1:12], circshift(lw_lat,shiftby,2), 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([ra sw lw], '$R_a$', '$\mathrm{Net SW}$', '$\mathrm{Net LW}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                dra=plot([1:12], circshift(dra_lat,shiftby,2), 'color', 0.5*[1 1 1]);
                dsw=plot([1:12], circshift(dsw_lat,shiftby,2), 'color', par.blue);
                dlw=plot([1:12], circshift(dlw_lat,shiftby,2), 'color', par.maroon);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, lat_center-lat_bound, lat_center+lat_bound));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound)); end;
                % xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([dra dsw dlw], '$\Delta R_a$', '$\Delta(\mathrm{Net SW})$', '$\Delta (\mathrm{Net LW})$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                    stf=plot([1:12], circshift(stf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(hfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(hfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound));
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\mathrm{LH+SH}$', '$\mathrm{LH}$', '$\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
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
                    stf=plot([1:12], circshift(dstf_lat ,shiftby,2), 'k');
                    lh=plot([1:12],  circshift(dhfls_lat,shiftby,2), 'color', par.blue);
                    sh=plot([1:12],  circshift(dhfss_lat,shiftby,2), 'color', par.orange);
                    title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, lat_center-lat_bound, lat_center+lat_bound));
                    % xlabel('Month');
                    ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                    legend([stf lh sh], '$\Delta(\mathrm{LH+SH})$', '$\Delta\mathrm{LH}$', '$\Delta\mathrm{SH}$', 'location', 'eastoutside');
                    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
                    print(sprintf('%s/0_mon_dstf', folder), '-dpng', '-r300');
                    close;

                end

            end

        end % for mse dse
    end % for land

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    % elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
    %         dlat = 0.25; % step size for standard lat grid
    %         if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
    %         else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dmse/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         ra = flux_z.lo.ra.(fw);
    %         ra_lat = interp1(grid.dim3.lat, ra, lat);
    %         ra_lat = nansum(ra_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.ra.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.ra.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RA lat x mon dependence of RCE and RAE
    %         var_text = '$R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(ra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$R_a$', '$R_{a,\mathrm{\,L+O}}$', '$R_{a,\mathrm{\,L}}$', '$R_{a,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_ra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         ra_ann = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         ra_ann_lat = interp1(grid.dim3.lat, ra_ann, lat);
    %         ra_ann_lat = nansum(ra_ann_lat.*clat_mon)/nansum(clat);
    %         dra = flux_z.lo.ra.(fw) - ra_ann;
    %         dra_lat = interp1(grid.dim3.lat, dra, lat);
    %         dra_lat = nansum(dra_lat.*clat_mon)/nansum(clat);
    %         ra_ann_l = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.ra.(fw) - ra_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         ra_ann_o = repmat(nanmean(flux_z.lo.ra.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.ra.(fw) - ra_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RA lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_a$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dra_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$\Delta(R_a)$', '$\Delta(R_{a,\mathrm{\,L+O}})$', '$\Delta(R_{a,\mathrm{\,L}})$', '$\Delta(R_{a,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_a$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dra_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res = flux_z.lo.res.(fw);
    %         res_lat = interp1(grid.dim3.lat, res, lat);
    %         res_lat = nansum(res_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.res.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.res.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % RES lat x mon dependence of RCE and RAE
    %         var_text = '$\nabla\cdot F_m$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(res_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$\nabla\cdot F_m$', '$\nabla\cdot F_{m,\mathrm{\,L+O}}$', '$\nabla\cdot F_{m,\mathrm{\,L}}$', '$\nabla\cdot F_{m,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\nabla\\cdot F_m$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_res_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         res_ann = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         res_ann_lat = interp1(grid.dim3.lat, res_ann, lat);
    %         res_ann_lat = nansum(res_ann_lat.*clat_mon)/nansum(clat);
    %         dres = flux_z.lo.res.(fw) - res_ann;
    %         dres_lat = interp1(grid.dim3.lat, dres, lat);
    %         dres_lat = nansum(dres_lat.*clat_mon)/nansum(clat);
    %         res_ann_l = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.res.(fw) - res_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         res_ann_o = repmat(nanmean(flux_z.lo.res.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.res.(fw) - res_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA RES lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta(\nabla\cdot F_m)$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dres_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\nabla\cdot F_m)$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L+O}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,L}})$', '$\Delta(\nabla\cdot F_{m,\mathrm{\,O}})$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta(\\nabla\\cdot F_m)$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dres_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf = flux_z.lo.stf.(fw);
    %         stf_lat = interp1(grid.dim3.lat, stf, lat);
    %         stf_lat = nansum(stf_lat.*clat_mon)/nansum(clat);
    %         comp1 = sftlf*1e-2.*flux_z.l.stf.(fw);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         comp2 = (1-sftlf*1e-2).*flux_z.o.stf.(fw);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % STF lat x mon dependence of RCE and STFE
    %         var_text = '$\mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(stf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':k', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$\mathrm{LH+SH}$', '$(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$(\mathrm{LH+SH})_{\mathrm{L}}$', '$(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_stf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %         stf_ann = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         stf_ann_lat = interp1(grid.dim3.lat, stf_ann, lat);
    %         stf_ann_lat = nansum(stf_ann_lat.*clat_mon)/nansum(clat);
    %         dstf = flux_z.lo.stf.(fw) - stf_ann;
    %         dstf_lat = interp1(grid.dim3.lat, dstf, lat);
    %         dstf_lat = nansum(dstf_lat.*clat_mon)/nansum(clat);
    %         stf_ann_l = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.stf.(fw) - stf_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);
    %         stf_ann_o = repmat(nanmean(flux_z.lo.stf.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.stf.(fw) - stf_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);
    %         % DELTA STF lat x mon dependence of RCE and STFE
    %         var_text = '$\Delta \mathrm{LH+SH}$';
    %         figure(); clf; hold all; box on;
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dstf_lat,shiftby,2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat,shiftby,2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat,shiftby,2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat,shiftby,2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_center-lat_bound, lat_center+lat_bound));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, var_text, lat_center-lat_bound, lat_center+lat_bound)); end;
    %         legend([tot c12 c1 c2], '$\Delta(\mathrm{LH+SH})$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L+O}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{L}}$', '$\Delta(\mathrm{LH+SH})_{\mathrm{O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\delta \\mathrm{LH+SH}$ (Wm$^{-2}$)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dstf_lo_decomp', folder), '-dpng', '-r300');
    %         close;

    %     end
    % end

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

                lwsfc = flux_z.(land).lwsfc;
                lwsfc_lat = interp1(grid.dim3.lat, lwsfc, lat);
                lwsfc_lat = nansum(lwsfc_lat.*clat_mon)/nansum(clat);

                % SWSFC and LWSFC lat x mon dependence of RCE and RAE
                figure(); clf; hold all; box on;
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                swsfc=plot([1:12], -swsfc_lat, 'color', 'b');
                lwsfc=plot([1:12], lwsfc_lat, 'color', 'r');
                max(lwsfc_lat)
                min(lwsfc_lat)
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                xlabel('Month');
                ylabel(sprintf('Energy flux (Wm$^{-2}$)'));
                legend([swsfc lwsfc], '$F_\mathrm{SW,\,SFC}$', '$F_\mathrm{LW,\,SFC}$', 'location', 'eastoutside');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_swsfc_lwsfc', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

end % for function
function plot_siced(type, par)
% plot sicedow depth
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
    load(sprintf('%s/siced.mat', prefix)); % read clear sky albedo data

    [mesh_lat, mesh_mon] = meshgrid([1:12], lat);

    siced = squeeze(nanmean(siced,1)); % zonal mean

    % mon x lat of sicedow depth
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(flip(cmp));
    [C,h] = contourf(mesh_lat, mesh_mon, siced, [0:0.5:10 30 40 100], 'linecolor', 'none');
    % clabel(C, h, 0:0.5:10, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-5 5]);
    ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 5], 'ytick', [0:0.5:5], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Sea ice depth (m)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/siced/siced_mon_lat', par.plotdir), '-dpng', '-r300');
    close;

end
function plot_sftlf(type, par)
% plot sicedow depth
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
    load(sprintf('%s/sftlf.mat', prefix)); % read clear sky albedo data

    sftlf = squeeze(nanmean(sftlf,1)); % zonal mean

    % mon x lat of sicedow depth
    figure(); clf; hold all; box on;
    plot(grid.dim2.lat, sftlf, '-k');
    xlabel('Latitude (deg)'); ylabel('Land fraction (\%)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
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

% SI BL PLOTS
function plot_ga_malr_diff_mon_lat(type, par)
% plot inversion strength
make_dirs_si_bl(type, par)

% load data
% [~, ~, ~, lat, par] = load_flux(type, par);
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
elseif any(strcmp(type, {'echam_ml', 'echam_pl'}))
    par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
    prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
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
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
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
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
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
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
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
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
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
function plot_ga_malr_diff_lon_lat(type, par)
% plot inversion strength
    make_dirs_si_bl(type, par)

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
            if strcmp(type, 'era5') | strcmp(type, 'erai')
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
            if strcmp(type, 'era5') | strcmp(type, 'erai')
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

% EP PLOTS
function plot_energy_lat(type, par)
% latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, par] = load_flux(type, par); % load data
    make_dirs(type, par)

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        if any(strcmp(type, {'erai', 'era5'})); f_vec = par.era.fw;
        elseif strcmp(type, 'merra2'); f_vec = par.merra2.fw;
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
                    elseif strcmp(type, 'merra2')
                        if contains(fw, 'mse'); plot(lat, flux.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
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
                    elseif strcmp(type, 'merra2'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
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
                    elseif strcmp(type, 'merra2')
                        if strcmp(fw, 'mse'); plot(lat, flux.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).EFLUX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).PRECTOT,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).HFLUX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
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
                    elseif strcmp(type, 'merra2'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
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
                ylim_lo = min(flux.(land).(time).r1.(fw)); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max(flux.(land).(time).r1.(fw)); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
                rcemax = par.ep;
                vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                raemin = par.ga;
                vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                % line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                % if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                % line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                else; plot(lat,flux.(land).(time).r1.(fw), '-k'); end
                if any(strcmp(type, {'era5', 'erai'}));
                    title(sprintf('%s, %s', upper(type), upper(time)));
                elseif strcmp(type, 'merra2'); title(sprintf('%s, %s', upper(type), upper(time)));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s', par.model, upper(time)));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s', upper(type), upper(time)));
                end
                xlabel('latitude (deg)');
                if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                else ylabel('$R_1$ (unitless)'); end
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                print(sprintf('%s/energy-flux/%s/%s/%s-r1', par.plotdir, land, time, fw), '-dpng', '-r300');
                close;

                % R1Z
                figure(); clf; hold all; box on;
                r1z=flux.(land).(time).res.(fw)./flux.(land).(time).ra.(fw);
                ylim_lo = min(r1z); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max(r1z); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
                rcemax = par.ep;
                vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                raemin = par.ga;
                vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                % line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                % line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                % if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                % line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                if strcmp(fw, 'dse'); plot(lat,r1z, '--k');
                else; plot(lat,r1z, '-k'); end
                if any(strcmp(type, {'era5', 'erai'}));
                    title(sprintf('%s, %s', upper(type), upper(time)));
                elseif strcmp(type, 'merra2'); title(sprintf('%s, %s', upper(type), upper(time)));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s', par.model, upper(time)));
                elseif strcmp(type, 'echam'); title(sprintf('%s, %s', upper(type), upper(time)));
                end
                xlabel('latitude (deg)');
                if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                else ylabel('$R_1$ (unitless)'); end
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                print(sprintf('%s/energy-flux/%s/%s/%s-r1z', par.plotdir, land, time, fw), '-dpng', '-r300');
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
                else
                    if any(strcmp(type, {'erai', 'era5', 'merra2'}))
                        title(sprintf('%s, Northward %s Transport, %s', upper(type), upper(fw), upper(time)));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, Northward %s Transport, %s', par.model, upper(fw), upper(time)));
                    end
                end;
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
function plot_r1z_lat(type, par)
% Compare GCM R1 with ERA5

    grid_era5 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', 'era5'));
    grid_gcm = load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    flux_era5 = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', 'era5', par.lat_interp));
    flux_gcm = load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/flux_zt.mat', type, par.model, par.gcm.clim, par.lat_interp));
    par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);

    % R1Z
    figure(); clf; hold all; box on;
    r1z_era5=flux_era5.flux_zt.lo.ann.res.mse./flux_era5.flux_zt.lo.ann.ra.mse;
    r1z_gcm=flux_gcm.flux_zt.lo.ann.res.mse./flux_gcm.flux_zt.lo.ann.ra.mse;
    ylim_lo = min([r1z_era5, r1z_gcm]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
    ylim_up = max([r1z_era5, r1z_gcm]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
    raemin = par.ga;
    vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    rcemax = par.ep;
    vertices = [-90 rcemax; 90 rcemax; 90 raemin; -90 raemin];
    patch(vertices(:,1), vertices(:,2), 0.75*[1 1 1], 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    % plot(grid_gcm.grid.dim2.lat,r1z_gcm, '-k');
    plot(grid_era5.grid.dim2.lat,r1z_era5, '-k');
    % legend('RAE', 'RCAE', 'RCE', sprintf('%s', par.model), 'ERA5', 'location', 'eastoutside')
    xlabel('latitude (deg)');
    ylabel('$R_1$ (unitless)')
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/energy-flux/%s/%s/%s-r1z-comp-era5', par.plotdir, 'lo', 'ann', 'mse'), '-dpng', '-r300');
    close;

end
function plot_flux(type, par)
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
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    % load(sprintf('%s/%s/ga_malr_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_malr_diff_orig = ga_malr_diff; clear ga_diff;
    % load(sprintf('%s/%s/ga_dalr_bl_diff_si_mon_lat.mat', prefix_proc, par.lat_interp)); ga_dalr_bl_diff_orig = ga_dalr_bl_diff; clear ga_bl_diff;
    if strcmp(type, 'gcm')
        load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    else
        landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
        par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;
    end

    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        % ga_malr_diff = ga_malr_diff_orig.(land);
        % ga_dalr_bl_diff = ga_dalr_bl_diff_orig.(land);

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
        if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
        elseif strcmp(type, 'merra2'); f_vec = par.merra2.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw; end
        for f = f_vec; fw = f{1};
            % lat x mon dependence of RCE and RAE
            if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
            else var_text = '$\nabla \cdot F_m$'; end
            figure(); clf; hold all; box on;
            cmp = colCog(12);
            colormap(cmp);
            contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw), [-300 -150:25:150 300], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw), [0 0], 'color', 0.75*[1 1 1]);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'merra2'); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text));
            elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, land_text)); end;
            caxis([-150 150]);
            cb = colorbar('limits', [-150 150], 'ytick', [-150:25:150], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (Wm$^{-2}$)', var_text));
            xlabel('Month'); ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_div_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % R1 lat x mon dependence of RCE and RAE
            var_text = '$R_1$';
            figure(); clf; hold all; box on;
            cmp = colCog(10);
            colormap(flipud(cmp));
            contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [0 0], 'color', 0.75*[1 1 1]);
            contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
            contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
            [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
            clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'merra2'); title(sprintf('%s, %s, %s', upper(type), var_text, land_text));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s', par.model, var_text, land_text));
            elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s', upper(type), var_text, land_text)); end;
            caxis([-1 1]);
            cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (unitless)', var_text));
            ylabel('Latitude (deg)');
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_r1_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % R1z lat x mon dependence of RCE and RAE
            var_text = '$R_1$';
            figure(); clf; hold all; box on;
            cmp = colCog(10);
            colormap(flipud(cmp));
            contourf(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [-16 -8 -4 -2:0.2:2 4 8 16], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [0 0], 'color', 0.75*[1 1 1]);
            contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
            contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
            [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
            clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s', upper(type)));
            elseif strcmp(type, 'merra2'); title(sprintf('%s', upper(type)));
            elseif strcmp(type, 'gcm'); title(sprintf('%s', par.model));
            elseif strcmp(type, 'echam'); title(sprintf('%s', upper(type))); end;
            caxis([-1 1]);
            cb = colorbar('limits', [-1 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (unitless)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_r1z_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % DEVIATION R1z lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            figure(); clf; hold all; box on;
            cmp = colCog(10);
            colormap(flipud(cmp));
            dev = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
            contourf(mesh_lat, mesh_mon, dev, 1/2*[-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'w');
            contour(mesh_lat, mesh_mon,  dev, [0 0], 'color', 0.75*[1 1 1]);
            contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ep*[1 1], 'linecolor', par.orange, 'linewidth', 1.5);
            contour(mesh_lat, mesh_mon,  flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), par.ga*[1 1], 'linecolor', par.cyan, 'linewidth', 1.5);
            [C, h] = contour(mesh_lat, mesh_mon, dev, [par.ep par.ga], 'linecolor', 'w', 'linewidth', 0.25);
            % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 'k');
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s', upper(type)));
            elseif strcmp(type, 'merra2'); title(sprintf('%s', upper(type)));
            elseif strcmp(type, 'gcm'); title(sprintf('%s', par.model));
            elseif strcmp(type, 'echam'); title(sprintf('%s', upper(type))); end;
            caxis([-0.5 0.5]);
            cb = colorbar('limits', [-0.5 0.5], 'ytick', [-0.5:0.1:0.5], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('%s (unitless)', var_text));
            ylabel('Latitude (deg)');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos);
            set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/flux/%s/%s/0_dr1z_mon_lat', par.plotdir, fw, land), '-dpng', '-r300');
            close;

            % % GA MALR and GA_BL DALR contour lat x mon dependence of RCE and RAE
            % figure(); clf; hold all;
            % cmp = colCog(10);
            % colormap(flipud(cmp));
            % contourf(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [-16 -8 -4 -2 -1:0.2:1 2 4 8 16], 'linecolor', 'w');
            % contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [0 0], 'color', 0.75*[1 1 1]);
            % [C, h] = contour(mesh_lat, mesh_mon, flux_z.(land).r1.(fw), [par.ep par.ga], 'linecolor', 0.5*[1 1 1], 'linewidth', 1);
            % clabel(C, h, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.5*[1 1 1]);
            % contour(mesh_lat, mesh_mon, ga_malr_diff, par.ga_thresh*[1 1], 'color', par.orange, 'linewidth', 2);
            % contour(mesh_lat, mesh_mon, ga_dalr_bl_diff, par.ga_bl_thresh*[1 1], 'color', par.blue, 'linewidth', 2);
            % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
            % elseif strcmp(type, 'merra2'); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh));
            % elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', par.model, par.ga_thresh, par.ga_bl_thresh));
            % elseif any(strcmp(type, {'echam'})); title(sprintf('%s, $\\Gamma_m - \\Gamma / \\Gamma_m < %g \\%%$, $(\\Gamma_d - \\Gamma )/ \\Gamma_d > %g \\%%$', upper(type), par.ga_thresh, par.ga_bl_thresh)); end;
            % caxis([-1 1]);
            % cb = colorbar('limits', [-0.4 1], 'ytick', [-1:0.2:1], 'location', 'eastoutside');
            % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
            % ylabel(cb, sprintf('%s (unitless)', var_text));
            % xlabel('Month'); ylabel('Latitude (deg)');
            % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabel, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            % print(sprintf('%s/flux/%s/%s/0_r1_mon_lat_ga_malr_ga_dalr_bl_overlay', par.plotdir, fw, land), '-dpng', '-r300');
            % close;

            [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim3.lon, lat);
            % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            for t = {'ann'}; time = t{1};
                % lat x lon of RCE and RAE
                if strcmp(fw, 'dse'); var_text = '$\nabla \cdot F_s$';
                else var_text = '$\nabla \cdot F_m$'; end
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).res.(fw)', 'linecolor', 'none');
                caxis([-150 150]);
                if strcmp(type, 'gcm')
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
                else
                    contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                end
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/div_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                var_text = '$R_1$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(flipud(cmp));
                contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).r1.(fw)', [-16 -8 -4 -2 -1:0.1:1 2 4 8 16], 'linecolor', 'none');
                if strcmp(type, 'gcm')
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
                else
                    contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                end
                caxis([-1 1]);
                cb = colorbar('limits', [0 1], 'ytick', [-1:0.1:1], 'location', 'eastoutside');
                cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s',upper(type), var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/r1_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

                var_text = '$R_a$';
                figure(); clf; hold all;
                cmp = colCog(20);
                colormap(cmp);
                contourf(mesh_ll_lat, mesh_ll_lon, flux_t.(land).(time).ra.(fw)', 'linecolor', 'none');
                if strcmp(type, 'gcm')
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
                    contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
                else
                    contour(par.landlon+180, par.landlat, par.land, [1 1], 'k');
                end
                caxis([-150 150]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, %s', upper(type), var_text, upper(time), land_text));
                elseif any(strcmp(type, 'gcm')); title(sprintf('%s, %s, %s, %s', par.model, var_text, upper(time), land_text));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, %s, %s',upper(type), var_text, upper(time), land_text)); end;
                xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/flux/%s/%s/%s/ra_lat_lon', par.plotdir, fw, land, time), '-dpng', '-r300');
                close;

            end % for time
        end % for mse dse
    end % for land
end
function plot_temp(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.model, par.gcm.clim, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    elseif strcmp(type, 'echam')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ta_si.mat', type, par.echam.clim, par.lat_interp, par.ep, par.ga));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g_ga_%g/ma_si.mat', type, par.echam.clim, par.lat_interp, par.ep, par.ga));
        plev = grid.dim3.plev/100; % convert Pa to hPa
    end

    make_dirs_ep(type, par)

    for f = {'mse'}; fw = f{1};
        for c = fieldnames(ta_si.rce.tp.(fw))'; crit = c{1};
            % for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            for t = {'ann'}; time = t{1};
                % for l = {'lo', 'l', 'o'}; land = l{1};
                for l = {'lo'}; land = l{1};
                    if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                    elseif strcmp(land, 'l'); land_text = 'Land';
                    elseif strcmp(land, 'o'); land_text = 'Ocean';
                    end

                % RCE and RAE separated into NH and SH
                    figure(); clf; hold all;
                    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.maroon);
                    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, '-', 'color', par.orange);
                    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, '--', 'color', par.orange);
                    h_rae_nh = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, '-', 'color', par.blue);
                    h_rae_sh = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, '--', 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('%s, %s', upper(time), land_text));
                    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE', 'NH ML RCE', 'SH ML RCE', 'NH RAE', 'SH RAE', 'location', 'northeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    line([0 0], [100 1000], 'linewidth', 0.5, 'color', 'k');
                    h_rce_all = plot(ta_si.rce.all.(fw).(crit).(land).(time)-ma_si.rce.all.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[200 1000], 'xtick', [-5:5:40], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_all_diff', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Difference between RCE temperature profile and moist adiabat
                    figure(); clf; hold all;
                    line([0 0], [100 1000], 'linewidth', 0.5, 'color', 'k');
                    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time)-ma_si.rce.tp.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.maroon);
                    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time)-ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, 'color', par.orange);
                    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time)-ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, '--', 'color', par.orange);
                    xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_tp, h_rce_nh, h_rce_sh], 'Tropics', 'NH ML', 'SH ML', 'location', 'southeast');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[200 1000], 'xtick', [-5:5:40], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_diff', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % All RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.all.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    % h_rae_ma_si = plot(ma_si.rae.all.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RAE, %s, %s', upper(time), land_text));
                    % if strcmp(type, 'era5') | strcmp(type, 'erai')
                    %     legend([h_rae, h_rae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    % elseif strcmp(type, 'gcm')
                    %     legend([h_rae, h_rae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    % end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % All RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.all.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.all.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_all', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % Tropical RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.tp.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.maroon);
                    h_rce_ma_si = plot(ma_si.rce.tp.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.maroon);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_tp', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RCE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RCE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
                        legend([h_rce, h_rce_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rce, h_rce_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RCAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rcae = plot(ta_si.rcae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RCAE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
                        legend([h_rcae, h_rcae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rcae, h_rcae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RCAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RCAE, %s, %s', upper(time), land_text));
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'merra2')
                        legend([h_rcae, h_rcae_ma_si], upper(type), 'Moist adiabat', 'location', 'southwest');
                    elseif strcmp(type, 'gcm')
                        legend([h_rcae, h_rcae_ma_si], par.model, 'Moist adiabat', 'location', 'southwest');
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rcae_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % NH RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH RAE, %s, %s', upper(time), land_text));
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % SH RAE compared with moist adiabat
                    figure(); clf; hold all;
                    h_rae = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH RAE, %s, %s', upper(time), land_text));
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rae_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % ALL NH compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.nh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    h_rae = plot(ta_si.rae.nh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('NH, %s', upper(time)));
                    legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                % ALL SH compared with moist adiabat
                    figure(); clf; hold all;
                    h_rce = plot(ta_si.rce.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.orange);
                    h_rce_ma_si = plot(ma_si.rce.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    h_rcae = plot(ta_si.rcae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', 0.25*[1 1 1]);
                    h_rcae_ma_si = plot(ma_si.rcae.sh.(fw).(crit).(land).(time).ta, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
                    h_rae = plot(ta_si.rae.sh.(fw).(crit).(land).(time), grid.dim3.si, 'color', par.blue);
                    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    title(sprintf('SH, %s', upper(time)));
                    legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xlim', [200 inf], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[50 1000], 'xminortick', 'on')
                    hline(0, '-k');
                    print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
                    close;

                end % land
                % land and ocean comp
                    % figure(); clf; hold all;
                    % h_rce_l = plot(ta_si.rce.nh.(fw).(crit).l.(time)-ma_si.rce.nh.(fw).(crit).l.(time).ta, grid.dim3.si, 'color', par.orange);
                    % h_rce_o = plot(ta_si.rce.nh.(fw).(crit).o.(time)-ma_si.rce.nh.(fw).(crit).o.(time).ta, grid.dim3.si, ':', 'color', par.orange);
                    % xlabel('$T-T_m$ (K)'); ylabel('$\sigma$ (unitless)');
                    % title(sprintf('Tropical RCE, %s, %s', upper(time), land_text));
                    % legend([h_rce_l, h_rce_o], 'Land', 'Ocean', 'location', 'southeast');
                    % axis('tight');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    % set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', 1e-3*[10 20 50 100 200 300 400:200:1000], 'ylim', 1e-3*[100 1000], 'xminortick', 'on')
                    % hline(0, '-k');
                    % print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/rce_lo_diff', par.plotdir, par.ep, par.ga, fw, crit, 'lo', time), '-dpng', '-r300');
                    % close;
            end % time avg
        end % RCE/RAE definition
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
    figure(); clf; hold all;
    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), plev, 'color', par.maroon);
    h_rce_nh = plot(ta_si.rce.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.orange);
    h_rce_sh = plot(ta_si.rce.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.orange);
    h_rae_nh = plot(ta_si.rae.nh.(fw).(crit).(land).(time), plev, '-', 'color', par.blue);
    h_rae_sh = plot(ta_si.rae.sh.(fw).(crit).(land).(time), plev, '--', 'color', par.blue);
    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
    legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE ($\pm 30^\circ$)', 'NH ML RCE ($>+30^\circ$)', 'SH ML RCE ($<-30^\circ$)', 'NH RAE', 'SH RAE', 'location', 'eastoutside');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [50 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp', par.plotdir), '-dpng', '-r300');
    close;
    % Legend for moist adiabat comparisons
    figure(); clf; hold all;
    h_rce_tp = plot(ta_si.rce.tp.(fw).(crit).(land).(time), plev, '-k');
    h_rce_tp_ma_si = plot(ma_si.rce.tp.(fw).(crit).(land).(time).ta, plev, ':k');
    xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        legend([h_rce_tp h_rce_tp_ma_si], upper(type), 'Moist adiabat', 'location', 'eastoutside');
    elseif strcmp(type, 'gcm')
        legend([h_rce_tp h_rce_tp_ma_si], par.model, 'Moist adiabat', 'location', 'eastoutside');
    end
    title(upper(sprintf('%s', time)));
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [50 1000], 'xminortick', 'on')
    hline(0, '-k');
    print(sprintf('%s/legends/temp_ma_si', par.plotdir), '-dpng', '-r300');
    close;
end
function plot_temp_ann(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/ta_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ma_mon_lat.mat', prefix_proc, par.lat_interp));

    make_dirs_ep(type, par)

    time = 'ann';
    crit = 'def';

    for f = {'mse'}; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            % take annual mean
            r1_ann = nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2);
            ta_ann = squeeze(nanmean(tasi.(land), 2));
            ma_ann = squeeze(nanmean(masi.(land), 2));

            % locate rce, rcae, and rae for NH and SH
            idx_rce_nh = find(r1_ann<=par.ep & grid.dim2.lat>0);
            idx_rcae_nh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim2.lat>0);
            idx_rae_nh = find(r1_ann>=par.ga & grid.dim2.lat>0);

            idx_rce_sh = find(r1_ann<=par.ep & grid.dim2.lat<0);
            idx_rcae_sh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim2.lat<0);
            idx_rae_sh = find(r1_ann>=par.ga & grid.dim2.lat<0);

            % take area averaged temperature profile
            ta_rce_nh = nansum( ta_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
            ma_rce_nh = nansum( ma_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
            ta_rcae_nh = nansum( ta_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
            ma_rcae_nh = nansum( ma_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
            ta_rae_nh = nansum( ta_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

            ta_rce_sh = nansum( ta_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ma_rce_sh = nansum( ma_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ta_rcae_sh = nansum( ta_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ma_rcae_sh = nansum( ma_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ta_rae_sh = nansum( ta_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

            % moist adiabats have nans below the initialization level so the nansums are 0 there. Make these spurious zeros nans.
            if ~strcmp(par.ma_init, 'surf')
                ma_rce_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rce_sh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_sh(grid.dim3.si>par.ma_init) = nan;
            end

            % ALL NH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_nh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_nh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_nh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_nh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('NH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh_ann', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % ALL SH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_sh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_sh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_sh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_sh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('SH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh_ann', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % NARROW ALL NH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_nh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_nh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_nh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_nh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('NH, %s', upper(time)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [190 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh_ann_vert', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % NARROW ALL SH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_sh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_sh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_sh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_sh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('SH, %s', upper(time)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [190 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh_ann_vert', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

        end % land
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
end
function plot_dr1_line(type, par)
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

    if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    for f = f_vec; fw = f{1};
        for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
            folder = sprintf('%s/dr1/%s/0_lat_%g', par.plotdir, fw, lat_eval);
            if ~exist(folder, 'dir'); mkdir(folder); end;

            r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat_eval);

            dr1 = flux_z.lo.r1.(fw) - r1_ann;
            dr1_lat = interp1(grid.dim3.lat, dr1, lat_eval);

            r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
            comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);

            r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
            comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
            comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);

            % DELTA R1 lat x mon dependence of RCE and RAE
            var_text = '$\Delta R_1$';
            ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]);
            ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]);
            figure(); clf; hold all; box on;
            if abs(lat_eval)<60
                rcemax = par.ep-r1_ann_lat(1);
                vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
            else
                raemin = par.ga-r1_ann_lat(1);
                vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            end
            line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
            tot=plot([1:12], dr1_lat, 'k');
            c12=plot([1:12], comp1_lat+comp2_lat, '-.', 'color', 0.5*[1 1 1]);
            c1=plot([1:12], comp1_lat, '--', 'color', par.maroon);
            c2=plot([1:12], comp2_lat, ':', 'color', par.blue);
            if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$', upper(type), var_text, lat_eval));
            elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$', par.model, var_text, lat_eval)); end;
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
            for le = 1:length(lat_eval_list); lat_eval = lat_eval_list(le);
                folder = sprintf('%s/dr1/%s/%s/0_lat_%g', par.plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat_eval);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat_eval);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1 = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1_lat = interp1(grid.dim3.lat, comp1, lat_eval);

                comp2_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2 = fm_ann./(ra_ann).^2.*delta_stf;
                comp2_lat = interp1(grid.dim3.lat, comp2, lat_eval);

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo); ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up); ylim_up = 1; end;
                figure(); clf; hold all; box on;
                if abs(lat_eval)<60
                    % rce=plot([1:12], par.ep-r1_ann_lat, 'color', par.orange);
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                else
                    % rae=plot([1:12], par.ga-r1_ann_lat, 'color', par.blue);
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper(type), var_text, land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', par.model, var_text, land_text, lat_eval)); end;
                xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                if abs(lat_eval)<60
                    rcemax = par.ep-r1_ann_lat(1);
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                else
                    raemin = par.ga-r1_ann_lat(1);
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], dr1_lat, 'k');
                c12=plot([1:12], comp1_lat+comp2_lat, '-.k');
                c1=plot([1:12], comp1_lat, '--k');
                c2=plot([1:12], comp2_lat, ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', upper(type), var_text, land_text, lat_eval));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$', par.model, var_text, land_text, lat_eval)); end;
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
function plot_dr1_midlatitude_line(type, par)
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
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [-5 5];

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
                if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
                else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                folder = sprintf('%s/dr1/%s/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;

                % R1 computed before zonal averaging
                r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                r1_lat = nansum(r1_lat.*clat_mon)/nansum(clat);
                comp1r_lat = interp1(grid.dim3.lat, flux_z.(land).comp1.(fw), lat);
                comp1r_lat = nansum(comp1r_lat.*clat_mon)/nansum(clat);
                comp2r_lat = interp1(grid.dim3.lat, flux_z.(land).comp2.(fw), lat);
                comp2r_lat = nansum(comp2r_lat.*clat_mon)/nansum(clat);

                % R1 computed after zonal averaging
                r1z_lat = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                r1z_lat = nansum(r1z_lat.*clat_mon)/nansum(clat);

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
                tot=plot([1:12], circshift(r1_lat,shiftby, 2), 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_r1', folder), '-dpng', '-r300');
                close;

                % R1 computed at each lat x lon RES and RA
                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                r1z_ann = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat = interp1(grid.dim3.lat, r1z_ann, lat);
                r1z_ann_lat = nansum(r1z_ann_lat.*clat_mon)/nansum(clat);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                dr1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - r1z_ann;
                dr1z_lat = interp1(grid.dim3.lat, dr1z, lat);
                dr1z_lat = nansum(dr1z_lat.*clat_mon)/nansum(clat);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1a = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1a_lat = interp1(grid.dim3.lat, comp1a, lat);
                comp1a_lat = nansum(comp1a_lat.*clat_mon)/nansum(clat);

                comp2a_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2a = fm_ann./(ra_ann).^2.*delta_stf;
                comp2a_lat = interp1(grid.dim3.lat, comp2a, lat);
                comp2a_lat = nansum(comp2a_lat.*clat_mon)/nansum(clat);

                % % DELTA R1 lat x mon dependence of RCE and RAE
                % var_text = '$\Delta R_1$';
                % figure(); clf; hold all; box on;
                % colororder({'k', 'k'});
                % yyaxis left
                % ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1_ann_lat(1)+max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                % ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % yyaxis right
                % ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % rcemax = par.ep-r1_ann_lat(1);
                % if rcemax > ylim_lo
                %     vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % raemin = par.ga-r1_ann_lat(1);
                % if raemin < ylim_up
                %     vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                % tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                % elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                % ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                % close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
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
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
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
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
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
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'southoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp', folder), '-dpng', '-r300');
                close;

                % NOLEGEND DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
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
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_noleg', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -0.2;
                ylim_up = 0.5;
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -r1z_ann_lat(1)+ylim_lo;
                ylim_up = -r1z_ann_lat(1)+ylim_up;
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2) - circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside', 'NumColumns', 2);
                legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'southoutside', 'orientation', 'horizontal');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp', folder), '-dpng', '-r300');
                close;

                % NO LEGEND DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -0.3;
                ylim_up = 0.5;
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = -r1z_ann_lat(1)+ylim_lo;
                ylim_up = -r1z_ann_lat(1)+ylim_up;
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2) - circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), -lat_bound+lat_center, lat_bound+lat_center));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), -lat_bound+lat_center, lat_bound+lat_center));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, -lat_bound+lat_center, lat_bound+lat_center)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside');
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_noleg', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    % elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
    %         dlat = 0.25; % step size for standard lat grid
    %         if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
    %         else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dr1/%s/0_midlatitude_pm_lat_%g', par.plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
    %         r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

    %         dr1 = flux_z.lo.r1.(fw) - r1_ann;
    %         dr1_lat = interp1(grid.dim3.lat, dr1, lat);
    %         dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

    %         % DELTA R1 lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_1$';
    %         ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
    %         ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
    %         figure(); clf; hold all; box on;
    %             rcemax = par.ep-r1_ann_lat(1);
    %             vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
    %             patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    %             raemin = par.ga-r1_ann_lat(1);
    %             vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    %             patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat, shiftby, 2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat, shiftby, 2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat, shiftby, 2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound+lat_center, lat_bound+lat_center));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, -lat_bound+lat_center, lat_bound+lat_center));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, -lat_bound+lat_center, lat_bound+lat_center)); end;
    %         legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_1$ (unitless)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
    %         close;
    %     end
    % end

end % for function
function plot_dr1_polar_line(type, par)
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
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE data
    load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    lat_bound_list = [-80 80];

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

                % R1 computed before zonal averaging
                r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                r1_lat = nansum(r1_lat.*clat_mon)/nansum(clat);
                comp1r_lat = interp1(grid.dim3.lat, flux_z.(land).comp1.(fw), lat);
                comp1r_lat = nansum(comp1r_lat.*clat_mon)/nansum(clat);
                comp2r_lat = interp1(grid.dim3.lat, flux_z.(land).comp2.(fw), lat);
                comp2r_lat = nansum(comp2r_lat.*clat_mon)/nansum(clat);

                % R1 computed after zonal averaging
                r1z_lat = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                r1z_lat = nansum(r1z_lat.*clat_mon)/nansum(clat);

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
                tot=plot([1:12], circshift(r1_lat,shiftby, 2), 'k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_r1', folder), '-dpng', '-r300');
                close;

                % R1 computed at each lat x lon RES and RA
                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                r1z_ann = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat = interp1(grid.dim3.lat, r1z_ann, lat);
                r1z_ann_lat = nansum(r1z_ann_lat.*clat_mon)/nansum(clat);

                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                dr1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - r1z_ann;
                dr1z_lat = interp1(grid.dim3.lat, dr1z, lat);
                dr1z_lat = nansum(dr1z_lat.*clat_mon)/nansum(clat);

                delta_fm = flux_z.(land).res.(fw) - fm_ann;
                comp1a = -stf_ann./(ra_ann).^2.*delta_fm;
                comp1a_lat = interp1(grid.dim3.lat, comp1a, lat);
                comp1a_lat = nansum(comp1a_lat.*clat_mon)/nansum(clat);

                comp2a_text = '$\frac{\nabla\cdot F_m}{R_a^2}\Delta (\mathrm{SH+LH})$';
                delta_stf = flux_z.(land).stf.(fw) - stf_ann;
                comp2a = fm_ann./(ra_ann).^2.*delta_stf;
                comp2a_lat = interp1(grid.dim3.lat, comp2a, lat);
                comp2a_lat = nansum(comp2a_lat.*clat_mon)/nansum(clat);

                % % DELTA R1 lat x mon dependence of RCE and RAE
                % var_text = '$\Delta R_1$';
                % figure(); clf; hold all; box on;
                % colororder({'k', 'k'});
                % yyaxis left
                % ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1_ann_lat(1)+max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
                % ylabel(sprintf('$R_1$ (unitless)'));
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % yyaxis right
                % ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                % rcemax = par.ep-r1_ann_lat(1);
                % if rcemax > ylim_lo
                %     vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                %     patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % raemin = par.ga-r1_ann_lat(1);
                % if raemin < ylim_up
                %     vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                %     patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                % end
                % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                % tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
                % if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                % elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                % elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                % ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
                % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                % print(sprintf('%s/0_mon_dr1', folder), '-dpng', '-r300');
                % close;

                % DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1a_lat comp2a_lat comp1a_lat+comp2a_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby, 2), 'k');
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
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
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
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
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

                % NOLEGEND DELTA R1 lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                ylim_lo = r1_ann_lat(1)+min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = r1_ann_lat(1)+max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                tot=plot([1:12], circshift(r1_lat, shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                ylim_lo = min([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                ylim_up = max([dr1_lat comp1r_lat comp2s_lat comp1r_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
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
                tot=plot([1:12], circshift(dr1_lat,shiftby,2), 'k');
                c12=plot([1:12], circshift(comp1r_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1r_lat,shiftby,2), '--k');
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), ':k');
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, land_text, lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, land_text, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1_decomp_noleg', folder), '-dpng', '-r300');
                close;

                % DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = 0.2;
                ylim_up = 2;
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = ylim_lo - r1z_ann_lat(1);
                ylim_up = ylim_up - r1z_ann_lat(1);
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2)-circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                legend([tot res c1 c2], '$\Delta R_1$', 'Residual', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside', 'numcolumns', 2);
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp', folder), '-dpng', '-r300');
                close;

                % NOLEG DELTA R1Z lat x mon dependence of RCE and RAE
                var_text = '$\Delta R_1$';
                figure(); clf; hold all; box on;
                colororder({'k', 'k'});
                yyaxis left
                % ylim_lo = r1z_ann_lat(1)+min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = r1z_ann_lat(1)+max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = 0.2;
                ylim_up = 2;
                tot=plot([1:12], circshift(r1z_lat,shiftby,2), 'k');
                ylabel(sprintf('$R_1$ (unitless)'));
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                yyaxis right
                % ylim_lo = min([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_lo) | ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max([dr1z_lat comp1s_lat comp2s_lat comp1s_lat+comp2s_lat]); if isnan(ylim_up) | ylim_up==0; ylim_up = 1; end;
                ylim_lo = ylim_lo - r1z_ann_lat(1);
                ylim_up = ylim_up - r1z_ann_lat(1);
                rcemax = par.ep-r1z_ann_lat(1);
                if rcemax > ylim_lo
                    vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                raemin = par.ga-r1z_ann_lat(1);
                if raemin < ylim_up
                    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                end
                line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
                tot=plot([1:12], circshift(dr1z_lat,shiftby,2), 'k');
                % c12=plot([1:12], circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                res=plot([1:12], circshift(dr1z_lat,shiftby,2)-circshift(comp1s_lat+comp2s_lat,shiftby,2), '-.k');
                c1=plot([1:12],  circshift(comp1s_lat,shiftby,2), '-', 'color', par.maroon);
                c2=plot([1:12],  circshift(comp2s_lat,shiftby,2), '-', 'color', 0.5*[1 1 1]);
                if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_bound, lat_pole));
                elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_bound, lat_pole));
                elseif strcmp(type, 'gcm'); title(sprintf('%s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, lat_bound, lat_pole)); end;
                % xlabel('Month');
                ylabel(sprintf('$\\Delta R_1$ (unitless)'));
                % legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1\mathrm{\,linear}}$', '$\Delta (\nabla\cdot F_m)$', '$\Delta R_a$', 'location', 'eastoutside');
                % legend([tot res c1 c2], '$\Delta R_1$', '$\Delta R_{1}-\left(\frac{\Delta\left(\nabla\cdot F_m\right)}{\overline{R_a}}-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a\right)$', '$\frac{\Delta (\nabla\cdot F_m)}{\overline{R_a}}$', '$-\frac{\overline{\nabla\cdot F_m}}{\overline{R_a^2}}\Delta R_a$', 'location', 'eastoutside');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/0_mon_dr1z_decomp_noleg', folder), '-dpng', '-r300');
                close;

            end

        end % for mse dse
    end % for land

    % if any(strcmp(type, {'era5', 'erai'})); f_vec = par.era.fw;
    % elseif any(strcmp(type, 'merra2')); f_vec = par.merra2.fw;
    % elseif any(strcmp(type, 'echam')); f_vec = par.echam.fw;
    % elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw; end
    % for f = f_vec; fw = f{1};
    %     for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
    %         dlat = 0.25; % step size for standard lat grid
    %         if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; monlabel=par.monlabel; shiftby=0;
    %         else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; monlabel=par.monlabelsh; shiftby=6; end;
    %         clat = cosd(lat); % cosine of latitude for cosine weighting
    %         clat_mon = repmat(clat', [1 12]);

    %         folder = sprintf('%s/dr1/%s/0_poleward_of_lat_%g', par.plotdir, fw, lat_bound);
    %         if ~exist(folder, 'dir'); mkdir(folder); end;

    %         r1_ann = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
    %         r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

    %         dr1 = flux_z.lo.r1.(fw) - r1_ann;
    %         dr1_lat = interp1(grid.dim3.lat, dr1, lat);
    %         dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_l = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp1 = sftlf*1e-2.*(flux_z.l.r1.(fw) - r1_ann_l);
    %         comp1_lat = interp1(grid.dim3.lat, comp1, lat);
    %         comp1_lat = nansum(comp1_lat.*clat_mon)/nansum(clat);

    %         r1_ann_o = repmat(nanmean(flux_z.lo.r1.(fw), 2), [1 12]);
    %         comp2 = (1-sftlf*1e-2).*(flux_z.o.r1.(fw) - r1_ann_o);
    %         comp2_lat = interp1(grid.dim3.lat, comp2, lat);
    %         comp2_lat = nansum(comp2_lat.*clat_mon)/nansum(clat);

    %         % DELTA R1 lat x mon dependence of RCE and RAE
    %         var_text = '$\Delta R_1$';
    %         ylim_lo = min([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
    %         ylim_up = max([dr1_lat comp1_lat comp2_lat comp1_lat+comp2_lat]); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
    %         figure(); clf; hold all; box on;
    %             rcemax = par.ep-r1_ann_lat(1);
    %             vertices = [1 ylim_lo; 12 ylim_lo; 12 rcemax; 1 rcemax];
    %             patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    %             raemin = par.ga-r1_ann_lat(1);
    %             vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    %             patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    %         line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    %         tot=plot([1:12], circshift(dr1_lat, shiftby, 2), 'k');
    %         c12=plot([1:12], circshift(comp1_lat+comp2_lat, shiftby, 2), '-.', 'color', 0.5*[1 1 1]);
    %         c1=plot([1:12],  circshift(comp1_lat, shiftby, 2), '--', 'color', par.maroon);
    %         c2=plot([1:12],  circshift(comp2_lat, shiftby, 2), ':', 'color', par.blue);
    %         if any(strcmp(type, {'era5', 'erai'})); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif any(strcmp(type, 'merra2')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), var_text, lat_bound, lat_pole));
    %         elseif strcmp(type, 'gcm'); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $$%g^\\circ$', par.model, var_text, lat_bound, lat_pole)); end;
    %         legend([tot c12 c1 c2], '$\Delta R_1$', '$\Delta R_{1,\mathrm{\,L+O}}$', '$\Delta R_{1,\mathrm{\,L}}$', '$\Delta R_{1,\mathrm{\,O}}$', 'location', 'eastoutside');
    %         xlabel('Month');
    %         ylabel(sprintf('$\\Delta R_1$ (unitless)'));
    %         set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    %         print(sprintf('%s/0_mon_dr1_lo_decomp', folder), '-dpng', '-r300');
    %         close;
    %     end
    % end

end % for function

% HELPER functions
function [flux_zt, vh, vh_mon, lat, par] = load_flux(type, par)
    % load processed data/proc
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_zt.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh.mat', type, par.lat_interp));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/vh_mon.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'merra2')
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
    elseif strcmp(type, 'merra2')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    elseif strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
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
        for m = [1 6 7]; month = m(1);
            if ~exist(sprintf('%s/temp_zon_sel/%s/%g', par.plotdir, land, month), 'dir')
                mkdir(sprintf('%s/temp_zon_sel/%s/%g', par.plotdir, land, month));
            end
            if ~exist(sprintf('%s/thetaeq_zon_sel/%s/%g', par.plotdir, land, month), 'dir')
                mkdir(sprintf('%s/thetaeq_zon_sel/%s/%g', par.plotdir, land, month));
            end
            if ~exist(sprintf('%s/dtdz_zon_sel/%s/%g', par.plotdir, land, month), 'dir')
                mkdir(sprintf('%s/dtdz_zon_sel/%s/%g', par.plotdir, land, month));
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
            if par.do_surf; v_vec = {'p', 'z', 'si', 'pi'};
            else v_vec = {'p', 'z', 'si'}; end
            for v = v_vec; vert = v{1};
                if ~exist(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zon/%s/%s/%s', par.plotdir, land, time, vert));
                end
                if ~exist(sprintf('%s/temp_zonmean/%s/%s/%s', par.plotdir, land, time, vert), 'dir')
                    mkdir(sprintf('%s/temp_zonmean/%s/%s/%s', par.plotdir, land, time, vert));
                end
            end
            if any(strcmp(type, {'era5', 'erai'})); fw_vec = par.era.fw;
            elseif strcmp(type, 'merra2'); fw_vec = par.merra2.fw;
            elseif strcmp(type, 'gcm'); fw_vec = par.gcm.fw;
            elseif contains(type, 'echam'); fw_vec = par.echam.fw; end;
            for f = fw_vec; fw = f{1};
                if ~exist(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/flux/%s/%s/%s', par.plotdir, fw, land, time));
                end
                if ~exist(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time), 'dir')
                    mkdir(sprintf('%s/dr2/%s/%s/%s', par.plotdir, fw, land, time));
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
    if ~exist(sprintf('%s/siced', par.plotdir), 'dir')
        mkdir(sprintf('%s/siced', par.plotdir));
    end
    if ~exist(sprintf('%s/sftlf', par.plotdir), 'dir')
        mkdir(sprintf('%s/sftlf', par.plotdir));
    end
    if ~exist(sprintf('%s/legends', par.plotdir), 'dir')
        mkdir(sprintf('%s/legends', par.plotdir));
    end
end
function make_dirs_si_bl(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'merra2')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    elseif strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl')
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s', par.plotdir, par.si_bl, land, time), 'dir')
                mkdir(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s', par.plotdir, par.si_bl, land, time));
            end
        end
    end
end
function make_dirs_ep(type, par)
    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        if any(strcmp(type, {'erai'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'era5'})); f_vec = par.era.fw;
        elseif any(strcmp(type, {'merra2'})); f_vec = par.merra2.fw;
        elseif strcmp(type, 'gcm'); f_vec = par.gcm.fw;
        elseif strcmp(type, 'echam'); f_vec = par.echam.fw;
        elseif strcmp(type, 'echam_ml') | strcmp(type, 'echam_pl'); f_vec = par.echam.fw; end
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

