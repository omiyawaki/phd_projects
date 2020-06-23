clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.ep_swp = 0.3; % threshold value for determining RCE and RAE
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget (option for ERA-Interim)
% if == stf, use SH and LH data from ERA-Interim to close energy budget (option for ERA-Interim)
% if == era5, use ERA5 radiative cooling and surface turbulent fluxes to close energy budget (only option for ERA5)
par.closure = 'era5';
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
type = 'era5';
% choose_plots(type, par);
for k = 1:length(par.gcm_models); par.model = par.gcm_models{k};
    type = 'gcm';
    % choose_plots(type, par);
end

% sweep through various threshold values
for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
    type = 'era5';
    % choose_plots_ep(type, par)
    for k=1:length(par.gcm_models); par.model = par.gcm_models{k};
        type = 'gcm';
        choose_plots_ep(type, par)
    end
end

%% define functions
function choose_plots(type, par)
    plot_energy_lat(type, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
end
function plot_energy_lat(type, par) % latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, lat, par] = load_flux(type, par); % load data
    figure(); clf; hold all;
    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
    plot(lat,nanmean(flux.ra,2), 'color', par.gray)
    plot(lat,nanmean(flux.res,2), 'color', par.maroon)
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        plot(lat, -nanmean(flux.slhf,2), 'color', par.blue)
        plot(lat, -nanmean(flux.sshf,2), 'color', par.orange)
    elseif strcmp(type, 'gcm')
        plot(lat, nanmean(flux.hfls,2), 'color', par.blue)
        plot(lat, nanmean(flux.hfss,2), 'color', par.orange)
    end
    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'ylim', [-inf 150], 'xminortick', 'on', 'yminortick', 'on')
    print(sprintf('%s/energy-flux', par.plotdir), '-dpng', '-r300');
    close;
% northward MSE transport
    figure(); clf; hold all;
    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
    line([0 0], [min(vh) max(vh)]*10^-15, 'linewidth', 0.5, 'color', 'k');
    plot(lat, vh*10^-15, 'color', par.maroon);
    xlabel('latitude (deg)'); ylabel('PW')
    title('Northward MSE Energy Transport');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
    print([par.plotdir '/mse-transport'], '-dpng', '-r300');
    close;
end

function choose_plots_ep(type, par)
    % plot_rcae_mon_lat(type, par) % plot RCE/RAE regimes
    plot_temp(type, par) % plot temperature profiles
    % plot_ma_diff(type, par) % plot difference of temperature profile from moist adiabat
end
function plot_rcae_mon_lat(type, par)
    % load data
    [~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.model);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.model);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/eps_%g/rcae_z.mat', prefix_proc, par.lat_interp, par.ep)); % load lat x mon RCAE data
    load(sprintf('%s/%s/eps_%g/rcae_t.mat', prefix_proc, par.lat_interp, par.ep)); % load lat x lon RCAE data

    for fn = fieldnames(rcae_z)'; crit = fn{1};
        % lat x mon dependence of RCE and RAE
        figure(); clf; hold all;
        cmp = colCog(10);
        colormap(cmp);
        imagesc([1 12], [lat(1) lat(end)], rcae_z.(crit));
        caxis([-2 2]);
        xlabel('Month'); ylabel('Latitude (deg)');
        set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/eps_%g/%s/rcae_mon_lat', par.plotdir, par.ep, crit), '-dpng', '-r300');
        close;

        for l = {'lo', 'l', 'o'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
                % lat x lon of RCE and RAE
                figure(); clf; hold all;
                cmp = colCog(10);
                colormap(cmp);
                imagesc([grid.dim3.lon(1) grid.dim3.lon(end)], [lat(1) lat(end)], rcae_t.(land).(time).(crit)');
                caxis([-2 2]);
                title(sprintf('%s, %s', upper(time), land_text));
                xlabel('Month'); ylabel('Latitude (deg)');
                set(gca, 'xlim', [0 360], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
                print(sprintf('%s/eps_%g/%s/%s/%s/rcae_lat_lon', par.plotdir, par.ep, crit, land, time), '-dpng', '-r300');
                close;
            end
        end

    end
end
function plot_temp(type, par)
    % load data
    [~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/ta.mat', type, par.lat_interp, par.ep));
        plev = grid.dim3.plev;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ta.mat', type, par.model, par.lat_interp, par.ep));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ma.mat', type, par.model, par.lat_interp, par.ep));
        plev = grid.dim3.plev/100;
    end
    for i = {'def', 'pe', 'cp'}; crit = i{1};
        for j = {'ann', 'djf', 'jja', 'mam', 'son'}; time = j{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                if strcmp(land, 'lo'); land_text = 'Land + Ocean';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

            % RCE and RAE separated into NH and SH
                figure(); clf; hold all;
                h_rce_tp = plot(ta.rce.tp.(crit).(land).(time), plev, 'color', par.maroon);
                h_rce_nh = plot(ta.rce.nh.(crit).(land).(time), plev, '-', 'color', par.orange);
                h_rce_sh = plot(ta.rce.sh.(crit).(land).(time), plev, '--', 'color', par.orange);
                h_rae_nh = plot(ta.rae.nh.(crit).(land).(time), plev, '-', 'color', par.blue);
                h_rae_sh = plot(ta.rae.sh.(crit).(land).(time), plev, '--', 'color', par.blue);
                xlabel('T (K)'); ylabel('p (hPa)');
                title(sprintf('%s, %s', upper(time), land_text));
                legend([h_rce_tp h_rce_nh h_rce_sh h_rae_nh h_rae_sh], 'Tropical RCE', 'NH ML RCE', 'SH ML RCE', 'NH RAE', 'SH RAE', 'location', 'northeast');
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [10 20 50 100 200 300 400:200:1000], 'ylim', [10 1000], 'xminortick', 'on')
                hline(0, '-k');
                print(sprintf('%s/eps_%g/%s/%s/%s/temp/rcae_all', par.plotdir, par.ep, crit, land, time), '-dpng', '-r300');
                close;
            % Tropical RCE compared with moist adiabat
                figure(); clf; hold all;
                h_rce = plot(ta.rce.tp.(crit).(land).(time), plev, 'color', par.maroon);
                h_rce_ma = plot(ma.rce.tp.(crit).(land).(time).ta, plev, ':', 'color', par.maroon);
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
                print(sprintf('%s/eps_%g/%s/%s/%s/temp/rce_tp', par.plotdir, par.ep, crit, land, time), '-dpng', '-r300');
                close;
            % NH RCE compared with moist adiabat
                figure(); clf; hold all;
                h_rce = plot(ta.rce.nh.(crit).(land).(time), plev, 'color', par.orange);
                h_rce_ma = plot(ma.rce.nh.(crit).(land).(time).ta, plev, ':', 'color', par.orange);
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
                print(sprintf('%s/eps_%g/%s/%s/%s/temp/rce_nh', par.plotdir, par.ep, crit, land, time), '-dpng', '-r300');
                close;
            % SH RCE compared with moist adiabat
                figure(); clf; hold all;
                h_rce = plot(ta.rce.sh.(crit).(land).(time), plev, 'color', par.orange);
                h_rce_ma = plot(ma.rce.sh.(crit).(land).(time).ta, plev, ':', 'color', par.orange);
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
                print(sprintf('%s/eps_%g/%s/%s/%s/temp/rce_sh', par.plotdir, par.ep, crit, land, time), '-dpng', '-r300');
                close;
            end
        end
    end
    % Legend for RCE and RAE separated into NH and SH
    figure(); clf; hold all;
    h_rce_tp = plot(ta.rce.tp.def.lo.ann, plev, 'color', par.maroon);
    h_rce_nh = plot(ta.rce.nh.def.lo.ann, plev, '-', 'color', par.orange);
    h_rce_sh = plot(ta.rce.sh.def.lo.ann, plev, '--', 'color', par.orange);
    h_rae_nh = plot(ta.rae.nh.def.lo.ann, plev, '-', 'color', par.blue);
    h_rae_sh = plot(ta.rae.sh.def.lo.ann, plev, '--', 'color', par.blue);
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
    h_rce_tp = plot(ta.rce.tp.(crit).(land).(time), plev, '-k');
    h_rce_tp_ma = plot(ma.rce.tp.(crit).(land).(time).ta, plev, ':k');
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
function plot_ma_diff(type, par)
    % load data
    [~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/eps_%g/ta_mon_lat.mat', type, par.lat_interp, par.ep));
        plev = grid.dim3.plev;
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.model));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ta_mon_lat.mat', type, par.model, par.lat_interp, par.ep));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/eps_%g/ma_mon_lat.mat', type, par.model, par.lat_interp, par.ep));
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
            xlabel('Month'); ylabel('Latitude (deg)');
            if strcmp(type, 'era5') | strcmp(type, 'erai')
                title(sprintf('%s, %s', upper(type), land_text));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, land_text));
            end
            cb = colorbar('limits', [-20 20], 'ytick', [-20:5:20], 'location', 'eastoutside');
            cb.TickLabelInterpreter = 'latex';
            cb.Label.Interpreter = 'latex';
            ylabel(cb, sprintf('$(T - T_m)_{%g \\,\\mathrm{hPa}}$ (K)', plev_eval));
            set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
            print(sprintf('%s/ma_diff/plev_%g/%s/ma_diff_lat_lon', par.plotdir, plev_eval, land), '-dpng', '-r300');
            close;
        end
    end
end

function [flux_z, vh, lat, par] = load_flux(type, par)
    % load processed data/proc
    if strcmp(type, 'era5') | strcmp(type, 'erai')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/flux_z.mat', type, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s', type, par.lat_interp);
    elseif strcmp(type, 'gcm')
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/flux_z.mat', type, par.model, par.lat_interp));
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.model, par.lat_interp);
    end

    % create a mon x lat meshgrid for contour plots
    [par.mesh_lat, par.mesh_mon] = meshgrid(1:12, lat);

    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        for j = {'def', 'pe', 'cp'}; crit = j{1};
            for l = {'lo', 'l', 'o'}; land = l{1};
                for k = {'ann', 'djf', 'jja', 'mam', 'son'}; time = k{1};
                    if ~exist(sprintf('%s/eps_%g/%s/%s/%s/temp', par.plotdir, par.ep, crit, land, time), 'dir')
                        mkdir(sprintf('%s/eps_%g/%s/%s/%s/temp', par.plotdir, par.ep, crit, land, time));
                    end
                end
            end
        end
    end
    for l = {'lo', 'l', 'o'}; land = l{1};
        for plev_eval = [300:100:500]
            if ~exist(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land), 'dir')
                mkdir(sprintf('%s/ma_diff/plev_%g/%s', par.plotdir, plev_eval, land));
            end
        end
    end

    if ~exist(sprintf('%s/legends', par.plotdir), 'dir')
        mkdir(sprintf('%s/legends', par.plotdir));
    end
end
