clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.ep = 0.1; % threshold value for determining RCE and RAE
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget
% if == stf, use SH and LH data from ERA-Interim to close energy budget
subdir = 'teten';

% load grid
load('/project2/tas1/miyawaki/projects/002/data/proc/era_grid.mat');
% load processed data/proc
load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/processed_rad_%s.mat', par.ep, par.lat_interp));
% create a mon x lat meshgrid for contour plots
[par.mesh_lat, par.mesh_mon] = meshgrid(1:12, lat);
% make figure directory if it does not exist
plotdir = sprintf('./figures/%s/eps_%g/', par.lat_interp, par.ep);
if ~exist([plotdir subdir], 'dir')
    mkdir([plotdir subdir]);
end
% define TETEN and stf depending on closure method
if strcmp(subdir, 'teten')
    TETEN = donz.TETEN; % use Donohoe TETEN
    stf = donz.stf_res; % stf is inferred as residual
    r1 = donz.r1_d; % corresponding non dimensional number r1
    rcae = donz.rcae_d;
elseif strcmp(subdir, 'stf')
    stf = donz.stf; % use ERA-Interim turbulent fluxes
    TETEN = donz.TETEN_res; % TETEN is inferred as residual
    r1 = donz.r1_e; % corresponding non dimensional number r1
    rcae = donz.rcae_e;
end
% set default figure parameters
if 1
    par.ppos = [0 0 10/3 7/3];
    par.ppos_sq = [0 0 10/3 10/3];
    par.ppos_wide = [0 0 13/3 7/3];
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
plot_rad_lat(subdir, plotdir, lat, donz, TETEN, stf, par);
plot_rad_mon_lat(subdir, plotdir, lat, donz, TETEN, stf, r1, rcae, par);
plot_temp(subdir, plotdir, plev_era, par)

%% define functions
function plot_rad_lat(subdir, plotdir, lat, donz, TETEN, stf, par)
% latitude vsenergy flux line plots, comparable to Hartmann (2016)
    figure();clf; hold all;
    plot(lat,nanmean(donz.ra,1), 'color', par.gray)
    plot(lat,nanmean(donz.TEDIV,1), 'color', par.maroon)
    plot(lat,nanmean(TETEN,1), 'color', par.green)
    if strcmp(subdir, 'stf')
        plot(lat, -nanmean(donz.slhf,1), 'color', par.blue)
        plot(lat, -nanmean(donz.sshf,1), 'color', par.orange)
    elseif strcmp(subdir, 'teten')
        plot(lat, nanmean(stf,1), '--', 'color', par.blue)
        plot(lat, nanmean(stf,1), ':', 'color', par.orange)
    end
    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    hline(0, '-k');
    print([plotdir subdir '/era-fig-6-1-hartmann'], '-dpng', '-r300');
    close;
end
function plot_rad_mon_lat(subdir, plotdir, lat, donz, TETEN, stf, r1, rcae, par)
    % seasonaity vs latitude plots for each term in the energy budget

    if 0
        % atmospheric radiative cooling
        figure();clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        contourf(par.mesh_lat, par.mesh_mon, donz.ra', -200:20:200, 'linecolor', 'none');
        xlabel('Month'); ylabel('Latitude (deg)');
        caxis([-200 200]);
        cb = colorbar('limits', [-180 0], 'ticks', [-200:20:0]);
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, '$R_a$ (Wm$^{-2}$)', 'fontsize', par.fs);
        set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print([plotdir subdir '/ra_mon_lat'], '-dpng', '-r300');
        close;

        % surfaceturbulent fluxes (SH + LH)
        figure(); clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        contourf(par.mesh_lat, par.mesh_mon, stf', -200:20:200, 'linecolor', 'none');
        xlabel('Month'); ylabel('Latitude (deg)');
        caxis([-200 200]);
        cb = colorbar('limits', [-60 160], 'ticks', [-60:20:160]);
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, 'SH + LH (Wm$^{-2}$)', 'fontsize', par.fs);
        set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print([plotdir subdir '/stf_mon_lat'], '-dpng', '-r300');
        close;

        % atmosphric moist static energy storage (TETEN)
        figure();clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        contourf(par.mesh_lat, par.mesh_mon, TETEN', -200:20:200, 'linecolor', 'none');
        xlabel('Month'); ylabel('Latitude (deg)');
        caxis([-200 200]);
        cb = colorbar('limits', [-60 160], 'ticks', [-60:20:160]);
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, '$\frac{\partial h}{\partial t}$ (Wm$^{-2}$)', 'fontsize', par.fs);
        set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print([plotdir subdir '/teten_mon_lat'], '-dpng', '-r300');
        close;

        % atmosphric moist static energy divergence (TEDIV)
        figure(); clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        contourf(par.mesh_lat, par.mesh_mon, donz.TEDIV', -200:20:200, 'linecolor', 'none');
        xlabel('Month'); ylabel('Latitude (deg)');
        caxis([-200 200]);
        cb = colorbar('limits', [-120 120], 'ticks', [-120:20:120]);
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, '$\nabla\cdot(\vec{v}h)$ (Wm$^{-2}$)', 'fontsize', par.fs);
        set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print([plotdir subdir '/tediv_mon_lat'], '-dpng', '-r300');
        close;

        % non-dimnsional number R1
        figure(); clf; hold all;
        cmp = colCog(20);
        colormap(cmp);
        contourf(par.mesh_lat, par.mesh_mon, r1', -1:0.1:1, 'linecolor', 'none');
        xlabel('Month'); ylabel('Latitude (deg)');
        caxis([-1 1]);
        cb = colorbar('limits', [-1 1], 'ticks', [-1:0.2:1]);
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, '$R_1$ (unitless)', 'fontsize', par.fs);
        set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        print([plotdir subdir '/r1_mon_lat'], '-dpng', '-r300');
        close;
    end

    % spatio-emporal dependence of RCE and RAE
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    imagesc([1 12], [lat(1) lat(end)], rcae');
    caxis([-2 2]);
    xlabel('Month'); ylabel('Latitude (deg)');
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/rce_rae_mon_lat'], '-dpng', '-r300');
    close;
end
function plot_temp(subdir, plotdir, plev_era, par)
% vertical temperature profile
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/eps_%g/processed_vert_filt_%s.mat', par.ep, par.lat_interp));
    figure(); clf; hold all;
    plot(vert_filt.rce_d, plev_era, 'color', par.orange);
    plot(vert_filt.rae_d, plev_era, 'color', par.blue);
    xlabel('T (K)'); ylabel('p (hPa)');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log')
    hline(0, '-k');
    print([plotdir subdir '/temp'], '-dpng', '-r300');
    close;
end
