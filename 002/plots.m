clc; close all; clear variables;

addpath(genpath('/project2/tas1/miyawaki/matlab'));

%% set parameters
% lat grid type
par.lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
par.ep_swp = [0.1:0.05:0.35]; % threshold value for determining RCE and RAE
% set how to close energy budget
% if == teten, use TETEN data from Donohoe to close energy budget
% if == stf, use SH and LH data from ERA-Interim to close energy budget
subdir = 'teten';

% load grid
load('/project2/tas1/miyawaki/projects/002/data/read/era_grid.mat');
% load processed data/proc
load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/fluxes.mat', par.lat_interp));
% create a mon x lat meshgrid for contour plots
[par.mesh_lat, par.mesh_mon] = meshgrid(1:12, lat);
% make figure directory if it does not exist
plotdir = sprintf('./figures/%s/', par.lat_interp);
if ~exist([plotdir subdir], 'dir')
    mkdir([plotdir subdir]);
end
% define TETEN and stf depending on closure method
if strcmp(subdir, 'teten')
    TETEN = donz.TETEN; % use Donohoe TETEN
    stf = donz.stf_res; % stf is inferred as residual
    r1 = donz.r1_d; % corresponding non dimensional number r1
elseif strcmp(subdir, 'stf')
    stf = donz.stf; % use ERA-Interim turbulent fluxes
    TETEN = donz.TETEN_res; % TETEN is inferred as residual
    r1 = donz.r1_e; % corresponding non dimensional number r1
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
% plot_ra_tediv_mon_lat(subdir, plotdir, lat, donz, par); % plot R_a and flux divergence profiles
plot_energy_lat(subdir, plotdir, lat, donz, TETEN, stf, par); % plot all energy fluxes vs latitude a la Fig. 6.1 in Hartmann (2016)
% plot_teten_stf_r1_mon_lat(subdir, plotdir, lat, TETEN, stf, r1, par); % plot remaining energy fluxes that depends on energy closure method

% sweep through various threshold values
% for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
%     if ~exist([plotdir subdir sprintf('/rcae_%g', par.ep)])
%         mkdir([plotdir subdir sprintf('/rcae_%g', par.ep)]);
%     end
%     plot_rcae_mon_lat(subdir, plotdir, lat, par); % plot RCAE regimes, depends on choice of threshold epsilon
%     plot_temp(subdir, plotdir, plev_era, par); % plot temperature profiles in RCAE regimes
% end

%% define functions
function plot_ra_tediv_mon_lat(subdir, plotdir, lat, donz, par)
    % seasonaity vs latitude plots for each term in the energy budget
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
    print([plotdir '/ra_mon_lat'], '-dpng', '-r300');
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
    print([plotdir '/tediv_mon_lat'], '-dpng', '-r300');
    close;

end
function plot_energy_lat(subdir, plotdir, lat, donz, TETEN, stf, par)
% latitude vsenergy flux line plots, comparable to Hartmann (2016)
    figure();clf; hold all;
    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
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
    set(gca, 'fontsize', par.fs, 'xminortick', 'on', 'yminortick', 'on')
    print([plotdir subdir '/era-fig-6-1-hartmann'], '-dpng', '-r300');
    close;
end
function plot_teten_stf_r1_mon_lat(subdir, plotdir, lat, TETEN, stf, r1, par)
    % surface turbulent fluxes (SH + LH)
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

    % non-dimensional number R1
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
function plot_rcae_mon_lat(subdir, plotdir, lat, par)
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/eps_%g/rcae.mat', par.lat_interp, par.ep));
    if strcmp(subdir, 'teten'); rcae_plot = rcae.d; end;
    if strcmp(subdir, 'stf'); rcae_plot = rcae.e; end;
    % spatio-emporal dependence of RCE and RAE
    figure(); clf; hold all;
    cmp = colCog(10);
    colormap(cmp);
    imagesc([1 12], [lat(1) lat(end)], rcae_plot');
    caxis([-2 2]);
    xlabel('Month'); ylabel('Latitude (deg)');
    set(gca, 'xlim', [0.5 12.5], 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir sprintf('/rcae_%g/rcae_mon_lat', par.ep)], '-dpng', '-r300');
    close;
end
function plot_temp(subdir, plotdir, plev_era, par)
% vertical temperature profile
    load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/eps_%g/vert_filt.mat', par.lat_interp, par.ep));
    figure(); clf; hold all;
    h_rce = plot(vert_filt.rce_d, plev_era, 'color', par.orange);
    h_rae = plot(vert_filt.rae_d, plev_era, 'color', par.blue);
    xlabel('T (K)'); ylabel('p (hPa)');
    legend([h_rce h_rae], 'RCE', 'RAE');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'log', 'ytick', [100 200 300 400:200:1000], 'xminortick', 'on')
    hline(0, '-k');
    print([plotdir subdir sprintf('/rcae_%g/temp', par.ep)], '-dpng', '-r300');
    close;
end
