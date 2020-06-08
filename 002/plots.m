clc; close all; clear variables;

% load useful MATLAB tools
addpath(genpath('/project2/tas1/miyawaki/matlab'));

% lat grid type
lat_interp = 'std'; % don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
% load processed data
load(['/project2/tas1/miyawaki/projects/002/data/processed_data_' lat_interp '.mat']);

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

% create a mon x lat meshgrid for contour plots
[par.mesh_lat, par.mesh_mon] = meshgrid(1:12, lat);

% set how to close energy budget
% if == TETEN, use TETEN data from Donohoe to close energy budget
% if == tf, use SH and LH data from ERA-Interim to close energy budget
subdir = 'tf';
% make figure directory if it does not exist
plotdir = ['./figures_' lat_interp '/'];
if ~exist([plotdir subdir], 'dir')
    mkdir([plotdir subdir]);
end

if strcmp(subdir, 'teten')

    TETEN = don.TETEN; % use Donohoe TETEN
    tf = don.tf_res; % tf is inferred as residual
    r1 = don.r1_d; % corresponding non dimensional number r1
    rcae = don.rcae_d;
    %plots_mon_lat(subdir, plotdir, lat, don, TETEN, tf, r1, rcae, par);
    plots_lat(subdir, plotdir, lat, don, TETEN, tf, par);

elseif strcmp(subdir, 'tf')

    tf = don.tf; % use ERA-Interim turbulent fluxes
    TETEN = don.TETEN_res; % TETEN is inferred as residual
    r1 = don.r1_e; % corresponding non dimensional number r1
    rcae = don.rcae_e;
    %plots_mon_lat(subdir, plotdir, lat, don, TETEN, tf, r1, rcae, par);
    plots_lat(subdir, plotdir, lat, don, TETEN, tf, par);

end

function plots_lat(subdir, plotdir, lat, don, TETEN, tf, par)
% latitude vs energy flux line plots, comparable to Hartmann (2016)
    figure(); clf; hold all;
    plot(lat, nanmean(nanmean(don.ra,1),3), 'color', par.gray)
    plot(lat, nanmean(nanmean(don.TEDIV,1),3), 'color', par.maroon)
    plot(lat, nanmean(nanmean(TETEN,1),3), 'color', par.green)
    plot(lat, -nanmean(nanmean(don.slhf,1),3), 'color', par.blue)
    plot(lat, -nanmean(nanmean(don.sshf,1),3), 'color', par.orange)
    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
    axis('tight');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
    set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'ylim', [-150 150], 'xminortick', 'on', 'yminortick', 'on')
    hline(0, '-k');
    print([plotdir subdir '/era-fig-6-1-hartmann'], '-dpng', '-r300');
    close;
end

function plots_mon_lat(subdir, plotdir, lat, don, TETEN, tf, r1, rcae, par)
    % seasonality vs latitude plots for each term in the energy budget

    % atmospheric radiative cooling
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(cmp);
    contourf(par.mesh_lat, par.mesh_mon, nanmean(don.ra,3)', -200:20:200, 'linecolor', 'none');
    xlabel('Month'); ylabel('Latitude (deg)');
    caxis([-200 200]);
    cb = colorbar('limits', [-180 0], 'ticks', [-200:20:0]);
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$R_a$ (Wm$^{-2}$)', 'fontsize', par.fs);
    set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/ra_mon_lat'], '-dpng', '-r300');
    close;

    % surface turbulent fluxes (SH + LH)
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(cmp);
    contourf(par.mesh_lat, par.mesh_mon, nanmean(tf,3)', -200:20:200, 'linecolor', 'none');
    xlabel('Month'); ylabel('Latitude (deg)');
    caxis([-200 200]);
    cb = colorbar('limits', [-60 160], 'ticks', [-60:20:160]);
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, 'SH + LH (Wm$^{-2}$)', 'fontsize', par.fs);
    set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/tf_mon_lat'], '-dpng', '-r300');
    close;

    % atmospheric moist static energy storage (TETEN)
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(cmp);
    contourf(par.mesh_lat, par.mesh_mon, nanmean(TETEN,3)', -200:20:200, 'linecolor', 'none');
    xlabel('Month'); ylabel('Latitude (deg)');
    caxis([-200 200]);
    cb = colorbar('limits', [-60 160], 'ticks', [-60:20:160]);
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\frac{\partial h}{\partial t}$ (Wm$^{-2}$)', 'fontsize', par.fs);
    set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/teten_mon_lat'], '-dpng', '-r300');
    close;

    % atmospheric moist static energy divergence (TEDIV)
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(cmp);
    contourf(par.mesh_lat, par.mesh_mon, nanmean(don.TEDIV,3)', -200:20:200, 'linecolor', 'none');
    xlabel('Month'); ylabel('Latitude (deg)');
    caxis([-200 200]);
    cb = colorbar('limits', [-120 120], 'ticks', [-120:20:120]);
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$\nabla\cdot(\vec{v}h)$ (Wm$^{-2}$)', 'fontsize', par.fs);
    set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/tediv_mon_lat'], '-dpng', '-r300');
    close;

    % non-dimensional number R1
    figure(); clf; hold all;
    cmp = colCog(20);
    colormap(cmp);
    contourf(par.mesh_lat, par.mesh_mon, nanmean(r1,3)', -1:0.1:1, 'linecolor', 'none');
    xlabel('Month'); ylabel('Latitude (deg)');
    caxis([-1 1]);
    cb = colorbar('limits', [-1 1], 'ticks', [-1:0.2:1]);
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, '$R_1$ (unitless)', 'fontsize', par.fs);
    set(gca, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print([plotdir subdir '/r1_mon_lat'], '-dpng', '-r300');
    close;

    % spatio-temporal dependence of RCE and RAE
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
